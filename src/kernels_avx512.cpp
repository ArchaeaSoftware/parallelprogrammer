// Compiled with AVX-512 flags. Nothing here may run before the CPU check in
// kernels_dispatch.cpp -- no static initializers, no unconditional entry.
#include <immintrin.h>

#include <limits>

#include "cbfp/limb_column.hpp"
#include "kernels.hpp"

namespace cbfp {
namespace kernels {
namespace {

constexpr std::uint64_t kFracMask = (std::uint64_t{1} << 52) - 1;

// Eight doubles from a contiguous column. `k` masks off the tail so a partial
// block never reads past the input. There is deliberately no gather here: the
// caller stages a strided column instead, because gathering costs more than
// the decomposition it would feed.
inline __m512d load_column(const double* values, std::size_t row, __mmask8 k)
{
    return _mm512_maskz_loadu_pd(k, values + row);
}

struct Split8 {
    __m512i mantissa;  // odd
    __m512i exponent;  // true-ulp exponent
    __m512i top;       // one past the highest bit position reached
    __mmask8 negative;
    __mmask8 live;  // nonzero mantissa
    __mmask8 nonfinite;
};

// Vector form of the scalar `split`: the same IEEE-754 field extraction.
inline Split8 split8(__m512d v)
{
    const __m512i bits = _mm512_castpd_si512(v);
    const __m512i biased =
        _mm512_and_si512(_mm512_srli_epi64(bits, 52), _mm512_set1_epi64(0x7FF));
    const __m512i frac = _mm512_and_si512(bits, _mm512_set1_epi64(kFracMask));

    const __mmask8 is_normal = _mm512_test_epi64_mask(biased, biased);
    const __m512i m = _mm512_mask_or_epi64(
        frac, is_normal, frac, _mm512_set1_epi64(std::uint64_t{1} << 52));
    const __m512i e =
        _mm512_sub_epi64(_mm512_max_epu64(biased, _mm512_set1_epi64(1)),
                         _mm512_set1_epi64(1075));

    Split8 s;
    s.live = _mm512_test_epi64_mask(m, m);
    s.negative = _mm512_movepi64_mask(bits) & s.live;
    s.nonfinite = _mm512_cmpeq_epi64_mask(biased, _mm512_set1_epi64(0x7FF));
    // top = e + the significand's bit width, which odd-normalization leaves
    // unchanged: it raises e and lowers the width by the same amount. That is
    // why the scan never needs a trailing-zero count.
    s.top = _mm512_add_epi64(
        e, _mm512_sub_epi64(_mm512_set1_epi64(64), _mm512_lzcnt_epi64(m)));

    // tz = popcount((m & -m) - 1); srlv by 64 gives 0, so zero lanes stay zero.
    const __m512i lowbit =
        _mm512_and_si512(m, _mm512_sub_epi64(_mm512_setzero_si512(), m));
    const __m512i tz =
        _mm512_popcnt_epi64(_mm512_sub_epi64(lowbit, _mm512_set1_epi64(1)));
    s.mantissa = _mm512_srlv_epi64(m, tz);
    s.exponent = _mm512_add_epi64(e, tz);
    return s;
}

inline __mmask8 tail_mask(std::size_t row, std::size_t rows)
{
    return row + 8 <= rows ? static_cast<__mmask8>(0xFF)
                           : static_cast<__mmask8>((1u << (rows - row)) - 1u);
}

}  // namespace

Scan scan_column_avx512(const double* values, std::size_t rows,
                        long long floor_exponent)
{
    // First pass touches only the exponent field.
    //
    // The maximum needs nothing more: for a normal, top = e + 53 =
    // biased - 1022, monotonic in that field, and every normal outranks every
    // subnormal (worst normal -1021, best subnormal -1022). With the sign
    // cleared the whole bit pattern orders by magnitude, so one max over it
    // yields both the exponent and, for the all-subnormal case, the
    // significand that decides the width.
    //
    // The minimum is a lower bound only -- the raw exponent, not the true ulp,
    // which sits above it by the significand's trailing zero count. That bound
    // is enough whenever it clears the caller's floor, because then no rescale
    // is due and the exact value would change nothing.
    __m512i vmin = _mm512_set1_epi64(std::numeric_limits<long long>::max());
    __m512i vmax_abs = _mm512_setzero_si512();
    __mmask8 any = 0;
    __mmask8 bad = 0;

    const __m512i kAbs = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFFll);
    const __m512i kExp = _mm512_set1_epi64(0x7FF);
    const __m512i kOne = _mm512_set1_epi64(1);

    for (std::size_t row = 0; row < rows; row += 8) {
        const __mmask8 k = tail_mask(row, rows);
        const __m512i bits = _mm512_castpd_si512(load_column(values, row, k));
        const __m512i abs_bits = _mm512_and_si512(bits, kAbs);
        const __m512i biased = _mm512_srli_epi64(abs_bits, 52);

        const __mmask8 live = _mm512_test_epi64_mask(abs_bits, abs_bits) & k;
        bad |= _mm512_cmpeq_epi64_mask(biased, kExp) & k;
        any |= live;

        vmin = _mm512_mask_min_epu64(vmin, live, vmin,
                                     _mm512_max_epu64(biased, kOne));
        vmax_abs = _mm512_mask_max_epu64(vmax_abs, live, vmax_abs, abs_bits);
    }

    Scan out{0, 0, any != 0, bad != 0};
    if (!out.any || out.nonfinite) return out;

    const long long min_raw =
        static_cast<long long>(_mm512_reduce_min_epu64(vmin)) - 1075;
    const std::uint64_t max_abs =
        static_cast<std::uint64_t>(_mm512_reduce_max_epu64(vmax_abs));
    const std::uint64_t max_biased = max_abs >> 52;
    out.max_top = max_biased != 0 ? static_cast<long long>(max_biased) - 1022
                                  : -1074 + (64 - __builtin_clzll(max_abs));

    if (min_raw >= floor_exponent) {
        out.min_exponent = min_raw;
        return out;
    }

    // Only now, when the exponent could actually drop, does the significand
    // matter -- the true ulp needs its trailing zero count.
    const __m512i kFrac = _mm512_set1_epi64(kFracMask);
    const __m512i kImplicit = _mm512_set1_epi64(std::uint64_t{1} << 52);
    const __m512i k1075 = _mm512_set1_epi64(1075);
    __m512i vulp = _mm512_set1_epi64(std::numeric_limits<long long>::max());

    for (std::size_t row = 0; row < rows; row += 8) {
        const __mmask8 k = tail_mask(row, rows);
        const __m512i bits = _mm512_castpd_si512(load_column(values, row, k));
        const __m512i biased =
            _mm512_and_si512(_mm512_srli_epi64(bits, 52), kExp);
        const __m512i frac = _mm512_and_si512(bits, kFrac);
        const __mmask8 is_normal = _mm512_test_epi64_mask(biased, biased);
        const __m512i m =
            _mm512_mask_or_epi64(frac, is_normal, frac, kImplicit);
        const __mmask8 live = _mm512_test_epi64_mask(m, m) & k;

        // tz = popcount((m & -m) - 1)
        const __m512i e =
            _mm512_sub_epi64(_mm512_max_epu64(biased, kOne), k1075);
        const __m512i lowbit =
            _mm512_and_si512(m, _mm512_sub_epi64(_mm512_setzero_si512(), m));
        const __m512i tz = _mm512_popcnt_epi64(_mm512_sub_epi64(lowbit, kOne));
        vulp = _mm512_mask_min_epi64(vulp, live, vulp, _mm512_add_epi64(e, tz));
    }
    out.min_exponent = _mm512_reduce_min_epi64(vulp);
    return out;
}

// Eight rows per pass. Rows are the vector axis and limb positions the
// sequential one, so the carry chain runs inside a lane and never crosses
// lanes -- no shuffles, no cross-lane propagation.
//
// Decomposition is fused in: mantissa and exponent go straight from the loaded
// doubles into the addend and are never written to memory.
//
// Add and subtract lanes share one loop: subtracting A is adding ~A + 1, so a
// negative lane complements its addend and enters its first limb with carry 1.
// Limbs below a lane's offset stay untouched because the complement is applied
// only where the lane is active.
void accumulate_avx512(std::uint64_t* const* limbs, std::size_t nlimbs,
                       const double* values, std::size_t rows,
                       std::int32_t column_exponent, std::size_t first_limb)
{
    const __m512i kOne = _mm512_set1_epi64(1);
    const __m512i kOnes = _mm512_set1_epi64(-1);
    const __m512i k63 = _mm512_set1_epi64(63);
    const __m512i k64 = _mm512_set1_epi64(64);
    const __m512i kColExp = _mm512_set1_epi64(column_exponent);

    for (std::size_t row = 0; row < rows; row += 8) {
        const __mmask8 k = tail_mask(row, rows);
        const Split8 s = split8(load_column(values, row, k));
        const __mmask8 live = s.live & k;
        if (live == 0) continue;

        const __mmask8 negative = s.negative & live;
        // Dead lanes get shift 0 so their limb offset cannot go negative.
        const __m512i shift = _mm512_maskz_sub_epi64(live, s.exponent, kColExp);

        const __m512i off = _mm512_srli_epi64(shift, 6);
        const __m512i bit = _mm512_and_si512(shift, k63);

        // The 53-bit mantissa lands in at most two limbs. srlv by 64 yields 0,
        // which is exactly what the bit == 0 case wants.
        const __m512i lo = _mm512_maskz_sllv_epi64(live, s.mantissa, bit);
        const __m512i hi = _mm512_maskz_srlv_epi64(live, s.mantissa,
                                                   _mm512_sub_epi64(k64, bit));
        const __m512i off1 = _mm512_add_epi64(off, kOne);

        __m512i carry = _mm512_setzero_si512();
        for (std::size_t p = first_limb; p < nlimbs; ++p) {
            const __m512i pv = _mm512_set1_epi64(static_cast<long long>(p));
            const __mmask8 active = _mm512_cmple_epi64_mask(off, pv) & live;
            // Below every lane's first limb there is nothing to add and no
            // carry can exist yet, so the position costs only this compare.
            if (active == 0) continue;

            const __mmask8 at_lo = _mm512_cmpeq_epi64_mask(off, pv) & live;
            const __mmask8 at_hi = _mm512_cmpeq_epi64_mask(off1, pv) & live;

            __m512i addend = _mm512_maskz_mov_epi64(at_lo, lo);
            addend = _mm512_mask_mov_epi64(addend, at_hi, hi);
            addend =
                _mm512_mask_xor_epi64(addend, negative & active, addend, kOnes);

            // A negative lane enters its first limb with the +1 of ~A + 1.
            carry = _mm512_mask_mov_epi64(carry, at_lo & negative, kOne);

            std::uint64_t* dst = limbs[p] + row;
            const __m512i x = _mm512_loadu_si512(dst);
            const __m512i sum = _mm512_add_epi64(x, addend);
            const __mmask8 c1 = _mm512_cmplt_epu64_mask(sum, x);
            const __m512i sum2 = _mm512_add_epi64(sum, carry);
            const __mmask8 c2 = _mm512_cmplt_epu64_mask(sum2, sum);
            _mm512_storeu_si512(dst, sum2);

            carry = _mm512_maskz_set1_epi64(c1 | c2, 1);

            // Exact per-block termination without a shuffle reduction: stop
            // once no lane has addend left above p and nothing is propagating.
            // The ~A + 1 form inverts the carry's meaning, so the test differs
            // by sign -- an adding lane is done when its carry is 0, a
            // subtracting lane when its carry is 1.
            const __mmask8 pending = _mm512_test_epi64_mask(carry, carry);
            const __mmask8 more = _mm512_cmpgt_epi64_mask(off1, pv) & live;
            if (more == 0 && ((pending ^ negative) & live) == 0) break;
        }
    }
}

}  // namespace kernels
}  // namespace cbfp
