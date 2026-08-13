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

// Eight doubles from a contiguous column. There is deliberately no gather
// here: the caller stages a strided column instead, because gathering costs
// more than the decomposition it would feed.
inline __m512d load_column(const double* values, std::size_t row)
{
    return _mm512_loadu_pd(values + row);
}

// The final partial block, which is the only one that could read past the
// caller's array. Masked-off lanes read as zero, which is already inert
// everywhere downstream -- a zero has no mantissa, so it never goes live.
//
// Every loop below peels this block out rather than masking all of them: the
// mask is loop-invariant for all but the last block, and threading it through
// the body costs ~4% of the accumulate and ~15% of the survey, whose pass-1
// body is short enough for three extra instructions to matter.
inline __m512d load_column_tail(const double* values, std::size_t row,
                                __mmask8 m_k)
{
    return _mm512_maskz_loadu_pd(m_k, values + row);
}

struct Split8 {
    __m512i v_mantissa;  // odd
    __m512i v_exponent;  // true-ulp exponent
    __m512i v_top;       // one past the highest bit position reached
    __mmask8 m_negative;
    __mmask8 m_live;  // nonzero mantissa
    __mmask8 m_nonfinite;
};

// Vector form of the scalar `split`: the same IEEE-754 field extraction.
inline Split8 split8(__m512d v_val)
{
    const __m512i v_bits = _mm512_castpd_si512(v_val);
    const __m512i v_biased = _mm512_and_si512(_mm512_srli_epi64(v_bits, 52),
                                              _mm512_set1_epi64(0x7FF));
    const __m512i v_frac =
        _mm512_and_si512(v_bits, _mm512_set1_epi64(kFracMask));

    const __mmask8 m_is_normal = _mm512_test_epi64_mask(v_biased, v_biased);
    const __m512i v_m = _mm512_mask_or_epi64(
        v_frac, m_is_normal, v_frac, _mm512_set1_epi64(std::uint64_t{1} << 52));
    const __m512i v_e =
        _mm512_sub_epi64(_mm512_max_epu64(v_biased, _mm512_set1_epi64(1)),
                         _mm512_set1_epi64(1075));

    Split8 s;
    s.m_live = _mm512_test_epi64_mask(v_m, v_m);
    s.m_negative = _mm512_movepi64_mask(v_bits) & s.m_live;
    s.m_nonfinite = _mm512_cmpeq_epi64_mask(v_biased, _mm512_set1_epi64(0x7FF));
    // top = e + the significand's bit width, which odd-normalization leaves
    // unchanged: it raises e and lowers the width by the same amount. That is
    // why the survey never needs a trailing-zero count.
    s.v_top = _mm512_add_epi64(
        v_e, _mm512_sub_epi64(_mm512_set1_epi64(64), _mm512_lzcnt_epi64(v_m)));

    // v_tz = popcount((v_m & -v_m) - 1); srlv by 64 gives 0, so zero lanes
    // stay zero.
    const __m512i v_lowbit =
        _mm512_and_si512(v_m, _mm512_sub_epi64(_mm512_setzero_si512(), v_m));
    const __m512i v_tz =
        _mm512_popcnt_epi64(_mm512_sub_epi64(v_lowbit, _mm512_set1_epi64(1)));
    s.v_mantissa = _mm512_srlv_epi64(v_m, v_tz);
    s.v_exponent = _mm512_add_epi64(v_e, v_tz);
    return s;
}

inline __mmask8 tail_mask(std::size_t row, std::size_t rows)
{
    return row + 8 <= rows ? static_cast<__mmask8>(0xFF)
                           : static_cast<__mmask8>((1u << (rows - row)) - 1u);
}

}  // namespace

Survey survey_column_avx512(const double* values, std::size_t rows,
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
    __m512i v_min = _mm512_set1_epi64(std::numeric_limits<long long>::max());
    __m512i v_max_abs = _mm512_setzero_si512();
    __mmask8 m_any = 0;
    __mmask8 m_bad = 0;

    const __m512i kAbs = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFFll);
    const __m512i kExp = _mm512_set1_epi64(0x7FF);
    const __m512i kOne = _mm512_set1_epi64(1);

    // `m_k` is all-ones for every block but the last, so it is passed as a
    // constant here and the whole body folds down to the unmasked form.
    const auto exponent_step = [&](__m512i v_bits, __mmask8 m_k) {
        const __m512i v_abs_bits = _mm512_and_si512(v_bits, kAbs);
        const __m512i v_biased = _mm512_srli_epi64(v_abs_bits, 52);

        const __mmask8 m_live =
            _mm512_test_epi64_mask(v_abs_bits, v_abs_bits) & m_k;
        m_bad |= _mm512_cmpeq_epi64_mask(v_biased, kExp) & m_k;
        m_any |= m_live;

        v_min = _mm512_mask_min_epu64(v_min, m_live, v_min,
                                      _mm512_max_epu64(v_biased, kOne));
        v_max_abs =
            _mm512_mask_max_epu64(v_max_abs, m_live, v_max_abs, v_abs_bits);
    };

    std::size_t row = 0;
    for (; row + 8 <= rows; row += 8) {
        exponent_step(_mm512_castpd_si512(load_column(values, row)), 0xFF);
    }
    if (row < rows) {
        const __mmask8 m_k = tail_mask(row, rows);
        exponent_step(_mm512_castpd_si512(load_column_tail(values, row, m_k)),
                      m_k);
    }

    Survey out{0, 0, m_any != 0, m_bad != 0};
    if (!out.any || out.nonfinite) return out;

    const long long min_raw =
        static_cast<long long>(_mm512_reduce_min_epu64(v_min)) - 1075;
    const std::uint64_t max_abs =
        static_cast<std::uint64_t>(_mm512_reduce_max_epu64(v_max_abs));
    const std::uint64_t max_biased = max_abs >> 52;
    out.max_top = 0 != max_biased ? static_cast<long long>(max_biased) - 1022
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
    __m512i v_ulp = _mm512_set1_epi64(std::numeric_limits<long long>::max());

    const auto ulp_step = [&](__m512i v_bits, __mmask8 m_k) {
        const __m512i v_biased =
            _mm512_and_si512(_mm512_srli_epi64(v_bits, 52), kExp);
        const __m512i v_frac = _mm512_and_si512(v_bits, kFrac);
        const __mmask8 m_is_normal = _mm512_test_epi64_mask(v_biased, v_biased);
        const __m512i v_m =
            _mm512_mask_or_epi64(v_frac, m_is_normal, v_frac, kImplicit);
        const __mmask8 m_live = _mm512_test_epi64_mask(v_m, v_m) & m_k;

        // v_tz = popcount((v_m & -v_m) - 1)
        const __m512i v_e =
            _mm512_sub_epi64(_mm512_max_epu64(v_biased, kOne), k1075);
        const __m512i v_lowbit = _mm512_and_si512(
            v_m, _mm512_sub_epi64(_mm512_setzero_si512(), v_m));
        const __m512i v_tz =
            _mm512_popcnt_epi64(_mm512_sub_epi64(v_lowbit, kOne));
        v_ulp = _mm512_mask_min_epi64(v_ulp, m_live, v_ulp,
                                      _mm512_add_epi64(v_e, v_tz));
    };

    row = 0;
    for (; row + 8 <= rows; row += 8) {
        ulp_step(_mm512_castpd_si512(load_column(values, row)), 0xFF);
    }
    if (row < rows) {
        const __mmask8 m_k = tail_mask(row, rows);
        ulp_step(_mm512_castpd_si512(load_column_tail(values, row, m_k)), m_k);
    }
    out.min_exponent = _mm512_reduce_min_epi64(v_ulp);
    return out;
}

// Eight rows per pass. Rows are the vector axis and limb positions the
// sequential one, so the carry chain runs inside a lane and never crosses
// lanes -- no shuffles, no cross-lane propagation.
//
// Decomposition is fused in: mantissa and exponent go straight from the loaded
// doubles into the addend and are never written to memory.
//
// Two blocks are decomposed before either one is applied. Consecutive blocks
// are independent and the decomposition is a long dependency chain, so this
// lets the chains interleave. Measured in isolation it is worth nothing at one
// limb, ~7% at two and ~10% at four; end to end it is not distinguishable from
// one block at a time, because the accumulate is only part of the pipeline. It
// is kept for the wide-column case, where it also cuts run-to-run variance
// noticeably. Unrolling further does not pay: 4x matches 2x within noise and
// 8x regresses, as the staged addends start costing registers.
//
// Add and subtract lanes share one loop: subtracting A is adding ~A + 1, so a
// negative lane complements its addend and enters its first limb with carry 1.
// Limbs below a lane's offset stay untouched because the complement is applied
// only where the lane is active.
namespace {

// One block's worth of decomposed values, ready to be added.
struct Addend8 {
    __m512i v_off, v_off1, v_lo, v_hi;
    __mmask8 m_live, m_negative;
    std::size_t row;
};

}  // namespace

void accumulate_avx512(std::uint64_t* const* limbs, std::size_t nlimbs,
                       const double* values, std::size_t rows,
                       std::int32_t column_exponent, std::size_t first_limb)
{
    const __m512i kOne = _mm512_set1_epi64(1);
    const __m512i kOnes = _mm512_set1_epi64(-1);
    const __m512i k63 = _mm512_set1_epi64(63);
    const __m512i k64 = _mm512_set1_epi64(64);
    const __m512i kColExp = _mm512_set1_epi64(column_exponent);

    const auto prepare = [&](std::size_t row, __m512d v_val, __mmask8 m_k) {
        const Split8 s = split8(v_val);

        Addend8 q;
        q.row = row;
        q.m_live = s.m_live & m_k;
        q.m_negative = s.m_negative & q.m_live;
        // Dead lanes get shift 0 so their limb offset cannot go negative.
        const __m512i v_shift =
            _mm512_maskz_sub_epi64(q.m_live, s.v_exponent, kColExp);
        q.v_off = _mm512_srli_epi64(v_shift, 6);
        const __m512i v_bit = _mm512_and_si512(v_shift, k63);
        // The 53-bit mantissa lands in at most two limbs. srlv by 64 yields 0,
        // which is exactly what the v_bit == 0 case wants.
        q.v_lo = _mm512_maskz_sllv_epi64(q.m_live, s.v_mantissa, v_bit);
        q.v_hi = _mm512_maskz_srlv_epi64(q.m_live, s.v_mantissa,
                                         _mm512_sub_epi64(k64, v_bit));
        q.v_off1 = _mm512_add_epi64(q.v_off, kOne);
        return q;
    };

    const auto apply = [&](const Addend8& q) {
        if (0 == q.m_live) return;
        __m512i v_carry = _mm512_setzero_si512();

        for (std::size_t p = first_limb; p < nlimbs; ++p) {
            const __m512i v_pv = _mm512_set1_epi64(static_cast<long long>(p));
            const __mmask8 m_active =
                _mm512_cmple_epi64_mask(q.v_off, v_pv) & q.m_live;
            // Below every lane's first limb there is nothing to add and no
            // carry can exist yet, so the position costs only this compare.
            if (0 == m_active) continue;

            const __mmask8 m_at_lo =
                _mm512_cmpeq_epi64_mask(q.v_off, v_pv) & q.m_live;
            const __mmask8 m_at_hi =
                _mm512_cmpeq_epi64_mask(q.v_off1, v_pv) & q.m_live;

            __m512i v_addend = _mm512_maskz_mov_epi64(m_at_lo, q.v_lo);
            v_addend = _mm512_mask_mov_epi64(v_addend, m_at_hi, q.v_hi);
            v_addend = _mm512_mask_xor_epi64(v_addend, q.m_negative & m_active,
                                             v_addend, kOnes);

            // A negative lane enters its first limb with the +1 of ~A + 1.
            v_carry =
                _mm512_mask_mov_epi64(v_carry, m_at_lo & q.m_negative, kOne);

            std::uint64_t* dst = limbs[p] + q.row;
            const __m512i v_x = _mm512_loadu_si512(dst);
            const __m512i v_sum = _mm512_add_epi64(v_x, v_addend);
            const __mmask8 m_c1 = _mm512_cmplt_epu64_mask(v_sum, v_x);
            const __m512i v_sum2 = _mm512_add_epi64(v_sum, v_carry);
            const __mmask8 m_c2 = _mm512_cmplt_epu64_mask(v_sum2, v_sum);
            _mm512_storeu_si512(dst, v_sum2);

            v_carry = _mm512_maskz_set1_epi64(m_c1 | m_c2, 1);

            // Exact per-block termination without a shuffle reduction: stop
            // once no lane has addend left above p and nothing is propagating.
            // The ~A + 1 form inverts the carry's meaning, so the test differs
            // by sign -- an adding lane is done when its carry is 0, a
            // subtracting lane when its carry is 1.
            const __mmask8 m_pending = _mm512_test_epi64_mask(v_carry, v_carry);
            const __mmask8 m_more =
                _mm512_cmpgt_epi64_mask(q.v_off1, v_pv) & q.m_live;
            if (0 == m_more && 0 == ((m_pending ^ q.m_negative) & q.m_live))
                break;
        }
    };

    std::size_t row = 0;
    for (; row + 16 <= rows; row += 16) {
        const Addend8 a = prepare(row, load_column(values, row), 0xFF);
        const Addend8 b = prepare(row + 8, load_column(values, row + 8), 0xFF);
        apply(a);
        apply(b);
    }
    for (; row + 8 <= rows; row += 8) {
        apply(prepare(row, load_column(values, row), 0xFF));
    }
    if (row < rows) {
        const __mmask8 m_k = tail_mask(row, rows);
        apply(prepare(row, load_column_tail(values, row, m_k), m_k));
    }
}

}  // namespace kernels
}  // namespace cbfp
