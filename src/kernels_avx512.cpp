// Compiled with AVX-512 flags. Nothing here may run before the CPU check in
// kernels_dispatch.cpp -- no static initializers, no unconditional entry.
#include <immintrin.h>

#include "cbfp/limb_column.hpp"
#include "kernels.hpp"

namespace cbfp {
namespace kernels {

// Eight rows per pass. Rows are the vector axis and limb positions the
// sequential one, so the carry chain runs inside a lane and never crosses
// lanes -- no shuffles, no cross-lane propagation.
//
// Add and subtract lanes share one loop: subtracting A is adding ~A + 1, so a
// negative lane complements its addend and enters its first limb with carry 1.
// Limbs below a lane's offset stay untouched because the complement is applied
// only where the lane is active.
void accumulate_avx512(std::uint64_t* const* limbs, std::size_t nlimbs,
                       const Batch& b)
{
    const std::size_t n = padded_rows(b.rows);
    const __m512i kOne = _mm512_set1_epi64(1);
    const __m512i kOnes = _mm512_set1_epi64(-1);
    const __m512i k63 = _mm512_set1_epi64(63);
    const __m512i k64 = _mm512_set1_epi64(64);

    for (std::size_t row = 0; row < n; row += 8) {
        const __m512i mant = _mm512_loadu_si512(b.mantissa + row);
        if (_mm512_test_epi64_mask(mant, mant) == 0) continue;  // all zero

        const __mmask8 live = _mm512_test_epi64_mask(mant, mant);
        // Zero-mantissa lanes carry no meaningful exponent, so clamp their
        // shift to 0 rather than let it go negative and wrap the limb offset.
        const __m512i shift = _mm512_maskz_mov_epi64(
            live, _mm512_cvtepi32_epi64(_mm256_sub_epi32(
                      _mm256_loadu_si256(
                          reinterpret_cast<const __m256i*>(b.exponent + row)),
                      _mm256_set1_epi32(b.column_exponent))));
        const __m512i negv = _mm512_cvtepu8_epi64(_mm_loadl_epi64(
            reinterpret_cast<const __m128i*>(b.negative + row)));
        const __mmask8 negative = _mm512_test_epi64_mask(negv, negv);

        const __m512i off = _mm512_srli_epi64(shift, 6);
        const __m512i bit = _mm512_and_si512(shift, k63);

        // The 53-bit mantissa lands in at most two limbs. srlv by 64 yields 0,
        // which is exactly what the bit == 0 case wants.
        const __m512i lo = _mm512_sllv_epi64(mant, bit);
        const __m512i hi = _mm512_srlv_epi64(mant, _mm512_sub_epi64(k64, bit));
        const __m512i off1 = _mm512_add_epi64(off, kOne);

        __m512i carry = _mm512_setzero_si512();
        for (std::size_t p = b.first_limb; p < nlimbs; ++p) {
            const __m512i pv = _mm512_set1_epi64(static_cast<long long>(p));
            const __mmask8 active = _mm512_cmple_epi64_mask(off, pv);
            // Below every lane's first limb there is nothing to add and no
            // carry can exist yet, so the position costs only this compare.
            if (active == 0) continue;

            const __mmask8 at_lo = _mm512_cmpeq_epi64_mask(off, pv);
            const __mmask8 at_hi = _mm512_cmpeq_epi64_mask(off1, pv);

            __m512i addend = _mm512_maskz_mov_epi64(at_lo, lo);
            addend = _mm512_mask_mov_epi64(addend, at_hi, hi);
            addend =
                _mm512_mask_xor_epi64(addend, negative & active, addend, kOnes);

            // A negative lane enters its first limb with the +1 of ~A + 1.
            carry = _mm512_mask_mov_epi64(carry, at_lo & negative, kOne);

            std::uint64_t* dst = limbs[p] + row;
            const __m512i x = _mm512_loadu_si512(dst);
            const __m512i s = _mm512_add_epi64(x, addend);
            const __mmask8 c1 = _mm512_cmplt_epu64_mask(s, x);
            const __m512i s2 = _mm512_add_epi64(s, carry);
            const __mmask8 c2 = _mm512_cmplt_epu64_mask(s2, s);
            _mm512_storeu_si512(dst, s2);

            carry = _mm512_maskz_set1_epi64(c1 | c2, 1);

            // Past every lane's addend, a lane is finished once its limbs stop
            // changing. The ~A + 1 form inverts the carry's meaning, so the
            // test differs by sign: an adding lane is done when its carry is 0,
            // a subtracting lane when its carry is 1 (carry 0 there means a
            // borrow is still propagating through the sign extension).
            // Exact per-block termination without a shuffle reduction: stop
            // once no lane has addend left above p and nothing is propagating.
            // The ~A + 1 form inverts the carry's meaning, so the test differs
            // by sign -- an adding lane is done when its carry is 0, a
            // subtracting lane when its carry is 1.
            const __mmask8 pending = _mm512_test_epi64_mask(carry, carry);
            const __mmask8 more = _mm512_cmpgt_epi64_mask(off1, pv);
            if (more == 0 && (pending ^ negative) == 0) break;
        }
    }
}

}  // namespace kernels
}  // namespace cbfp
