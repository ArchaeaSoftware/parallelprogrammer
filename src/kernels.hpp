// The flat kernel surface: free functions over raw pointers, one variant per
// instruction set, selected once per column rather than per element.
//
// Every kernel assumes the caller has already made the column wide enough and
// its exponent low enough, so none of them allocates or reallocates.
#pragma once

#include <cstddef>
#include <cstdint>

namespace cbfp {
namespace kernels {

// One column's worth of incoming values, already decomposed. A zero mantissa
// means the row contributes nothing.
struct Batch {
    const std::uint64_t* mantissa;  // magnitude, at most 53 significant bits
    const std::int32_t*
        exponent;  // absolute; the kernel subtracts the column's
    const std::uint8_t* negative;  // 1 = subtract
    std::int32_t column_exponent;
    std::size_t rows;
};

using AccumulateFn = void (*)(std::uint64_t* const*, std::size_t, const Batch&);

void accumulate_scalar(std::uint64_t* const* limbs, std::size_t nlimbs,
                       const Batch& batch);

#if defined(CBFP_HAVE_AVX512)
void accumulate_avx512(std::uint64_t* const* limbs, std::size_t nlimbs,
                       const Batch& batch);
#endif

// Chosen once, on first use, from the running CPU's capabilities.
AccumulateFn accumulate();

// Name of the selected variant, for reporting.
const char* accumulate_name();

// dst = src << shift, sign-extended. dst and src are distinct allocations.
void shift_left(std::uint64_t* const* dst, std::size_t ndst,
                std::uint64_t* const* src, std::size_t nsrc, std::size_t rows,
                unsigned shift);

// dst[i] = sign of top[i], i.e. all ones when negative and zero otherwise.
void sign_fill(std::uint64_t* dst, const std::uint64_t* top, std::size_t rows);

}  // namespace kernels
}  // namespace cbfp
