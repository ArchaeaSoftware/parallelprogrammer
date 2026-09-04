// The flat kernel surface: free functions over raw pointers, one variant per
// instruction set, selected once per column rather than per element.
//
// Accumulation is two passes over a column's doubles. The first only learns
// the exponent range, because the rescale/widen decision has to be made before
// any value can be added. The second fuses decomposition into the add, so no
// decomposed form is ever written to memory.
//
// Every kernel takes a *contiguous* column. Strided input is the caller's
// problem: a strided vector gather costs more than the work it feeds, so the
// caller stages such a column once rather than letting the kernels gather.
#pragma once

#include <cstddef>
#include <cstdint>

#include "truesum/survey.hpp"

namespace truesum {
namespace kernels {

// The public Survey is what the kernels produce; nothing here needs its own.
using Survey = ::truesum::Survey;

// `floor_exponent` is the column's current exponent. The survey needs the exact
// minimum only when the incoming values could drop below it; otherwise a cheap
// lower bound settles that no rescale is due and the significand is never
// touched. Pass LLONG_MAX to force the exact value, which is what a column
// with no scale yet requires, since it adopts whatever the survey returns.
using SurveyFn = Survey (*)(const double *, std::size_t, long long);

// Fused decompose-and-add. `column_exponent` already absorbs any power-of-two
// scale, so the shift is simply value_exponent - column_exponent.
// `first_limb` is a starting hint from the survey; a row block that has not
// reached its own first limb skips the position outright.
// `flags` reports where the batch contradicted what the column was sized for:
// 1 a non-finite value, 2 an exponent below the column's, 4 an addend reaching
// the sign bit or past the width, 8 a sum carried past the sign bit. Bits are
// OR-ed in, never cleared. A column sized from its own survey cannot raise any
// of them, so checking it there costs nothing and asserts the survey agreed
// with the accumulate; a column sized from metadata supplied by a producer is
// where they earn their keep.
//
// The first three are found in the addend before it is applied. The fourth is
// signed overflow at the top limb: the operands agree in sign and the result
// does not. It is checked only at that limb, and the carry chain only reaches
// that limb when an addend lands in it or a carry climbs to it, so the check
// runs exactly where overflow is possible. Bit 4 guarantees that every addend
// applied at the top limb has its sign bit clear, which is what makes the
// sign rule exact there.
using AccumulateFn = void (*)(std::uint64_t *const *, std::size_t,
                              const double *, std::size_t, std::int32_t,
                              std::size_t, unsigned *);

inline constexpr unsigned kBadNonFinite = 1;
inline constexpr unsigned kBadExponent = 2;
inline constexpr unsigned kBadWidth = 4;
inline constexpr unsigned kBadOverflow = 8;

// Folds `count` contiguous columns, one per matrix, into a single pass over
// the accumulator. Same arithmetic as AccumulateFn applied `count` times; what
// differs is the traffic. One batch at a time reads and writes every limb it
// touches once per batch, so K batches pay 8 bytes of input and 16*nlimbs of
// accumulator per element. Folded, the limbs are read once, all K addends
// applied to them in registers, and written once: 8*K + 16*nlimbs for the same
// work. At two limbs and K=8 that is 96 bytes an element against 320.
//
// The saving is in accumulator traffic, so it only shows where that traffic is
// the bound -- past L3, where the CPU is waiting on DRAM rather than issue.
using AccumulateFoldFn = void (*)(std::uint64_t *const *, std::size_t,
                                  const double *const *, std::size_t,
                                  std::size_t, std::int32_t, unsigned *);

// A row's limbs are held in registers across the fold, so the width has to be
// bounded. Wider columns fall back to one batch at a time, which is no loss:
// the single-batch kernel stops as soon as a carry dies, while the fold has to
// write back every limb it loaded, so past this width folding would move more
// memory rather than less.
inline constexpr std::size_t kMaxFoldLimbs = 8;

Survey
survey_column_scalar(const double *values, std::size_t rows,
                     long long floor_exponent);

void
accumulate_scalar(std::uint64_t *const *limbs, std::size_t nlimbs,
                  const double *values, std::size_t rows,
                  std::int32_t column_exponent, std::size_t first_limb,
                  unsigned *flags);

void
accumulate_fold_scalar(std::uint64_t *const *limbs, std::size_t nlimbs,
                       const double *const *columns, std::size_t count,
                       std::size_t rows, std::int32_t column_exponent,
                       unsigned *flags);

#if defined(TRUESUM_HAVE_AVX512)
Survey
survey_column_avx512(const double *values, std::size_t rows,
                     long long floor_exponent);

void
accumulate_avx512(std::uint64_t *const *limbs, std::size_t nlimbs,
                  const double *values, std::size_t rows,
                  std::int32_t column_exponent, std::size_t first_limb,
                  unsigned *flags);

void
accumulate_fold_avx512(std::uint64_t *const *limbs, std::size_t nlimbs,
                       const double *const *columns, std::size_t count,
                       std::size_t rows, std::int32_t column_exponent,
                       unsigned *flags);
#endif

// Chosen once, on first use, from the running CPU's capabilities.
SurveyFn
survey();
AccumulateFn
accumulate();
AccumulateFoldFn
accumulate_fold();
const char *
accumulate_name();

// Adds one already-decomposed value to a single row. Not hot; used by the
// element-at-a-time API. Returns true if the sum overflowed its width; the
// caller is responsible for the addend itself fitting.
bool
accumulate_one(std::uint64_t *const *limbs, std::size_t nlimbs, std::size_t row,
               std::uint64_t mantissa, std::size_t shift, bool negative);

// dst = src << shift, sign-extended. dst and src are distinct allocations.
void
shift_left(std::uint64_t *const *dst, std::size_t ndst,
           std::uint64_t *const *src, std::size_t nsrc, std::size_t rows,
           unsigned shift);

// dst[i] = sign of top[i], i.e. all ones when negative and zero otherwise.
void
sign_fill(std::uint64_t *dst, const std::uint64_t *top, std::size_t rows);

}  // namespace kernels
}  // namespace truesum
