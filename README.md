# column-block-fp

Exact accumulation of double-valued matrices into an arbitrary-precision
fixed-point matrix. No intermediate precision is ever lost: summing a stream of
`double` matrices produces the mathematically exact result, which is rounded
only once, when you read a cell back out.

Header-only interface over a small static library. C++17, no dependencies.

## The idea

A finite `double` is exactly `m · 2^e` for an integer `m` of at most 53 bits.
So a sum of doubles is exactly representable in fixed point, given enough bits.
The catch is *how many* bits: covering the whole `double` range at once takes
~2100 bits per entry, which is wasteful when a column's values actually live in
a narrow band of magnitudes.

This structure gives every **column** its own scale. Entry `(i, j)` is stored as
a two's complement integer `V` of `column_limbs(j)` 64-bit limbs, denoting

```
A(i, j) = V · 2^column_exponent(j)
```

A column pays only for the dynamic range its own data spans. In the bundled
demo, a 4×4 accumulator holding values from 10⁻²³ to 10²⁶ occupies 584 bytes
total, with per-column widths of 128–192 bits rather than a flat 2100.

Both column parameters adapt automatically, and only ever in the safe direction:

- **The exponent only decreases.** A value whose ulp falls below the column's
  current LSB triggers a rescale: every entry in the column is shifted left and
  the exponent drops to match. Existing bits move, none are discarded.
- **The width only increases.** Rather than detect overflow, the required width
  is *derived*: no entry can exceed `count · 2^max_addend_bits`, so that bound
  plus a sign bit says how many limbs are needed before a batch is applied. The
  kernels can then skip overflow checks entirely, which is what makes them
  vectorizable.

Incoming mantissas are normalized to odd before use, so the column exponent
tracks each value's *true* ulp rather than its `frexp` exponent. This matters a
lot in practice: a column of small integers stays at exponent 2⁰ and one limb
wide, instead of immediately sinking to 2⁻⁵².

## Usage

```cpp
#include "cbfp/column_accumulator.hpp"

cbfp::ColumnBlockMatrix acc(rows, cols);

acc.reserve_for(first_batch.data(), n_batches);   // optional: avoids rescaling
for (const auto& batch : batches) {
    acc.add_matrix_col_major(batch.data());       // column-major: preferred
    // acc.add_matrix(batch.data());              // row-major also works
}

double x = acc.to_double(i, j);                   // correctly rounded, once
std::string exact = acc.to_exact_decimal(i, j);   // every digit, never rounds
```

Also available: `add` / `sub` per element, `add_matrix_scaled_pow2` (exact
scaling by any power of two), `to_matrix` for bulk readback,
`is_exactly_representable`, and `describe()` for a per-column exponent/width
report.

`add_column(j, v)` accumulates a single column from `rows()` contiguous
doubles, so a caller can stream columns instead of holding a whole input
matrix. The accumulator only ever touches one column at a time, and this is
what lets it be far larger than any input that could be resident alongside it —
on a 46 GB machine, a square accumulator goes from 42,000 x 42,000 to
60,000 x 60,000 (verified) simply by feeding it a column at a time.

`reserve_for` is purely an optimization — it pre-sizes each column from a
representative matrix and a batch count so accumulation never rescales. Results
are bit-identical with or without it; the test suite asserts this.

### Readback semantics

`to_double` is correctly rounded, round-to-nearest ties-to-even, over the full
range: it overflows to `±inf` and passes through the subnormal range without
double rounding (the subnormal path rounds the exact value directly to a
multiple of 2⁻¹⁰⁷⁴ rather than rounding twice via 53 bits).

`to_exact_decimal` is exact and finite — a binary fixed-point value always has a
terminating decimal expansion, since `V · 2⁻ᵏ = (V · 5ᵏ) / 10ᵏ`.

### Errors and edge cases

- `inf` / `NaN` throw `std::domain_error`; out-of-range indices throw
  `std::out_of_range`.
- `-0.0` accumulates as an exact zero. The sign of zero is not tracked, so a
  cell that sums to zero reads back as `+0.0`.
- `set_zero()` clears the values but keeps each column's learned scale, so a
  reused accumulator does not re-pay for rescaling.

### Storage layout

Storage is **limb-major**: one allocation per limb position, holding that limb
of every row (`limbs[k][i]` is entry `i`'s `k`-th limb). Rows are the vector
axis and limb positions the sequential one, so a carry chain runs inside a lane
and never crosses lanes. Widening is then a pure append — existing limb arrays
are untouched — and each allocation is 64-byte aligned and skewed off cache-set
congruence. See [docs/simd-design.md](docs/simd-design.md).

Kernels are flat free functions selected once per column from the running CPU's
capabilities; `active_kernel()` reports which. Set `CBFP_KERNEL=scalar` to force
the portable path, which is how the two are cross-checked in testing.

Accumulating a column is two passes over its values. The first learns only the
exponent range, because the rescale and widen decisions have to be made before
anything can be added; the second fuses decomposition into the add, so no
decomposed form is ever written to memory.

Both passes want the column contiguous. `add_matrix_col_major` gives them that
directly, so it is the faster entry point — **1.6x** over row-major on a
4096x64 accumulation, with bit-identical results. A row-major input has its
columns strided by `cols` doubles, and a strided vector gather costs more on
Zen 4 than the decomposition it feeds, so such a column is staged into a
contiguous buffer once rather than gathered.

Note which layout is which: the *accumulator's* storage always keeps a column
contiguous, since `limbs[k]` holds limb `k` of every row. Only the caller's
matrix can be strided.

### Threading

Columns are independent in storage, but `add_matrix` stages each column through
a shared buffer, so a single `ColumnBlockMatrix` is **not** reentrant. Use one
accumulator per thread and merge, or serialize calls.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/test_cbfp     # or: ctest --test-dir build
./build/demo
```

## Style

[.clang-format](.clang-format) captures the house style: four-space indent, and
the opening brace on its own line for function definitions (classes,
namespaces, and control statements keep theirs on the same line). Short
accessors defined inside a class may stay on one line; the few that wrap do so
only because they exceed the 80-column limit.

This machine has clang-format at `/usr/lib/llvm-18/bin/clang-format` (it is not
on `PATH`):

```sh
/usr/lib/llvm-18/bin/clang-format -i include/cbfp/*.hpp src/*.cpp \
    tests/*.cpp examples/*.cpp
```

Loops are written in `init/check/step` form with a bound that is visibly
finite; `while (1)` / `for (;;)` are avoided entirely, and conditional `break`
or `return` inside a bounded loop is the preferred exit. The two grid-shaped
matrix literals in the tests are fenced with `// clang-format off` so their
row/column layout survives formatting.

## Tests

`tests/test_cbfp.cpp` is a dependency-free assertion runner covering:

- exact decomposition and round-trip of random bit patterns, subnormals,
  `DBL_MAX`, `DBL_MIN`, and both zeros;
- catastrophic cancellation (`1e300 + 1 - 1e300` → exactly `1`);
- exact decimal expansions checked digit-for-digit against known values;
- correct rounding at ties, including ties-to-even at `1 + 2⁻⁵³`, the
  subnormal boundary, and the `2⁻¹⁰⁷⁵` round-to-zero tie;
- exponent-lowering and width-growth across the full 2⁻¹⁰⁷⁴…2¹⁰²³ span;
- add-then-subtract-in-shuffled-order returning to exact zero;
- absorption of a million values whose ulp is far below the running total,
  where naive summation stalls completely.

At size, where indexing and scale bookkeeping are what can break:

- **257 x 9 cross-checked against scalars.** A 257-row column shares one
  exponent across every row, while a 1x1 accumulator picks the exponent that
  suits its single cell, so the two hold different integers. Every one of the
  2313 cells must still agree digit-for-digit on its exact decimal — which pins
  down both the row-stride arithmetic and the claim that a column's shared
  scale never changes what it holds.
- **257 x 9 cancelled back to zero** with the batches and the cells within them
  walked in reverse, so nothing is undone in the order it was applied.
- **512 x 64 with a distinct power-of-two scale per column.** Each cell sums
  integers below 2^53 at a fixed scale, so plain double addition is itself
  exact and serves as an independent reference for all 32768 cells.
- **Padded row strides** on `add_matrix`, `add_matrix_scaled_pow2`,
  `reserve_for`, and `to_matrix`. The padding is NaN, which `add` rejects, so
  reading it throws rather than quietly passing; `to_matrix` must leave the
  padding it writes around untouched.

The suite runs in 25ms (140ms under sanitizers), and everything passes under
ASan and UBSan.

The large tests were checked by mutation: forcing `row_stride` to be ignored is
caught only by the stride test — the rest of the suite passes it clean — and
corrupting the per-row offset in `rescale` fails 3370 checks across the three
large tests.

## Layout

| path | contents |
| --- | --- |
| [include/cbfp/limbs.hpp](include/cbfp/limbs.hpp) | low-level two's complement limb arithmetic |
| [include/cbfp/column_accumulator.hpp](include/cbfp/column_accumulator.hpp) | the `ColumnBlockMatrix` interface |
| [src/limbs.cpp](src/limbs.cpp) | shifts, widening, carry/borrow chains, decimal conversion |
| [src/column_accumulator.cpp](src/column_accumulator.cpp) | scale tracking, rescaling, rounding |
| [tests/test_cbfp.cpp](tests/test_cbfp.cpp) | test suite |
| [examples/demo.cpp](examples/demo.cpp) | exact vs. naive summation comparison |
