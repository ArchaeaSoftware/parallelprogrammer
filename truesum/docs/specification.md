# truesum — exact accumulation of `double` matrices

**Specification.** For the reasoning behind these choices, and the
measurements that decided them, see `simd-design.md`. This document states what
the library guarantees and what it requires; that one states why.

---

## 1. What it guarantees

Given a sequence of `double` values, the library computes their sum **exactly**
— no intermediate rounding, at any point, for any input sequence — and rounds
once, on readback.

Three properties follow, and they are the reason to use it:

- **Exactness.** `to_exact_decimal` returns the true sum as a decimal string.
  `to_double` returns the correctly rounded (round-to-nearest, ties-to-even)
  `double` nearest that sum.
- **Order independence.** Any permutation or regrouping of the same addends
  gives bit-identical results. Threading, batching and folding cannot change
  the answer.
- **Reproducibility.** The CPU and CUDA implementations produce identical
  stored limbs for identical inputs. This is asserted by the test suite, which
  compares limbs rather than rounded values.

What it does not do is make inaccurate inputs accurate. If the summands are
already rounded — `float` products, TF32 matmul results — exact accumulation
preserves that error exactly. It is a determinism and exactness tool, not a
precision recovery tool.

---

## 2. Representation

A matrix is stored as **column blocks**. Every column has its own scale.
Entry `(i, j)` is a two's-complement integer `V` of `column_limbs(j)` 64-bit
limbs, denoting

```
A(i, j) = V × 2^column_exponent(j)
```

A `double` is exactly `m × 2^e` for integer `m`, so every accumulation is exact
provided the column exponent is low enough and the block wide enough. Both
adapt, and only in the safe direction: **the exponent never rises and the width
never shrinks.** No stored bit is ever discarded.

Limbs are stored **limb-major**: one allocation per limb position, each holding
that limb for every row of the column. The carry chain therefore runs along
limb positions inside one lane, while rows are the vector axis. Widening a
column appends an allocation and sign-fills it; nothing already allocated
moves.

---

## 3. Sizing

A column must be wide enough that no entry can overflow. The bound is
**derived, not detected**:

```
bits = max_addend_bits + ceil_log2(count + 1) + 1
```

where `max_addend_bits` spans the lowest bit any addend reaches to the highest
bit the largest one occupies, `count` is the number of addends, and the final
`+1` is the sign. This is what lets the carry chain run without an overflow
test at every limb, which is what lets it vectorize. The top limb alone
includes a sign test (§7), so a width the bound did not in fact cover is
reported rather than wrapped.

For a caller sizing an accumulation matrix in advance, the same rule reads:

```
bits = significand + exponent_spread + ceil_log2(count + 1) + 1
```

**Exponent spread and count are interchangeable, bit for bit.** One binade of
spread costs exactly what a doubling of the count costs. Worked examples:

| addends | significand | spread | bits at count = 2^14 |
| --- | --- | --- | --- |
| e4m3 products | 8 | 36 (full range) | 51 — one limb |
| TF32 products | 22 | 0 (well scaled) | 37 — one limb |
| TF32 products | 22 | 508 (full range) | 545 — nine limbs |
| bf16 products | 16 | ~394 | ~425 — seven limbs |
| `double` | 53 | up to 2098 | up to ~2100 |

The dominant term is the exponent range of the inputs, not their mantissa
width. Scaling inputs into a narrow band is therefore a *sizing* decision as
much as a numerical one.

---

## 4. Vocabulary types

### `Survey` — 12 bytes

```cpp
struct Survey {
    int min_exponent;   // lowest true-ulp exponent among nonzero values
    int max_top;        // one past the highest bit position occupied
    bool any;           // false if every value was zero
    bool nonfinite;     // true if any value was inf or NaN
};
```

`min_exponent` is the exponent after normalizing each mantissa to odd — the
value's true unit in the last place, not what `frexp` reports. A column of
integers surveys as exponent 0, not −52.

Both extents are `int`: a `double`'s true-ulp exponent lies in [−1074, 1023]
and its top in [−1073, 1024], so 32 bits covers six orders of magnitude more
range than the format can produce.

A survey is 12 bytes per column against 8 bytes per matrix element — **one part
in 43690** at 65536×64. It is small enough to travel ahead of the data it
describes.

**Producers:**

```cpp
Survey survey_column(const double *values, std::size_t rows);
void survey_matrix_col_major(Survey *out, const double *b, std::size_t rows,
                             std::size_t cols, std::size_t col_stride = 0);
void survey_matrix(Survey *out, const double *b, std::size_t rows,
                   std::size_t cols, std::size_t row_stride = 0);
// both pointers must be device or mapped pinned memory; asynchronous
void survey_matrix_col_major_device(Survey *out, const double *b,
                                    std::size_t rows, std::size_t cols,
                                    std::size_t col_stride = 0,
                                    CUstream_st *stream = nullptr);
```

The host and device forms agree field for field. Which side computes a survey
cannot change how an accumulation matrix is sized from it.

### `Uplo`

```cpp
enum class Uplo { Lower, Upper };
```

Following LAPACK. `Lower` stores entries with `i >= j`, so column `j` holds
rows `j..n-1`; `Upper` stores `i <= j`, so column `j` holds rows `0..j`.

### `InputRead` — device only

A move-only handle owning a CUDA event, returned by every device submission.
It reports when the device has finished reading the submitted buffer. See §7.

---

## 5. `AccumulationMatrix` — host

### Construction

```cpp
AccumulationMatrix(rows, cols);                    // adaptive
AccumulationMatrix(n, uplo);                       // symmetric, adaptive
AccumulationMatrix(rows, cols, surveys);           // pre-sized
AccumulationMatrix(n, uplo, surveys);              // symmetric, pre-sized
```

Move-only. `surveys` is indexed by matrix, then by column; its outer size
declares how many matrices will arrive and supplies the count term.

Pre-sizing reduces the surveys to **one reservation per column** — minimum of
the minima, maximum of the maxima. Because that is an aggregate, the
accumulation matrix never learns which matrix it is being handed, and **submission
order is irrelevant**. Sized this way it never rescales, never widens, and
never surveys an incoming batch.

### Accumulation

```cpp
void add_matrix(b, row_stride = 0);         // row-major
void add_matrix_col_major(b, col_stride = 0);
void add_matrix_scaled_pow2(b, log2_scale, row_stride = 0);
void add_matrix_col_major_scaled_pow2(b, log2_scale, col_stride = 0);
void add_matrices_col_major(b, count, col_stride = 0);   // fold
void add_column(j, v);
void add_column_scaled_pow2(j, v, log2_scale);
void set_zero();
```

`*_scaled_pow2` multiplies by a power of two. This is the only scaling that
stays exact without widening the mantissa; it is folded into the column
exponent rather than applied to the values.

`add_matrices_col_major` applies `count` matrices in a single pass over the
accumulation matrix. Results are identical to a loop; only the traffic differs. It is
**column-major only** — folding row-major input would require staging `count`
columns per worker, which is the traffic it exists to avoid. Columns wider than
8 limbs fall back to one matrix at a time.

### Readback

```cpp
double to_double(i, j) const;
double to_double(double *residual, i, j) const;
void to_matrix(double *out, row_stride = 0) const;
void to_matrix_with_residual(double *out, double *residual,
                             row_stride = 0) const;
double to_double_mean(i, j, std::uint64_t count) const;
double to_double_mean(double *residual, i, j, std::uint64_t count) const;
void to_matrix_mean(double *out, std::uint64_t count, row_stride = 0) const;
void to_matrix_mean_with_residual(double *out, double *residual,
                                  std::uint64_t count, row_stride = 0) const;
std::string to_exact_decimal(i, j) const;
bool is_exactly_representable(i, j) const;
bool is_zero(i, j) const;
std::vector<limb_t> entry_limbs(i, j) const;
```

`to_double` overflows to ±inf and underflows through the denormal range
without double rounding.

The **residual** is exactly `(stored value − returned double)`, itself
correctly rounded. It is computed by limb subtraction in the accumulation matrix's own
fixed point, never in floating point. Guaranteed properties:

- the returned `double` is unaffected by asking for the residual;
- `|residual| <= ulp(result) / 2`, negative exactly when the rounding went up,
  so the pair is a non-overlapping two-term expansion of about 106 bits;
- the residual is `+0.0` exactly when `is_exactly_representable` is true;
- a value outside `double`'s range returns ±inf with a residual of `0.0`.

The pair must be kept **unevaluated**. `hi + lo` in `double` arithmetic returns
`hi` unchanged in almost every case, because `|lo| <= ulp(hi)/2` by
construction.

`to_double_mean` is the correctly rounded value of `(stored value / count)`,
the exact rational, with the same range behaviour as `to_double`. The count is
supplied, never inferred: the accumulation matrix does not know how many values a sum
stands for. `count == 0` raises `std::domain_error`. A power-of-two count is
bit-identical, residual included, to having accumulated the input scaled by
that power. The mean's residual is exactly `(stored value / count − returned
double)`, correctly rounded, formed as an integer numerator over the same
count and divided once; every property listed above holds for it, with "the
residual is `+0.0` exactly when the quotient is a double" in place of the
`is_exactly_representable` clause. One boundary is worth knowing: a quotient
that falls exactly halfway between two denormals leaves a residual of half
the smallest denormal, which no double can represent, and that residual reads as
`+0.0`. There is no exact decimal of a mean, because there is none in general.

### Shape, sizing and tuning

```cpp
std::size_t rows() const, cols() const;
bool symmetric() const;  Uplo uplo() const;
std::size_t column_rows(j) const;        // rows() unless symmetric
std::size_t column_first_row(j) const;   // 0 unless symmetric and Lower
int column_exponent(j) const;
std::size_t column_limbs(j) const, column_bit_width(j) const;
void reserve_column(j, exponent, bits);
void reserve_for(b, count = 1, row_stride = 0);
std::size_t memory_bytes() const;
std::string describe() const;
void set_threads(unsigned n);   unsigned threads() const;
```

`reserve_*` is an optimization only; results are identical without it.

### Symmetric storage

`A(i, j)` and `A(j, i)` are **the same stored entry**, not two that agree.
This matters beyond halving the memory: every column has its own exponent and
width, so under full storage the two halves would hold equal values in
different limbs, and symmetry would be a property the data had to keep earning.
Stored once, it cannot drift — `entry_limbs(i, j) == entry_limbs(j, i)`
identically.

Accumulation reads only the stored triangle of the input, so input traffic
halves with the storage. **The input is assumed symmetric and this is not
checked**, because checking it means reading the half the design exists to
avoid reading.

`add_column(j, v)` takes `column_rows(j)` values beginning at logical row
`column_first_row(j)`.

### Threading

`set_threads(n)` spreads accumulation across `n` workers, partitioned by
column. Columns are independent in storage, so no locking is required and
results are **bit-identical to single-threaded**. For symmetric matrices the
partition is weighted by stored entries, since a triangular column's work runs
from `n` down to 1.

A `AccumulationMatrix` is not internally synchronized: one writer at a time.

---

## 6. `CudaAccumulationMatrix` — device

Not a third kernel variant behind the host dispatch. What the two
implementations share is the *algorithm* — when to rescale, when to widen, how
the exponent moves — not the memory. `cuda_available()` reports whether a
usable device exists; everything else throws if not.

The header forward-declares CUDA's opaque handle types rather than including
`cuda_runtime.h`, so it is usable from a translation unit with no CUDA toolkit
on its include path. Only the link needs `cudart`.

### Construction and reservation

```cpp
CudaAccumulationMatrix(rows, cols);
CudaAccumulationMatrix(n, uplo);
CudaAccumulationMatrix(rows, cols, surveys);
CudaAccumulationMatrix(n, uplo, surveys);

void reserve_column(j, exponent, bits);
void reserve_like(const AccumulationMatrix &cpu);
void reserve_for(b, count = 1, col_stride = 0);          // host input
void reserve_for_device(b, count = 1, col_stride = 0);   // device input
```

Non-copyable and non-movable. A device column may be reserved **once**;
`reserve_like` requires matching dimensions *and* matching symmetry and uplo,
since a triangular column is a different length.

Every column must be reserved before any accumulation.

### Accumulation

```cpp
double *acquire_input();
InputRead add_matrix_col_major(b, col_stride = 0);              // host input
InputRead add_matrix_col_major_device(b, col_stride = 0);       // device input
InputRead add_matrices_col_major_device(b, count, col_stride = 0);
void synchronize() const;
```

**Host input must be page-locked and device-mapped** — `cudaHostAlloc` with
`cudaHostAllocMapped`, `cudaHostRegister` with `cudaHostRegisterMapped`, or a
buffer from `acquire_input()`. Pageable memory is rejected, not copied. The
kernel reads the buffer in place over PCIe; there is no staging copy and no
device-side copy of the input.

`add_matrices_col_major_device` folds up to **16** device-resident matrices
into one pass over the accumulation matrix.

Accumulation is asynchronous. `synchronize()` blocks until every submission has
completed, and throws if any column reported a contradiction (§7).

### Readback

```cpp
std::vector<limb_t> entry_limbs(i, j) const;
std::vector<limb_t> download_column(j) const;   // out[k*rows + i]
int column_exponent(j) const;  std::size_t column_limbs(j) const;
int column_occupancy(j) const;
unsigned column_contradictions(j) const;
std::size_t memory_bytes() const;
```

Readback synchronizes. `column_occupancy` is the highest limb position the
accumulate has ever disturbed, measured rather than derived, monotonic across
batches, `-1` before anything is added. It has no consumer inside the library
and exists as a diagnostic: a caller can use it to see whether the surveys they
supplied were over-conservative.

There is no `to_double` on the device container. Convert by reading limbs, or
accumulate a comparison on the host.

---

## 7. Contracts

### Lifetime of device input

The kernel streams the submitted buffer over PCIe and **is still reading it
when the call returns**. Do not write to that buffer until the returned
`InputRead` reports the read has finished:

```cpp
InputRead h = acc.add_matrix_col_major_device(buf);
h.wait();           // now buf may be overwritten
```

Dropping the handle without waiting is a statement that the buffer will not be
written again. Ignoring this corrupts results — measured at about 63% of
entries.

`acquire_input()` returns a buffer the accumulation matrix owns and blocks only if the
device is still reading it, so two in rotation let the host fill one while the
device streams the other. That buffer stays valid until the next
`acquire_input()`.

**Fill on your own CUDA stream, not the default one.** The accumulation matrix's stream
is a blocking stream and therefore synchronizes with the legacy null stream: a
producer using plain `cudaMemcpy` serializes against the accumulate and gets no
overlap at all.

### Trust in supplied surveys

Pre-sizing takes the metadata on trust. A matrix reaching outside what was
declared, or more matrices than were described, voids the sizing.

- **The count is checked**, because that costs nothing.
- **The extents are not prevented**, because preventing them means surveying,
  which is what pre-sizing removed.

Instead the accumulate kernels **detect** contradictions using values already
in registers, at 0.8–4.8% of the kernel:

| bit | meaning |
| --- | --- |
| 1 | a non-finite value |
| 2 | an exponent below the column's |
| 4 | an addend whose top reaches the column's sign bit, or lies past its width |
| 8 | a sum carried past the sign bit |

On the host, any of these conditions cause the `std::runtime_error` exception to be thrown. On the
device, they are reported at the next `synchronize()`, or read without throwing via `column_contradictions(j)`.

**A detected contradiction means the sums are incorrect.** The accumulation matrix
is left inconsistent by design; the library reports rather than recovers. The
failure being guarded against is silent: an exponent below the column's makes
the shift negative, the conversion to unsigned puts the limb offset past the
width, and the value is dropped.

Bits 1, 2 and 4 are found in the addend before it is applied. Bit 8 is signed
overflow at the top limb — the operands agree in sign and the result does not
— found as it happens. It is tested only at that limb, and the carry chain
reaches that limb only when an addend lands in it or a carry climbs to it, so
the test runs exactly where overflow is possible. Bit 4 is what makes the test
exact: an addend whose top reaches the sign bit is rejected before it is
applied, so every addend that does reach the top limb has its sign bit clear.
The sum it catches is the one a survey understated by a few bits produces —
every addend fits, and the total does not — which was the one contradiction
the other three could not see.

Cost, on the CPU in §9, 4096×64, one thread: not measurable on AVX-512 (the
check is one ternary-logic op on the top limb and the run-to-run spread
swallows it). The scalar fallback pays 3% at two limbs and 27% at three under
GCC 9, whose carry loop is at the edge of what it keeps in registers; the
one-limb difference is code alignment, not the check. The device path is
bandwidth-bound and was not measured.

### Undefined and unchecked

- Passing a non-symmetric matrix to a symmetric accumulation matrix.
- Writing a device input buffer before its `InputRead` clears.
- Concurrent mutation of one accumulation matrix from several threads.

---

## 8. Errors

| exception | raised when |
| --- | --- |
| `std::invalid_argument` | malformed surveys; pageable memory where pinned is required; more than 16 folded matrices |
| `std::out_of_range` | row or column index out of range |
| `std::domain_error` | inf or NaN accumulated; `log2_scale` out of range; a mean over a count of zero |
| `std::runtime_error` | a detected contradiction; unreserved device column; more matrices than declared; any failed CUDA call |
| `std::bad_alloc` | allocation failure |

Every checked CUDA call reports the function and source location in its message.
Destructors do not throw.

---

## 9. Performance

Measured on a Ryzen 7 7700X (8 cores, 32 MiB L3, AVX-512 on a 256-bit
datapath) and an RTX 3060 (12 GB, 326.5 GB/s measured, 26.7 GB/s over PCIe).
Gelem/s is exact accumulations per second. **Indicative, not contractual.**

| | 4096×64 | 16384×64 | 65536×64 |
| --- | --- | --- | --- |
| CPU, 1 thread | 0.80 | 0.77 | 0.67 |
| CPU, 8 threads | 4.31 | 4.65 | 1.11 |
| CPU, 8 threads, folded K=8 | — | — | 3.01 |
| GPU, host input | 3.11 | 3.26 | 3.32 |
| GPU, input resident | 4.40 | 6.07 | 8.13 |
| GPU, resident, folded K=8 | — | — | 25.3 |

Which resource binds, and therefore which lever helps:

- The **CPU is issue-bound** while the working set fits in L3 — 2.32 IPC, 1%
  miss rate, scaling six-fold across eight cores — and **DRAM-bound** past it,
  where eight cores scale only 1.7×. Folding is the lever there.
- The **GPU is bandwidth-bound throughout**, at 97–99% of the card's streaming
  rate. Every gain comes from moving fewer bytes; none from faster arithmetic.
- **PCIe caps host input** at about 3.3 Gelem/s regardless of kernel quality.
  Device-resident input is 2.35–2.45× better, and that gap is the entire
  argument for delivering data over GPUDirect.

The crossover between the two targets is cache capacity, not compute.

Symmetric storage halves the memory exactly and is worth 2.0–2.2× on the CPU
(slightly over 2× where the triangle regains cache residency) and 1.4–1.5× on
the GPU host path.

---

## 10. Portability

- **C++17.** No exceptions to that; no C++20 features.
- **CPU kernels** are selected once, at first use, from the running CPU's
  capabilities: a scalar path everywhere, AVX-512 (F/BW/DQ/VL, VBMI2,
  VPOPCNTDQ) where present. Nothing compiled for AVX-512 runs before that
  check. `active_kernel()` reports which was chosen.
- **The public API has no templates** and exposes no SIMD or CUDA types.
- **64-bit limbs on every target.** 32-bit limbs measured slower on both.
- **Radix 64, canonical.** Carry-save at radix 52 was implemented and measured:
  1.33–1.55× at 2–4 limbs, but 0.81–0.96× at one limb and collapsing to 1.18×
  at eight once past L3. Not a global win, and the narrowest columns are the
  most common. The kernels remain templated on the radix so the hook survives.

