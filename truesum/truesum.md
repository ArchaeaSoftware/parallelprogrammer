`truesum`: Fast, Arbitrary-Precision Floating Point Summation

Floating point operations notoriously are not associative. In particular, (a+b)+c ≠ a+(b+c), introducing a degree of uncertainty about accuracy of workloads that rely on floating point summation.

To eliminate that uncertainty, our library `truesum` provides an *accumulation matrix* that can perform fast, elementwise exact summation of floating point matrices. Each column is represented as block floating point, with a shared exponent and the mantissas held in an SOA layout of 64-bit *limbs*[^limb], the segments of a multi-precision number that each fit in a single machine word.

A *survey* is a compact (12B per column) characterization of an input matrix that may be used to pre-size an accumulation matrix to receive its contents without overflow, rescaling, or the allocation of new column-limbs during the accumulation. Extending this idea, an accumulation matrix can be pre-sized for multiple input matrices by supplying one survey for each prospective input matrix.

A final feature of the accumulation matrix is that it can return the mean, exactly rounded, instead of requiring the client to divide a readback by the number of accumulations.

# Structure Of Arrays

The SOA layout enables both AVX-512 and CUDA implementations to process summations at the speed of memory bandwidth.

Within the matrix, a single element is stored as a two's complement integer spanning several 64-bit limbs. Written the usual way (most significant first), it reads `[ limb 2 | limb 1 | limb 0 ]`. The SOA layout provides for a separate, contiguous allocation per limb position, so limbs 0, 1 and 2 for a given element are stored at the same offset into three different arrays. For want of a better term, we refer to these arrays as *limb-columns*.

Due to the dynamic range of `double`, with its 10-bit exponent, the worst-case space requirements for a single column of `double` is an eye-popping 33 limbs, or 264 bytes per entry. At risk of stating the obvious, the library was designed with an eye toward better-conditioned inputs that drive more modest requirements.

```
                  row 0    row 1    row 2    row 3
                ┌────────┬────────┬────────┬────────┬─────┐
     limbs[2]   │ limb 2 │ limb 2 │ limb 2 │ limb 2 │ ... │   2^(e+128), sign
                ├────────┼────────┼────────┼────────┼─────┤
     limbs[1]   │ limb 1 │ limb 1 │ limb 1 │ limb 1 │ ... │   2^(e+64)
                ├────────┼────────┼────────┼────────┼─────┤
     limbs[0]   │ limb 0 │ limb 0 │ limb 0 │ limb 0 │ ... │   2^e
                └────────┴────────┴────────┴────────┴─────┘
                     ^
                     └─ limbs of one entry

                each limbs[k] is a limb-column: one allocation, every row

                └─ one vector instruction spans this way ─┘
```

Propagating carries along the limbs of a single entry is inherently sequential: limb `k + 1` cannot be computed until after limb `k`. Rows, by contrast, are completely independent of one another, so the limb-major memory layout enables rows to be processed by a single vector instruction (8 64-bit elements at a time on AVX-512), while the sequential carry chain becomes the loop over limb positions, one vector instruction per position.

# Per-Column Block Floating Point (Shared Exponents)

On the assumption that columns contain similarly scaled values, the library implements a block floating point scheme that shares an exponent per column: not the exponent of its largest value, but the position of the *least significant* bit of any of its values. This arrangement stands in contrast to IEEE standard representations whose exponents identify the position of the *most* significant set bit. Recall that floating point values have more than one exact decomposition into an integer times a power of two: 1.0 is exactly 1 × 2⁰, but just as exactly 2⁵² × 2⁻⁵².[^three] IEEE-754 uses the latter, because the most significant set bit is redundant: with the exception of denormals, the format encodes a hidden bit to increase its effective precision.

The power of two specified by a column's exponent, then, is each value's true unit in the last place, or ulp. Taking the minimum over true ulps rather than over standard exponents can save up to 52 bits per column: a column of small integers stays at exponent 2⁰ and one limb wide, whereas a `frexp`-based scale would reduce it to 2⁻⁵², forcing the allocation of a second limb holding nothing but zeros.

```
       limb 2                  limb 1                  limb 0
  ┌───────────────────────┬───────────────────────┬───────────────────────┐
  │S                      │                       │                       │
  └───────────────────────┴───────────────────────┴───────────────────────┘
   191                 128 127                  64 63                    0
   2^(e+191)                                                           2^e
```

A column maintains four parameters that must be updated with a survey before applying an accumulation[^presized]: `exponent`, `max_addend_bits`, `add_count`, and the number of limbs. `exponent` can only decrease; the other three only increase.

`max_addend_bits` is the greatest number of bits any addend in the column has spanned, counting from the column's current bit 0 up to that addend's most significant bit. `max_addend_bits` and `add_count`, a running count of the addends, both are updated from a survey of the incoming values, taken before any of those values is added. A rescale adjusts `max_addend_bits` too, adding the shift it applied, since moving the column's exponent down by that many bits widens every addend's span from the new bit 0 by the same amount. Both exist to bound the width the accumulated sum can reach, so a column whose entries are all zero needs neither. `max_addend_bits` and `add_count` are cleared when the caller empties the accumulation matrix for reuse, and when a rescale finds every entry already zero, which is also the one case where the exponent can change without shifting anything.

`exponent` decreases whenever a value arrives whose least significant set bit sits below the column's current `exponent`. Every entry is then left-shifted to match: existing bits move, and none are discarded.

An increase in the limb count, and allocation of one or more new limb columns, may be prompted by any of three events:

- A value whose most significant bit sits above the high water mark denoted by `max_addend_bits`,
- A lower exponent, since that prompts a rescale of `shift` bits, and
- More addends, since a longer run can sum to a larger total.

The sum of these considerations determines the number of bits needed for a column. `max_addend_bits` accounts for the largest single addend; `ceil_log2(add_count + 1)` bits are reserved above `max_addend_bits`, gaining one more each time the count doubles; and a final bit is reserved for the sign. Rounding the total to the next multiple of 64 gives the number of limbs needed for the allocation.

A batch of input matrices is surveyed before any are accumulated: for any given column, a first pass over the incoming values determines whether the column must rescaled and widened before any values are added.

Reallocation and rescaling of a column may be prompted by either very large or very small values.

Consider a column with 1.0, which sets `exponent=0` and allocates a single limb. Adding a large value cannot change `exponent`, because a `double` of magnitude 2ᵏ has its lowest set bit no lower than 2ᵏ⁻⁵² and so never reaches below it, but it may increase the number of limbs (2⁵⁰⁰ widens the column to 8 limbs and 10³⁰⁰ to 16, with `exponent` unchanged at 0). Adding a small value decreases `exponent` and also increases the limb count (2⁻⁵⁰⁰ drops `exponent` to -500 and widens the column to 8 limbs; 10⁻³⁰⁰ decreases `exponent` -1049 and widens the column to 17 limbs).

```
  Widen: append a limb position. Nothing already stored moves.

      before            [ limb 1 | limb 0 ]
      after    [ limb 2 | limb 1 | limb 0 ]
                  new    └── unchanged ──┘
                         └─ not read, not written, not moved

  Rescale: the window slides down, and every entry shifts left to match.

      before   [ 1 0 1 1 ]                 bit 0 has significance 2^0
      after    [ 1 0 1 1 0 0 0 ]           bit 0 has significance 2^-3
                         └── three new low bits. The stored 1 0 1 1 did
                             not change; only the significance of bit 0 did,
                             so the value is the same and nothing was lost.
```

A `double` arrives as a 53-bit odd mantissa placed at `shift`, the distance from its own true ulp down to the column's. Since 53<64, a `double` may span at most two limbs, so only `off` and `off + 1` need to be updated, with limbs above `off + 1` receiving carry propagation until it is absorbed.

```
  ┌───────────────────────┬───────────────────────┬───────────────────────┐
  │                       │        ┌──────────────┼──────────────┐        │
  │                       │        │   53-bit odd mantissa       │        │
  │                       │        └──────────────┼──────────────┘        │
  └───────────────────────┴───────────────────────┴───────────────────────┘
                                   └── off + 1 ───┴───── off ────┘

  bits   =  max_addend_bits  +  ceil_log2(add_count + 1)  +  1
            └ what one addend   └ headroom, so that          └ the sign bit S,
              spans, lowest bit   `add_count` of them          which no sum may
              to highest          cannot carry past            ever reach

  limbs  =  ceil(bits / 64)   <- and this is what gets allocated
```

The carry chain runs fully vectorized, without an overflow test until it reaches the most-significant limb, where a sign comparison of the operands is performed against the result.

A column need only allocate as many column-limbs as needed for the dynamic range spanned by its own data. In the bundled demo, a 4×4 accumulation matrix holding values from 10⁻²³ to 10²⁶ occupies 584 bytes, with per-column allocations of two and three limbs rather than the thirty-three that covering the whole `double` range at once would require.

# Surveys

Surveys supply the per-column metadata needed to know how many column-limbs are needed, and whether the column must be rescaled before the survey's matrix is submitted for accumulation.

A survey is 12B per column and contains four fields:

- `min_exponent`, the lowest true-ulp exponent among the column's nonzero values, taken after normalizing each mantissa to odd rather than from `frexp`, so a column of small integers reports 0 and not -52.
- `max_top`, one past the highest bit position any value occupies.
- `any`, true if any value was nonzero.
- `nonfinite`, true if any value was INF/NaN.

Both exponents are `int` rather than `long long`, since 32-bit exponents can represent six orders of magnitude more range than `double` can produce.

Computing a survey requires at most two sweeps over the column. The first masks off the sign bit, early-outs if any value is +-0.0 or INF/NaN, and computes the maximum of the remaining values, which later may be used to compute `max_top`. Additionally, the minimum biased exponent is computed, clamped at 1 so every denormal reports the same one. That minimum is a lower bound on the true ulp: a biased exponent locates the lowest place a significand could occupy, while the true ulp sits above it by the significand's trailing zero count. Unless the column held nothing but zeros, or a non-finite value ended the first sweep, or that bound is larger than a cutoff supplied by the caller, a second sweep of the column is needed to split each value to its odd mantissa and take the minimum in order to compute `min_exponent`.

An accumulation matrix constructed without surveys computes them itself, one per incoming matrix, and in that case the second sweep is often unnecessary because the accumulation matrix already knows the exponent set for each column. What it does with `min_exponent` is a single comparison: if the value is less than the column's exponent, it rescales the column down to the new minimum. So it passes the column's exponent to the kernel as the cutoff. A first-sweep bound at or above that cutoff answers the comparison on its own: no true ulp can fall below the bound, so no rescale is due, and the kernel returns the bound without splitting a single mantissa. The width does not depend on it either: `max_addend_bits` is `max_top` less the column's exponent, and the first sweep computes `max_top` exactly. A producer computing a survey to send elsewhere has no column to compare against, so it supplies no cutoff at all and the second sweep always runs. Neither does the accumulation matrix when it surveys a column that has no exponent yet: that column will adopt whatever the survey reports, and since a column's exponent only ever moves down, a bound would leave its scale below what its values require, for good.

```
  Computing a survey of one column, in two sweeps over its values:

    sweep 1   a shift and a compare on each value's raw bits
              |
              +--  largest |bits|       ->  max_top
              +--  biased == 0x7FF      ->  nonfinite, and the survey stops
              +--  lowest raw exponent  ->  a lower bound on the true ulp

    sweep 2   split each value to its odd mantissa
              |
              +--  the lowest bit set   ->  min_exponent

  An accumulation matrix surveying its own input hands the kernel the
  exponent already set for its column. When sweep 1's lower bound clears
  that cutoff, no rescale can be due and sweep 2 never runs.
```

There is an entry point for a single column, for a column-major matrix, and for a row-major one, which stages each column into contiguous memory first because a strided gather costs more than the work it feeds. A fourth surveys a matrix already resident on the device, for a producer whose data never touches host memory at all.

Sizing a column from one survey takes its exponent from `min_exponent` and its span from `max_top` minus `min_exponent`, then adds the headroom for however many matrices are coming. Sizing from several surveys of the same column takes the minimum of the minima and the maximum of the maxima. Reducing them to a single aggregate is what makes submission order irrelevant: any matrix inside those extents fits the reservation, so the accumulation matrix never has to know which one it is being handed.

Two adjustments may follow from a survey. A survey whose `min_exponent` falls below the column's `exponent` calls for a rescale: fresh limb-columns are allocated, every stored value is shifted left by the difference, and `max_addend_bits` is increased by that same shift. The exponent specifies the significance of the lowest bit position the column can hold, so increasing it would right-shift every stored value, dropping bits already accumulated; it only ever decreases. An all-zero column is the exception, having no information to lose, so its scale moves freely and its bookkeeping starts over.

A survey whose `max_top` exceeds what the column's limb-columns can hold will cause more limb-columns to be appended and initialized with a sign extension. In contrast, a rescale allocates a fresh set of limb-columns and shifts every stored value into them, shifting zeroes into the least significant bits.

The separation of concerns between matrices and their surveys enables the work to be done by different constituencies. In a data center, the nodes that computed the matrices can also compute their surveys and transmit them ahead to the accumulation node; once the accumulation matrix there has been pre-sized to receive them all, the matrices themselves follow, one matrix or even one column at a time. The final output is bit-identical regardless of they order they arrive or are processed.

```
  Ahead of the data, one survey per column -- 12 bytes each, whatever the
  row count:

      node 1  ── survey ──┐
      node 2  ── survey ──┤
        ...               ├──>  accumulation matrix sized once:
      node N  ── survey ──┘        no rescale, no widen, no survey, ever

  Then the matrices themselves, 8 bytes an element, in any order at all:

      node k  ══ matrix ═══════════════════>  accumulate  ══>  the same
      node 1  ══ matrix ═══════════════════>  accumulate  ══>  bits either
      node 3  ══ matrix ═══════════════════>  accumulate  ══>  way
```

# Accumulation

Since a survey has been performed and applied to the accumulation matrix ahead of time, the accumulation kernel need not concern itself with rescaling or allocation of column-limbs. Every column already has the correct exponent and limb count when the kernel starts, so can be updated in a single pass over the input. A value's lowest set bit sits some number of positions above the lowest position the column can hold. That distance, the value's true-ulp exponent minus the column's exponent, is the shift: divided by 64, it gives the limb where the addend starts, and the remainder gives its bit offset within that limb. A 53-bit significand occupies at most two limbs at any bit offset, and the carry propagates upward from the lower one.

A column allocated wide for its dynamic range therefore costs no more per addend than a narrow one, except as required by the input values.

On a pre-sized accumulation matrix, the exponents and limb counts were set at construction from the surveys, which are verified as the accumulation kernel runs: three tests on an addend before it is applied, and a fourth on the running total as it crosses into the top limb.

Between them, those tests detect any deficiencies in the surveys. A survey can understate either of a column's extents, or fail to report a non-finite value; separately, a caller can supply fewer surveys than the matrices it goes on to submit. Each of those throws `std::runtime_error` rather than quietly returning an incorrect sum. They are designed to err on the side of caution: a total that wraps past the sign bit and later wraps back is correct by modular arithmetic, but is reported anyway.

An addend below the column's exponent or too wide for it is skipped; a non-finite one is recorded and applied anyway; and an overflowed sum is left where it wrapped. The rest of the batch goes in either way, and the error is raised afterwards, by which point the accumulation matrix is inconsistent. What it is never left in is a torn state: the kernel loads an entry's limbs, applies whole addends to them, and writes them all back together, so an addend is either absent or present in full.

If no surveys are provided, the accumulation matrix will compute one before beginning the accumulation.

The accumulation matrix only ever touches one column at a time, and accordingly, the caller can submit individual columns for accumulation. This interface is useful when matrices are large and columns are much smaller, so only one column and the accumulation matrix have to be resident. A positive or negative power-of-two scale factor also may be specified, and the shift is applied before adding the input into the limb-column.

A symmetric n × n matrix can be accumulated into one stored triangle. That halves the memory and the input traffic, and runs up to 2.2× faster on the CPU, but the deeper reason is exactness. In full storage, the two halves would hold equal values in different limbs under different column scales. Stored once, it cannot drift. The input is assumed to be symmetric.

Accumulation can be spread over worker threads, partitioned by column. Columns are independent in storage and each is touched by exactly one worker, so there is no locking on the hot path and the result is bit-identical to a serial run. On a 4096 × 64 accumulation over 8 cores, multithreading improves performance by >5x, from 0.80 to 4.31 Gelem/s. Scaling stops when the working set leaves cache: at 65536 × 64, the input alone is 33.6 MB against 32 MiB of L3, and 8 threads reach only 1.11 Gelem/s because every core is waiting on the same DRAM. A single core in cache is issue-bound, at 2.32 IPC with a 1% miss rate, which is why it scales at all until it doesn't.

To conserve bandwidth to the accumulation matrix, batches of input matrices may be accumulated in a single pass: the limbs are read once, every addend applied in registers, and written once. On the DRAM-bound case, eight matrices per pass increases performance from 1.14 to 3.01 Gelem/s and is still rising.

# Reading Back The Result

Once the accumulation matrix is computed, the library provides the following services:
- Extraction of its contents at `double` precision, correctly rounded: round-to-nearest, ties-to-even, over the full range. `to_double` returns one entry and `to_matrix` the whole thing. A sum beyond the range of a `double` reads as ±inf. One that lands in the denormal range is rounded directly to a multiple of 2⁻¹⁰⁷⁴, rather than twice, first to 53 bits and then again on the way into the denormals.
- Extraction of the *residual*, i.e. the portion of each matrix element that rounding dropped, itself correctly rounded. The subtraction is done in the matrix's own fixed point, so the residual is the true difference and not a floating point estimate of it. It is at most half an ulp of the returned value and is negative exactly when the rounding went up, so the pair is a non-overlapping two-term expansion of about 106 bits. It is +0.0 exactly when the value was representable, which `is_exactly_representable` reports on its own. The residual is optional: `to_double` and `to_matrix` each have an overload that takes somewhere to put it, `to_double(residual, i, j)` and `to_matrix_with_residual`, and the forms without it do no extra work.
- Extraction of an entry in decimal form by `to_exact_decimal`, every digit, never rounded. A binary fixed-point value always has a terminating decimal expansion, since V · 2⁻ᵏ = (V · 5ᵏ) / 10ᵏ, so the expansion is finite however wide the sum. A single 0.1 accumulated once reads back as 0.1000000000000000055511151231257827021181583404541015625.
- Extraction of the *mean* over a count the caller supplies, correctly rounded as a quotient. Dividing the rounded sum would round twice; `to_double_mean` and `to_matrix_mean` round once, by long division of the exact sum, and offer the same residual overloads as the readbacks above, `to_double_mean(residual, i, j, count)` and `to_matrix_mean_with_residual`. The count is the caller's because the accumulation matrix does not know what a sum stands for: a scaled add is a weighted one, and per-element, per-column and whole-matrix adds all feed the same cell.

# CUDA Considerations

When CUDA is available, developers can avail themselves of `CudaAccumulationMatrix`, a separate class that keeps the accumulation matrix in device memory.[^cuda] The kernels are bandwidth-bound throughout, at 97 to 99 percent of the card's streaming rate, so performance is determined by how many times an input crosses a bus and how often the CPU and GPU wait on each other.

The survey runs on whichever processor already holds the matrix. A producer holding a matrix in host memory surveys it on the CPU. A producer whose matrix never leaves the device, delivered by GPUDirect or computed there, calls `survey_matrix_col_major_device` instead, a kernel that writes the same 12 bytes per column into device memory or mapped host memory.

The two forms enlist CPU/GPU concurrency in different ways. When running the survey on the CPU, that computation can run during an asynchronous transfer of the matrix data from host to device memory. When running the survey on the GPU, writing the output data to mapped host memory requires the CPU to synchronize to avoid a race condition.

If the surveys of every input matrix are provided to the `CudaAccumulationMatrix` constructor, it sets each column's exponent and limb count once. No survey kernel is launched for any submission after that, and with no verdict to wait for, the device accumulate stops synchronizing altogether.

With the extents known before the accumulate launches, the kernel reads each input exactly once, so the input matrix can be read straight out of mapped memory at bus-saturating speeds rather than copying it to device memory first. Input from the host must therefore be page-locked and device-mapped. The kernel dereferences the pointer it is given, so pageable memory faults rather than merely running slowly; it is rejected on submission, where the diagnostic can say what to allocate instead.

The kernel is still streaming an input over PCIe after the call that submitted it has returned, so a producer that refills that buffer immediately corrupts about 63 percent of entries. An earlier design inferred the buffer's lifetime from the pointer's memory type, which gave one function two contracts and nothing in the signature to tell them apart. Every submission now returns an `InputRead`, a move-only handle owning a CUDA event recorded on the compute stream: `wait()` blocks until the device has finished with the buffer and `ready()` asks the same question without blocking, so a producer with two buffers in rotation fills one while the device reads the other and never has to synchronize the stream. Destroying a handle without waiting is not an error, but a statement that the buffer will not be written again.

The accumulation matrix uses a blocking stream, which implicitly synchronizes with the legacy null stream, so a producer that fills its buffers with plain `cudaMemcpy` serializes against the accumulate. At 65536 × 64 a two-buffer rotation costs 1776 µs a batch that way, against 1301 µs with the producer on a stream of its own, where the fill alone is 1256 µs. On its own stream the producer hides 96 percent of the accumulate.

PCIe caps host input at about 3.3 Gelem/s regardless of kernel quality. Device-resident input is 2.35 to 2.45× better, and submitting eight device-resident matrices in one batch reaches 25.3 Gelem/s, because the limb-columns are read and written once for the whole batch instead of once per matrix. That gap is the whole argument for delivering data over GPUDirect. Measured on a Ryzen 7 7700X and an RTX 3060, in exact accumulations per second:

| | 4096×64 | 16384×64 | 65536×64 |
| --- | --- | --- | --- |
| CPU, 1 thread | 0.80 | 0.77 | 0.67 |
| CPU, 8 threads | 4.31 | 4.65 | 1.11 |
| CPU, 8 threads, batched K=8 | | | 3.01 |
| GPU, host input | 3.11 | 3.26 | 3.32 |
| GPU, input resident | 4.40 | 6.07 | 8.13 |
| GPU, resident, batched K=8 | | | 25.3 |

Every claim above is pinned by a test rather than remembered. Exact accumulation restores associativity, so 24 random permutations of ten matrices spanning 400 binades must agree digit for digit; whole-matrix, column-major, per-column and per-element accumulation must all produce the same total; the threaded, batched, symmetric, pre-sized and GPU paths must each match the plain serial one limb for limb; and the AVX-512 kernels are cross-checked against the scalar ones by forcing the portable path. Catastrophic cancellation (10³⁰⁰ + 1 − 10³⁰⁰ yields exactly 1), a million values whose ulp is far below the running total, and every rounding tie from 1 + 2⁻⁵³ down to the 2⁻¹⁰⁷⁵ round-to-zero tie are checked individually. The suite runs in 25 ms and is clean under ASan, UBSan and ThreadSanitizer.


[^cuda]: `CudaAccumulationMatrix` is not an acceleration that gets switched on when CUDA is present. Swapping the CPU's scalar kernel for its AVX-512 one leaves the memory and the object identical, so that is a true backend swap and it is chosen at runtime. The CUDA limb arrays live in device memory, and putting them behind the same interface would disguise host-to-device transfers as ordinary calls. What the two implementations share is the algorithm: when to rescale, when to widen, how the exponent moves.

[^ij]: `i` is the row and `j` the column, in the usual order. The two are not interchangeable here: a column owns an exponent, its own limb allocations and exactly one worker thread, while a row is only an offset into each of those allocations.

[^limb]: The term appears to have been coined by the authors of the GNU Multiple Precision Arithmetic Library (GMP). Its manual defines a limb as "the part of a multi-precision number that fits in a single machine word", and explains the choice: "We chose this word because a limb of the human body is analogous to a digit, only larger, and containing several digits." GMP allows 32 or 64 bits; `truesum` uses 64 bits. [GNU MP 6.3.0 manual, §3.2, Nomenclature and Types](https://gmplib.org/manual/Nomenclature-and-Types)

[^presized]: None of them changes on a pre-sized accumulation matrix. When the surveys of every incoming matrix are supplied at construction, `exponent` and the limb count are set once from those extents, `max_addend_bits` and `add_count` are never touched at all, and each accumulate skips the survey, the rescale, and the widen.

[^three]: 1.0 is chosen for legibility, and it is the tidy case: a mantissa of 1 means the value is a power of two. Nothing about odd normalization reduces every mantissa to a single bit. 3.0 is exactly 3 × 2⁰, and its 53-bit form is 6755399441055744 × 2⁻⁵¹, so `frexp` implies a last place fifty-one positions below the value's own rather than fifty-two.
