`truesum`: Fast, Arbitrary-Precision Floating Point Summation

Floating point operations notoriously are not associative. In particular, (a+b)+c ≠ a+(b+c), introducing a degree of uncertainty about accuracy of workloads that rely on floating point summation.

To eliminate that uncertainty, our library `truesum` provides an *accumulation matrix* that can perform fast, elementwise exact summation of floating point matrices. Each column is represented as block floating point, with a shared exponent and the mantissas held in an SOA layout of 64-bit *limbs*[^limb], the segments of a multi-precision number that each fit in a single machine word.

A *survey* is a compact (12B per column) characterization of an input matrix that may be used to pre-size an accumulation matrix to receive its contents without overflow, rescaling, or the allocation of new column-limbs during the accumulation. Extending this idea, an accumulation matrix can be pre-sized for multiple input matrices by supplying one survey for each prospective input matrix.

A final feature of the accumulation matrix is that it can return the mean, exactly rounded, instead of requiring the client to divide a readback by the number of accumulations.

# Structure Of Arrays

The SOA layout enables both AVX-512 and CUDA implementations to process summations at the speed of memory bandwidth.

Within the matrix, entry `(i, j)`[^ij] is a two's complement integer `V` of several 64-bit limbs. Written the usual way (most significant first), it reads `V = [ limb 2 | limb 1 | limb 0 ]`. The SOA layout provides for a separate, contiguous allocation per limb position, so limbs 0, 1 and 2 for `V` are stored at the same offset into three different arrays. For want of a better term, we refer to these arrays as *limb-columns*.

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

On the assumption that columns contain similarly scaled values, the library organizes a block floating point scheme that shares an exponent per column: not the exponent of its largest value, but the position of the *least significant* bit of any of its values. This arrangement stands in contrast to IEEE standard representations whose exponents identify the position of the *most* significant set bit. Recall that floating point values have more than one exact decomposition into an integer times a power of two: 1.0 is exactly 1 × 2⁰, but just as exactly 2⁵² × 2⁻⁵².[^three] IEEE-754 stores the latter, because the most significant set bit is redundant: with the exception of denormals, the format encodes a hidden bit to increase its effective precision.

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

`max_addend_bits` is the greatest number of bits any addend in the column has spanned, counting from the column's current bit 0 up to that addend's most significant bit. `max_addend_bits` and `add_count`, a running count of the addends, both are updated from a survey of the incoming values, taken before any of those values is added. A rescale adjusts `max_addend_bits` too, adding the shift it applied, since moving the floor down by that many bits widens every addend's span from the new bit 0 by the same amount. Both reset to zero when the column is emptied, or when a rescale finds it already empty, since a column holding nothing constrains nothing.

`exponent` decreases whenever a value arrives whose least significant set bit sits below the column's current floor. Every entry is then left-shifted to match: existing bits move, and none are discarded.

An increase in the limb count, and allocation of one or more new limb columns, may be prompted by any of three events:

- A value whose most significant bit sits above the high water mark denoted by `max_addend_bits`,
- A lower exponent, since that prompts a rescale of `shift` bits, and
- More addends, since a longer run can sum to a larger total.

The sum of these considerations determines the number of bits needed for a column. `max_addend_bits` covers the largest single addend; `ceil_log2(add_count + 1)` bits are reserved above `max_addend_bits`, gaining one more each time the count doubles; and a final bit carries the sign. Rounding the total up to a whole number of limbs gives the allocation.

A batch is surveyed before any of it is applied: for any given column, a first pass over the incoming values determines whether column must rescaled and widened before any values are added.

Reallocation and rescaling of the column may be prompted by either very large or very small values.

Consider a column with 1.0, which fixes its exponent at 2⁰ and its allocation at a single limb. Adding a large value cannot change `exponent`, because a `double` of magnitude 2ᵏ has its lowest set bit no lower than 2ᵏ⁻⁵² and so never reaches below the floor, but it may increase the number of limbs (2⁵⁰⁰ takes the column to 8 limbs and 10³⁰⁰ to 16, with `exponent` unchanged at 2⁰). Adding a small value decreases `exponent` and increases the limb count as well (2⁻⁵⁰⁰ drops `exponent` to 2⁻⁵⁰⁰ and takes the column to 8 limbs; 10⁻³⁰⁰ drops it to 2⁻¹⁰⁴⁹ and 17 limbs). Both extremes cost limbs, by opposite routes: a large value raises the top of the window, a small one lowers the bottom.

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

A `double` arrives as a 53-bit odd mantissa placed at `shift`, the distance from its own true ulp down to the column's. Since 53<64, a `double` may span at most two limbs, so only `off` and `off + 1` need to be updated.

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

A column need only allocate as many column-limbs as needed for the dynamic range spanned by its own data. In the bundled demo, a 4×4 accumulation matrix holding values from 10⁻²³ to 10²⁶ occupies 584 bytes, with per-column allocations of two and three limbs rather than the thirty-three that covering the whole `double` range at once would take.

Allocating only for the range spanned by a column's data spans implies the range must be known before undertaking any accumulations. Surveys supply the needed per-column metadata: the extents occupied by the incoming values, which the accumulation matrix can use t determine whether the column must be rescaled and whether it needs more limb columns.

Considering surveys separately from the matrix they describe enables the work to be done by different constituencies. In a data center, the nodes that computed the matrices can also compute their surveys and transmit them ahead to the accumulation node; once the accumulation matrix there has been pre-sized to receive them all, the matrices themselves follow, one matrix or even one column at a time. The final output is bit-identical regardless of they order they arrive or are processed.

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

A pre-sized accumulation matrix trusts its surveys and verifies them as they are processed. Each kernel tests, from values already in registers, whether the inputs fit the column as sized: three tests on an addend before it is applied, and a fourth on the running total as it crosses into the top limb.

Between them, those tests detect any deficiencies in the surveys. A survey can understate a column's floor or its ceiling, or fail to report a non-finite value; separately, a caller can supply fewer surveys than the matrices it goes on to submit. Each of those throws `std::runtime_error` rather than quietly returning an incorrect sum. They err toward alarm rather than silence: a total that wraps past the sign bit and later wraps back is correct by modular arithmetic, but is reported anyway.

An offending addend is skipped, an overflowed sum is left where it wrapped, the rest of the batch is applied, and the error is raised afterwards, by which point the accumulation matrix is inconsistent and the message says so. A column that surveyed itself cannot trip any of the tests, since its own rescale and widen have already happened; there they cost nothing and assert that the survey and the accumulation agree.

If no surveys are provided, the accumulation matrix will compute one  before beginning the accumulation.

Once the accumulation matrix is computed, the library provides the following services:
- Extraction of its contents at `double` precision, correctly rounded: round-to-nearest, ties-to-even, over the full range. A sum beyond the range of a `double` reads as ±inf. One that lands in the denormal range is rounded directly to a multiple of 2⁻¹⁰⁷⁴, rather than twice, first to 53 bits and then again on the way into the denormals.
- Extraction of the *residual*, i.e. the portion of each matrix element that rounding dropped, itself correctly rounded. The subtraction is done in the matrix's own fixed point, so the residual is the true difference and not a floating point estimate of it. It is at most half an ulp of the returned value and is negative exactly when the rounding went up, so the pair is a non-overlapping two-term expansion carrying about 106 bits. It is +0.0 exactly when the value was representable, which a separate query reports.
- Extraction of the matrix in decimal form, every digit, never rounded. A binary fixed-point value always has a terminating decimal expansion, since V · 2⁻ᵏ = (V · 5ᵏ) / 10ᵏ, so the expansion is finite however wide the sum. A single 0.1 accumulated once reads back as 0.1000000000000000055511151231257827021181583404541015625.
- Extraction of the *mean* over a count the caller supplies, correctly rounded as a quotient. Dividing the rounded sum would round twice; this rounds once, by long division of the exact sum, and returns the same kind of residual. The count is the caller's because the accumulation matrix does not know what a sum stands for: a scaled add is a weighted one, and per-element, per-column and whole-matrix adds all feed the same cell.

The accumulation matrix only ever touches one column at a time, and accordingly, the caller can submit individual columns for accumulation. This interface is useful when matrices are large and columns are much smaller, so only one column and the accumulation matrix have to be resident. Inputs can also be added scaled by any power of two, positive or negative, which is the only scaling that stays exact without widening the mantissa, and subtracted, so an accumulation matrix can be brought back to exact zero and reused with each column's learned scale intact.

A symmetric n × n matrix can be accumulated into one stored triangle. That halves the memory and the input traffic, and runs up to 2.2× faster on the CPU, but the deeper reason is exactness. In full storage, the two halves would hold equal values in different limbs under different column scales. Stored once, it cannot drift. The input is taken to be symmetric and that is not checked.

Accumulation can be spread over worker threads, partitioned by column. Columns are independent in storage and each is touched by exactly one worker, so there is no locking on the hot path and the result is bit-identical to a serial run. On a 4096 × 64 accumulation over 8 cores, multithreading improves performance by >5x, from 0.80 to 4.31 Gelem/s. Scaling stops when  the working set leaves cache: at 65536 × 64, the input alone is 33.6 MB against 32 MiB of L3, and 8 threads reach only 1.11 Gelem/s because every core is waiting on the same DRAM. A single core in cache is issue-bound, at 2.32 IPC with a 1% miss rate, which is why it scales at all until it doesn't.

To conserve bandwidth to the accumulation matrix, batches of input matrices may be accumulated in a single pass: the limbs are read once, every addend applied in registers, and written once. On the DRAM-bound case, eight matrices per pass increases performance from 1.14 to 3.01 Gelem/s and is still rising.

# CUDA Considerations

When CUDA is available, developers can avail themselves of `CudaAccumulationMatrix`, a separate class that keeps the accumulation matrix in device memory.[^cuda] The kernels are bandwidth-bound throughout, at 97 to 99 percent of the card's streaming rate, so performance is determined by how many times an input crosses a bus and how often the CPU and GPU wait on each other.

A producer holding a matrix in host memory surveys it on the CPU. A producer whose matrix never leaves the device, delivered by GPUDirect or computed there, calls `survey_matrix_col_major_device` instead, a kernel that writes the same 12 bytes a column into device memory, or into mapped host memory if that is where the caller wants to read it. The survey runs on whichever processor already holds the matrix.

The two forms enlist CPU/GPU concurrency in different ways. When running the survey on the CPU, that computation can run during an asynchronous transfer of the matrix data from host to device memory. When running the survey on the GPU, writing the output data to mapped host memory requires the CPU to synchronize to avoid a race condition.

Given the surveys of every matrix that will arrive, the constructor sets each column's exponent and limb count once, no survey kernel is launched, and the device accumulate stops synchronizing altogether.

With the extents settled before the accumulate launches, the kernel reads each input exactly once, and that is what lets it stream host input straight out of mapped memory instead of copying it to the device first. At 65536 × 64, copying the 33.6 MB and then reading it from device memory takes 1768 µs against 1257 µs to read it in place, and 1257 µs is what the transfer alone costs, so the arithmetic hides entirely behind the bus. Input from the host must therefore be page-locked and device-mapped. The kernel dereferences the pointer it is given, so pageable memory faults rather than merely running slowly; it is rejected on submission, where the diagnostic can name what to allocate instead.

The kernel is still streaming an input over PCIe after the call that submitted it has returned, so a producer that refills that buffer immediately corrupts about 63 percent of entries. An earlier design inferred the buffer's lifetime from the pointer's memory type, which gave one function two contracts and nothing in the signature to tell them apart. Every submission now returns an `InputRead`, a move-only handle owning a CUDA event recorded on the compute stream: `wait()` blocks until the device has finished with the buffer and `ready()` asks the same question without blocking, so a producer with two buffers in rotation fills one while the device reads the other and never has to synchronize the stream. Destroying a handle without waiting is not an error, but a statement that the buffer will not be written again.

Fill those buffers on your own stream. The accumulation matrix uses a blocking stream, which implicitly synchronizes with the legacy null stream, so a producer filling its buffers with plain `cudaMemcpy` serializes against the accumulate and gets no overlap at all. The trap is CUDA's rather than `truesum`'s, but the rotation is where it bites. At 65536 × 64 a two-buffer rotation costs 1776 µs a batch that way, against 1301 µs with the producer on a stream of its own and 1256 µs for the fill by itself. The accumulate is 96 percent hidden when the streams are kept apart and not hidden at all when they are not.

PCIe caps host input at about 3.3 Gelem/s regardless of kernel quality. Device-resident input is 2.35 to 2.45× better, and folded eight deep reaches 25.3 Gelem/s; that gap is the whole argument for delivering data over GPUDirect. Measured on a Ryzen 7 7700X and an RTX 3060, in exact accumulations per second:

| | 4096×64 | 16384×64 | 65536×64 |
| --- | --- | --- | --- |
| CPU, 1 thread | 0.80 | 0.77 | 0.67 |
| CPU, 8 threads | 4.31 | 4.65 | 1.11 |
| CPU, 8 threads, folded K=8 | | | 3.01 |
| GPU, host input | 3.11 | 3.26 | 3.32 |
| GPU, input resident | 4.40 | 6.07 | 8.13 |
| GPU, resident, folded K=8 | | | 25.3 |

Every claim above is pinned by a test rather than remembered. Exact accumulation restores associativity, so 24 random permutations of ten matrices spanning 400 binades must agree digit for digit; whole-matrix, column-major, per-column and per-element accumulation must all produce the same total; the threaded, folded, symmetric, pre-sized and GPU paths must each match the plain serial one limb for limb; and the AVX-512 kernels are cross-checked against the scalar ones by forcing the portable path. Catastrophic cancellation (10³⁰⁰ + 1 − 10³⁰⁰ yields exactly 1), a million values whose ulp is far below the running total, and every rounding tie from 1 + 2⁻⁵³ down to the 2⁻¹⁰⁷⁵ round-to-zero tie are checked individually. The suite runs in 25 ms and is clean under ASan, UBSan and ThreadSanitizer.


[^cuda]: `CudaAccumulationMatrix` is not an acceleration that gets switched on when CUDA is present. Swapping the CPU's scalar kernel for its AVX-512 one leaves the memory and the object identical, so that is a true backend swap and it is chosen at runtime. The CUDA limb arrays live in device memory, and putting them behind the same interface would disguise host-to-device transfers as ordinary calls. What the two implementations share is the algorithm: when to rescale, when to widen, how the exponent moves.

[^ij]: `i` is the row and `j` the column, in the usual order. The two are not interchangeable here: a column owns an exponent, its own limb allocations and exactly one worker thread, while a row is only an offset into each of those allocations.

[^limb]: The term is GMP's. Its manual defines a limb as "the part of a multi-precision number that fits in a single machine word", and explains the choice: "We chose this word because a limb of the human body is analogous to a digit, only larger, and containing several digits." GMP allows 32 or 64 bits; `truesum` uses 64 bits. [GNU MP 6.3.0 manual, §3.2, Nomenclature and Types](https://gmplib.org/manual/Nomenclature-and-Types)

[^presized]: None of them changes on a pre-sized accumulation matrix. When the surveys of every incoming matrix are supplied at construction, `exponent` and the limb count are set once from those extents, `max_addend_bits` and `add_count` are never touched at all, and each accumulate skips the survey, the rescale, and the widen.

[^three]: 1.0 is chosen for legibility, and it is the tidy case: a mantissa of 1 means the value is a power of two. Nothing about odd normalization reduces every mantissa to a single bit. 3.0 is exactly 3 × 2⁰, and its 53-bit form is 6755399441055744 × 2⁻⁵¹, so `frexp` implies a last place fifty-one positions below the value's own rather than fifty-two.
