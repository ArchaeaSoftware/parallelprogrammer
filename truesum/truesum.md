`truesum`: Fast, Arbitrary-Precision Floating Point Summation

Floating point operations notoriously are not associative. In particular, (a+b)+c ≠ a+(b+c), introducing a degree of uncertainty about accuracy of workloads that rely on floating point summation.

To eliminate that uncertainty, our library `truesum` provides an *accumulation matrix* that can perform fast, elementwise exact summation of floating point matrices. Each column is represented as block floating point, with a shared exponent and the mantissas held in an SOA layout of 64-bit *limbs*[^limb], the segments of a multi-precision number that each fit in a single machine word.

A *survey* is a compact (12B per column) characterization of an input matrix that may be used to pre-size an accumulation matrix to receive its contents without overflow, rescaling, or the allocation of new column-limbs during the accumulation. Extending this idea, an accumulation matrix can be pre-sized for multiple input matrices by supplying one survey for each prospective input matrix. Separating surveys from their matrices makes distributed accumulation practical: the nodes that computed the partial results also compute the surveys, the surveys travel ahead to the accumulation node and, once the accumulation matrix has been pre-sized to receive them all, the matrices may be submitted and accumulated in any order.

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

That transposition enables data-parallel computation across rows. Propagating carries along the limbs of a single entry is inherently sequential: limb `k + 1` cannot be computed until after limb `k`. Rows, by contrast, are completely independent of one another, so the limb-major memory layout enables rows to be processed by a single vector instruction (8 64-bit elements at a time on AVX-512), while the sequential carry chain becomes the loop over limb positions, one vector instruction per position.

# Per-Column Block Floating Point (Shared Exponents)

On the assumption that columns contain similarly scaled values, the library organizes a block floating point scheme that shares an exponent per column: not the exponent of its largest value, but the position of the *least significant* bit of any of its values. This arrangement stands in contrast to IEEE standard representations whose exponents identify the position of the *most* significant set bit. Recall that floating point values have more than one exact decomposition into an integer times a power of two: 1.0 is exactly 1 × 2⁰, but just as exactly 2⁵² × 2⁻⁵².[^three] IEEE-754 stores the latter, because the most significant set bit is redundant: with the exception of denormals, the format encodes a hidden bit to increase its effective precision.

The power of two specified by a column's exponent, then, is each value's true unit in the last place, or ulp: the significance of its lowest bit.Taking the minimum over true ulps rather than over standard exponents can save up to 52 bits per column: a column of small integers stays at exponent 2⁰ and one limb wide, whereas a `frexp`-based scale would reduce it to 2⁻⁵², forcing the allocation of a second limb holding nothing but zeros.

```
       limb 2                  limb 1                  limb 0
  ┌───────────────────────┬───────────────────────┬───────────────────────┐
  │S                      │                       │                       │
  └───────────────────────┴───────────────────────┴───────────────────────┘
   191                 128 127                  64 63                    0
   2^(e+191)                                                           2^e
```

A column maintains four parameters that must be updated with a survey before applying an accumulation[^presized]: `exponent`, `max_addend_bits`, `add_count`, and the number of limbs. `exponent` can only decrease; the other three only increase.

`max_addend_bits` is the greatest number of bits any addend in the column has spanned, counting from the column's current bit 0 up to that addend's most significant bit. `max_addend_bits` and `add_count`, a running count of the addends, both are updated from a pass over the incoming values, taken before any of those values is added. They reset to zero when the column is emptied, or when a rescale finds it already empty, since a column holding nothing constrains nothing.

`exponent` decreases whenever a value arrives whose least significant set bit sits below the column's current floor. Every entry is then left-shifted to match: existing bits move, and none are discarded.

An increase in the limb count, and allocation of one or more new limb columns, may be prompted by any of three events:

- A value whose most significant bit sits above the high water mark denoted by `max_addend_bits`,
- A lower exponent, since that prompts a rescale of `shift` bits, and
- More addends, since a longer run can sum to a larger total.

Those causes feed three quantities, and their sum is the number of bits a column needs. `max_addend_bits` covers the largest single addend. The headroom covers the sum outgrowing that addend: `add_count` of them can total `add_count` times the largest, so `ceil_log2(add_count + 1)` bits are reserved above `max_addend_bits`, gaining one more each time the count doubles. A final bit carries the sign. Rounding the total up to a whole number of limbs gives the allocation.

A batch is surveyed before any of it is applied: for any given column, a first pass over the incoming values determines whether column must rescaled and widened before any values are added.

Reallocation and rescaling of the column may be prompted by either very large or very small values.

Consider a column with 1.0, which fixes its exponent at 2⁰ and its allocation at a single limb. Adding a large value cannot change `exponent`, because a `double` of magnitude 2ᵏ has its lowest set bit no lower than 2ᵏ⁻⁵² and so never reaches below the floor, but it may increase the number of limbs (2⁵⁰⁰ takes the column to 8 limbs, 10³⁰⁰ to 16, both still at 2⁰). Adding a small value decreases `exponent` and increases the limb count as well (2⁻⁵⁰⁰ gives 2⁻⁵⁰⁰ and 8 limbs, 10⁻³⁰⁰ gives 2⁻¹⁰⁴⁹ and 17). Both extremes cost limbs, by opposite routes: a large value raises the top of the window, a small one lowers the bottom.

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

The carry chain runs fully vectorized, without an overflow test until it reaches the most-significant limb. At the top limb, a sign comparison of the operands is performed against the result. The chain reaches that limb only when an addend lands in it or a carry propagates to it, so the test runs exactly where overflow is possible and nowhere else.

A column pays only for the dynamic range its own data spans. In the bundled demo, a 4×4 accumulation matrix holding values from 10⁻²³ to 10²⁶ occupies 584 bytes, with per-column allocations of two and three limbs rather than the thirty-three that covering the whole `double` range at once would take.

To make that guarantee before anything is added, the library computes a *survey* of the input matrix to predetermine whether any columns will have to be re-scaled and any new limbs must be allocated to hold the output. One benefit of this separation of concerns is that if the input matrices are available for scanning before summation begins, providing all the input surveys at initialization time enables a one-time allocation of all memory needed to hold the final summation. This feature is particularly useful in data center applications where the input matrices may be computed on different nodes: the same nodes that computed the input matrices can pre-scan those matrices on behalf of the accumulation node. If the compact surveys are transferred to the accumulation node before summation begins, the matrices themselves can be transferred to the accumulation node for processing, one at a time, in any order; and for any order, the final output will be bit-identical.

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

A pre-sized accumulation matrix takes its surveys on trust: only the number of matrices is checked. A matrix that reaches outside its declared extents is not prevented, because preventing it would mean surveying, which is what pre-sizing removed. What the kernels do check, from values already in registers and at 0.8 to 4.8 percent of their running time, is whether each incoming value fits the column as sized: a value that is not finite, a value whose true ulp lies below the column's exponent and would lose a bit off the bottom, a value whose top reaches the column's sign bit, and, at the top limb, a sum that has crossed it. The first three are addends the column was never sized for and are found before they are applied; the fourth is the sum itself, found where it happens, and it is the case a survey understated by a few bits produces: every addend fits and the total does not. All are detected, not prevented. An offending addend is skipped, an overflowed sum is left where it wrapped, the rest of the batch is applied, and the error is raised afterward, by which point the accumulation matrix is inconsistent and the message says so. The checks are always on, but an accumulation matrix that surveyed its own column cannot trip them, since the rescale and widen have already happened; there they are a free assertion that the survey and the accumulate agree, and only a wrong survey from a producer can make them fire. A survey also records whether a column held an infinity or a NaN, which the accumulation matrix rejects outright rather than letting either poison a column.

Without surveys the accumulation matrix takes its own. Accumulating a column is then two passes over its values: the first learns only the exponent range, because the rescale and widen decisions have to be made before anything can be added; the second fuses decomposition into the add, so no decomposed form is ever written to memory. Either way, results are bit-identical, and the test suite asserts it.

Once the accumulation matrix is computed, the library provides the following services:
- Extraction of its contents at `double` precision, correctly rounded: round-to-nearest, ties-to-even, over the full range. A sum beyond the range of a `double` reads as ±inf. One that lands in the denormal range is rounded directly to a multiple of 2⁻¹⁰⁷⁴, rather than twice, first to 53 bits and then again on the way into the denormals.
- Extraction of the *residual*, i.e. the portion of each matrix element that rounding dropped, itself correctly rounded. The subtraction is done in the matrix's own fixed point, so the residual is the true difference and not a floating point estimate of it. It is at most half an ulp of the returned value and is negative exactly when the rounding went up, so the pair is a non-overlapping two-term expansion carrying about 106 bits. It is +0.0 exactly when the value was representable, which a separate query reports.
- Extraction of the matrix in decimal form, every digit, never rounded. A binary fixed-point value always has a terminating decimal expansion, since V · 2⁻ᵏ = (V · 5ᵏ) / 10ᵏ, so the expansion is finite however wide the sum. A single 0.1 accumulated once reads back as 0.1000000000000000055511151231257827021181583404541015625.
- Extraction of the *mean* over a count the caller supplies, correctly rounded as a quotient. Dividing the rounded sum would round twice; this rounds once, by long division of the exact sum, and returns the same kind of residual. The count is the caller's because the accumulation matrix does not know what a sum stands for: a scaled add is a weighted one, and per-element, per-column and whole-matrix adds all feed the same cell.

The accumulation matrix only ever touches one column at a time, and the caller can feed it the same way. A single column can be accumulated from `rows` contiguous values, so only the accumulation matrix has to be resident, not any input matrix. That is what lets it be far larger than any input that could sit beside it: on a 46 GB machine, a square accumulation matrix goes from 42,000 × 42,000 to 60,000 × 60,000 simply by feeding it a column at a time. Inputs can also be added scaled by any power of two, positive or negative, which is the only scaling that stays exact without widening the mantissa, and subtracted, so an accumulation matrix can be brought back to exact zero and reused with each column's learned scale intact.

A symmetric n × n matrix can be accumulated into one stored triangle. That halves the memory and the input traffic, and is worth 2.0 to 2.2× on the CPU, but the deeper reason is exactness. In full storage the two halves would hold equal values in different limbs under different column scales, and symmetry would become a property the data had to keep earning. Stored once, it cannot drift. The input is taken to be symmetric and that is not checked, because checking it means reading the half the design exists to avoid reading.

Accumulation can be spread over worker threads, partitioned by column. Columns are independent in storage and each is touched by exactly one worker, so there is no locking on the hot path and the result is bit-identical to a serial run. On a 4096 × 64 accumulation over eight cores this takes 0.80 to 4.31 Gelem/s. Scaling stops where the working set leaves cache: at 65536 × 64 the input alone is 33.6 MB against 32 MiB of L3, and eight threads reach only 1.11 Gelem/s because every core is waiting on the same DRAM. A single core in cache is issue-bound, at 2.32 IPC with a 1% miss rate, which is why it scales at all until it doesn't.

Past that point the lever is traffic, not threads. Several input matrices can be folded into one pass over the accumulation matrix: the limbs are read once, every addend applied in registers, and written once. The saving is in traffic to the accumulation matrix only, so it shows up where that traffic is the bound. On the DRAM-bound case, eight matrices per pass takes 1.14 to 3.01 Gelem/s and is still rising, while single-threaded and in cache the extra input streams cost more than they save.

The CUDA implementation is a separate container with device-resident storage, not a third kernel behind the CPU dispatch, and the GPU is bandwidth-bound throughout, at 97 to 99 percent of the card's streaming rate. Because the survey runs on the host, input is read exactly once and never re-read, and that is what lets the kernel stream it straight out of mapped host memory instead of copying it to the device first. At 65536 × 64, copying the 33.6 MB and then reading it from device memory takes 1768 µs against 1257 µs to read it in place, and 1257 µs is what the transfer alone costs, so the arithmetic hides entirely behind the bus. PCIe therefore caps host input at about 3.3 Gelem/s regardless of kernel quality. Device-resident input is 2.35 to 2.45× better, and folded eight deep reaches 25.3 Gelem/s; that gap is the whole argument for delivering data over GPUDirect. Measured on a Ryzen 7 7700X and an RTX 3060, in exact accumulations per second:

| | 4096×64 | 16384×64 | 65536×64 |
| --- | --- | --- | --- |
| CPU, 1 thread | 0.80 | 0.77 | 0.67 |
| CPU, 8 threads | 4.31 | 4.65 | 1.11 |
| CPU, 8 threads, folded K=8 | | | 3.01 |
| GPU, host input | 3.11 | 3.26 | 3.32 |
| GPU, input resident | 4.40 | 6.07 | 8.13 |
| GPU, resident, folded K=8 | | | 25.3 |

Every claim above is pinned by a test rather than remembered. Exact accumulation restores associativity, so 24 random permutations of ten matrices spanning 400 binades must agree digit for digit; whole-matrix, column-major, per-column and per-element accumulation must all produce the same total; the threaded, folded, symmetric, pre-sized and GPU paths must each match the plain serial one limb for limb; and the AVX-512 kernels are cross-checked against the scalar ones by forcing the portable path. Catastrophic cancellation (10³⁰⁰ + 1 − 10³⁰⁰ yields exactly 1), a million values whose ulp is far below the running total, and every rounding tie from 1 + 2⁻⁵³ down to the 2⁻¹⁰⁷⁵ round-to-zero tie are checked individually. The suite runs in 25 ms and is clean under ASan, UBSan and ThreadSanitizer.


[^ij]: `i` is the row and `j` the column, in the usual order. The two are not interchangeable here: a column owns an exponent, its own limb allocations and exactly one worker thread, while a row is only an offset into each of those allocations.

[^limb]: The term is GMP's. Its manual defines a limb as "the part of a multi-precision number that fits in a single machine word", and explains the choice: "We chose this word because a limb of the human body is analogous to a digit, only larger, and containing several digits." GMP allows 32 or 64 bits; `truesum` uses 64 bits. [GNU MP 6.3.0 manual, §3.2, Nomenclature and Types](https://gmplib.org/manual/Nomenclature-and-Types)

[^presized]: None of them changes on a pre-sized accumulation matrix. When the surveys of every incoming matrix are supplied at construction, `exponent` and the limb count are set once from those extents, `max_addend_bits` and `add_count` are never touched at all, and each accumulate skips the survey, the rescale, and the widen.

[^three]: 1.0 is chosen for legibility, and it is the tidy case: a mantissa of 1 means the value is a power of two. Nothing about odd normalization reduces every mantissa to a single bit. 3.0 is exactly 3 × 2⁰, and its 53-bit form is 6755399441055744 × 2⁻⁵¹, so `frexp` implies a last place fifty-one positions below the value's own rather than fifty-two.
