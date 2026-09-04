`truesum`: Fast, Arbitrary-Precision Floating Point Summation

Floating point operations notoriously are not associative. In particular, (a+b)+c ≠ a+(b+c), introducing a degree of uncertainty about accuracy of workloads that rely on floating point summation.

To eliminate that uncertainty, a library `truesum` provides fast, elementwise exact summation of floating point matrices using a block floating point representation that reserves one exponent per column. Accumulations are computed with no loss of intermediate precision. Each column is represented as block floating point, with a shared exponent and the mantissas held in 64-bit limbs with an SOA layout. 

# Structure Of Arrays

The SOA layout enables both AVX-512 and CUDA implementations to process summations at the speed of memory bandwidth.

Within the matrix, entry `(i, j)`[^ij] is a two's complement integer `V` of several 64-bit limbs. Written the usual way (most significant first), it reads `V = [ limb 2 | limb 1 | limb 0 ]`. The SOA layout provides for a separate, contiguous allocation per limb *position*, so limbs 0, 1 and 2 for `V` are stored at the same offset into three different arrays.

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

                └─ one vector instruction spans this way ─┘
```

That transposition enables data-parallel computation across rows. Propagating carries along the limbs of a single entry is inherently sequential: limb `k + 1` cannot be computed until after limb `k`. Rows, by contrast, are completely independent of one another, so the limb-major memory layout enables rows to be processed by a single vector instruction (8 64-bit elements at a time on AVX-512), while the sequential carry chain becomes the loop over limb positions, one vector instruction per position.

# Per-Column Block Floating Point (Shared Exponents)

On the assumption that columns contain similarly scaled values, the library organizes a block floating point scheme that shares an exponent per column.

A column's exponent is not the exponent of its largest value, but the weight of the least significant bit any of its values reaches. The distinction rests on a `double` having more than one exact decomposition into an integer times a power of two: 1.0 is exactly 1 × 2⁰, but just as exactly 2⁵² × 2⁻⁵². IEEE-754 stores the latter, because normalizing every significand until its leading bit is set makes that bit redundant: the format leaves it out and spends the space on precision. Denormals are the exception, having no leading 1 to hide. What the hidden bit costs is that the exponent must then track the leading bit, and that is the exponent `frexp` reports, so the last place it implies for 1.0 sits at 2⁻⁵², fifty-two positions below anything the value actually holds.

An accumulator has no fixed-width significand to economize on, so it normalizes from the other end. Every incoming mantissa is shifted right until it is odd, which pins the exponent to the trailing bit rather than the leading one, and that power of two is the value's true unit in the last place, or ulp: the weight of its lowest set bit, and the smallest step the value is really built from. Nothing is lost in the move, since the surplus only shifts between mantissa and exponent.

A column must sit at or below the true ulp of everything it will hold, because it stores integers scaled by that one power of two and anything finer has nowhere to go. Taking the minimum over true ulps rather than over `frexp` exponents can save up to 52 bits of width per column. A column of small integers stays at exponent 2⁰ and one limb wide, where a `frexp`-based scale would immediately sink it to 2⁻⁵², causing the allocation of a second limb holding nothing but zeros.

```
       limb 2                  limb 1                  limb 0
  ┌───────────────────────┬───────────────────────┬───────────────────────┐
  │S                      │                       │                       │
  └───────────────────────┴───────────────────────┴───────────────────────┘
   191                 128 127                  64 63                    0
   2^(e+191)                                                           2^e
```

A column carries two parameters: its exponent and its width, which can only decrease and increase, respectively. When the exponent is decreased, every entry in the column is left-shifted to match, so existing bits move and none are discarded.

The width is counted in bits and only grows, and it comes from a bound rather than from the data. No entry can exceed `count · 2^max_addend_bits`, so that bound plus one bit for the sign is the width the column needs, rounded up to whole 64-bit limbs before anything is allocated. The figure is settled before a batch is applied and is never read back from what is stored.

The two parameters answer to different kinds of value, and it is worth walking through both. Start a column with 1.0, which fixes its exponent at 2⁰ and its width at a single limb, then add something enormous. Adding 2⁵⁰⁰ leaves the exponent exactly where it was and takes the column to eight limbs; adding 10³⁰⁰ on top of that takes it to sixteen, still at 2⁰. A large value cannot move the floor, because a `double` of magnitude 2ᵏ carries at most 53 significant bits and so has its lowest set bit no lower than 2ᵏ⁻⁵². The true ulp of 10³⁰⁰ is 2⁹⁴⁶, far above the floor. All a large value can do is raise the ceiling.

Small values are what move the floor. Into that same column, 2⁻⁵⁰⁰ drops the exponent to 2⁻⁵⁰⁰ and the width to eight limbs, and 10⁻³⁰⁰ drops it to 2⁻¹⁰⁴⁹ and seventeen. Both extremes cost limbs, then, but by opposite routes: a large value raises the top of the window and a small one lowers the bottom. Lowering the bottom widens the column by exactly as much as the floor fell, since a rescale of `shift` bits adds `shift` to `max_addend_bits` as well.

```
  Widen: append a limb position. Nothing already stored moves.

      before   [ limb 1 | limb 0 ]
      after    [ limb 2 | limb 1 | limb 0 ]
                  new    └── unchanged ──┘
                         └─ not read, not written, not moved

  Rescale: the window slides down, and every entry shifts left to match.

      before   [ 1 0 1 1 ]                 bit 0 weighs 2^0
      after    [ 1 0 1 1 0 0 0 ]           bit 0 weighs 2^-3
                         └── three new low bits. The stored 1 0 1 1 did
                             not change; only what its bit 0 is worth did,
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

  width  =  max_addend_bits  +  ceil_log2(count + 1)  +  1
            └ what one addend    └ headroom, so that     └ the sign bit S,
              spans, lowest bit     `count` of them        which no sum may
              to highest            cannot carry past      ever reach
```

The carry chain therefore runs without an overflow test at any limb, which is what lets it vectorize. The one test it carries is at the top limb, a sign comparison of the operands against the result. The chain reaches that limb only when an addend lands in it or a carry climbs to it, so the test runs exactly where overflow is possible and nowhere else.

A column pays only for the dynamic range its own data spans. In the bundled demo, a 4×4 accumulator holding values from 10⁻²³ to 10²⁶ occupies 584 bytes, with per-column widths of 128 to 192 bits rather than the 2100 that covering the whole `double` range at once would take.

To make that guarantee before anything is added, the library computes a *survey* of the input matrix to predetermine whether any columns will have to be re-scaled and any new limbs must be allocated to hold the output. One benefit of this separation of concerns is that if the input matrices are available for scanning before summation begins, providing all the input surveys at initialization time enables a one-time allocation of all memory needed to hold the final summation. This feature is particularly useful in data center applications where the input matrices may be computed on different nodes: the same nodes that computed the input matrices can pre-scan those matrices on behalf of the accumulation node. If the compact surveys are transferred to the accumulation node before summation begins, the matrices themselves can be transferred to the accumulation node for processing, one at a time, in any order; and for any order, the final output will be bit-identical.

```
  Ahead of the data, one survey per column -- 12 bytes each, whatever the
  row count:

      node 1  ── survey ──┐
      node 2  ── survey ──┤
        ...               ├──>  accumulator sized once, exactly:
      node N  ── survey ──┘        no rescale, no widen, no survey, ever

  Then the matrices themselves, 8 bytes an element, in any order at all:

      node k  ══ matrix ═══════════════════>  accumulate  ══>  the same
      node 1  ══ matrix ═══════════════════>  accumulate  ══>  bits either
      node 3  ══ matrix ═══════════════════>  accumulate  ══>  way
```

A pre-sized accumulator takes its surveys on trust: only the number of matrices is checked. A matrix that reaches outside its declared extents is not prevented, because preventing it would mean surveying, which is what pre-sizing removed. What the kernels do check, from values already in registers and at 0.8 to 4.8 percent of their running time, is whether each incoming value fits the column as sized: a value that is not finite, a value whose true ulp lies below the column's exponent and would lose a bit off the bottom, a value whose top reaches the column's sign bit, and, at the top limb, a sum that has crossed it. The first three are addends the column was never sized for and are found before they are applied; the fourth is the sum itself, found where it happens, and it is the case a survey understated by a few bits produces: every addend fits and the total does not. All are detected, not prevented. An offending addend is skipped, an overflowed sum is left where it wrapped, the rest of the batch is applied, and the error is raised afterward, by which point the accumulator is inconsistent and the message says so. The checks are always on, but an accumulator that surveyed its own column cannot trip them, since the rescale and widen have already happened; there they are a free assertion that the survey and the accumulate agree, and only a wrong survey from a producer can make them fire. A survey also records whether a column held an infinity or a NaN, which the accumulator rejects outright rather than letting either poison a column.

Without surveys the accumulator takes its own. Accumulating a column is then two passes over its values: the first learns only the exponent range, because the rescale and widen decisions have to be made before anything can be added; the second fuses decomposition into the add, so no decomposed form is ever written to memory. Either way, results are bit-identical, and the test suite asserts it.

Once the accumulation matrix is computed, the library provides the following services:
- Extraction of its contents at `double` precision, correctly rounded: round-to-nearest, ties-to-even, over the full range. A sum beyond the range of a `double` reads as ±inf. One that lands in the denormal range is rounded directly to a multiple of 2⁻¹⁰⁷⁴, rather than twice, first to 53 bits and then again on the way into the denormals.
- Extraction of the *residual*, i.e. the portion of each matrix element that rounding dropped, itself correctly rounded. The subtraction is done in the accumulator's own fixed point, so the residual is the true difference and not a floating point estimate of it. It is at most half an ulp of the returned value and is negative exactly when the rounding went up, so the pair is a non-overlapping two-term expansion carrying about 106 bits. It is +0.0 exactly when the value was representable, which a separate query reports.
- Extraction of the matrix in decimal form, every digit, never rounded. A binary fixed-point value always has a terminating decimal expansion, since V · 2⁻ᵏ = (V · 5ᵏ) / 10ᵏ, so the expansion is finite however wide the sum. A single 0.1 accumulated once reads back as 0.1000000000000000055511151231257827021181583404541015625.
- Extraction of the *mean* over a count the caller supplies, correctly rounded as a quotient. Dividing the rounded sum would round twice; this rounds once, by long division of the exact sum, and returns the same kind of residual. The count is the caller's because the accumulator does not know what a sum stands for: a scaled add is a weighted one, and per-element, per-column and whole-matrix adds all feed the same cell.

The accumulator only ever touches one column at a time, and the caller can feed it the same way. A single column can be accumulated from `rows` contiguous values, so only the accumulator has to be resident, not any input matrix. That is what lets it be far larger than any input that could sit beside it: on a 46 GB machine, a square accumulator goes from 42,000 × 42,000 to 60,000 × 60,000 simply by feeding it a column at a time. Inputs can also be added scaled by any power of two, positive or negative, which is the only scaling that stays exact without widening the mantissa, and subtracted, so an accumulator can be brought back to exact zero and reused with each column's learned scale intact.

A symmetric n × n matrix can be accumulated into one stored triangle. That halves the memory and the input traffic, and is worth 2.0 to 2.2× on the CPU, but the deeper reason is exactness. In full storage the two halves would hold equal values in different limbs under different column scales, and symmetry would become a property the data had to keep earning. Stored once, it cannot drift. The input is taken to be symmetric and that is not checked, because checking it means reading the half the design exists to avoid reading.

Accumulation can be spread over worker threads, partitioned by column. Columns are independent in storage and each is touched by exactly one worker, so there is no locking on the hot path and the result is bit-identical to a serial run. On a 4096 × 64 accumulation over eight cores this takes 0.80 to 4.31 Gelem/s. Scaling stops where the working set leaves cache: at 65536 × 64 the input alone is 33.6 MB against 32 MiB of L3, and eight threads reach only 1.11 Gelem/s because every core is waiting on the same DRAM. A single core in cache is issue-bound, at 2.32 IPC with a 1% miss rate, which is why it scales at all until it doesn't.

Past that point the lever is traffic, not threads. Several input matrices can be folded into one pass over the accumulator: the limbs are read once, every addend applied in registers, and written once. The saving is in accumulator traffic only, so it shows up where accumulator traffic is the bound. On the DRAM-bound case, eight matrices per pass takes 1.14 to 3.01 Gelem/s and is still rising, while single-threaded and in cache the extra input streams cost more than they save.

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


[^ij]: `i` is the row and `j` the column, in the usual order. The two are not interchangeable here: a column owns an exponent, a width, its own allocations and exactly one worker thread, while a row is only an offset into each of those allocations.
