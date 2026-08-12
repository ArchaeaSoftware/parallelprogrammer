# SIMD and portability design record

**Status: design, not implemented.** The only pieces of this that have landed
are the bit-manipulation `decompose` (`8565474`) and the `kLimbBits` cleanup
(`60711df`). Everything below the "Implemented" section is a decision record,
not a description of the code.

Targets both AVX-512 on the CPU and a future CUDA implementation. CUDA turned
out to impose no special constraint on the representation — measurement found
the two targets want the same thing — so the abstraction is driven by the radix
choice below rather than by either machine's word size.

Measurements are from a Ryzen 7 7700X (Zen 4: AVX512F/BW/DQ/VL, VBMI2,
VPOPCNTDQ; 256-bit datapath, so AVX-512 ops are double-pumped).

## The layout: limb-major `LimbColumn`

Today a column stores each entry's limbs contiguously — entry `i` occupies
`data[i*L .. i*L+L-1]`. The proposed layout inverts that: one array per limb
position, so structure `k` holds `rows` 64-bit words, the `k`-th limb of every
entry in the column.

```cpp
struct Column {
    int exponent;
    std::vector<LimbColumn> limbs;   // limbs[k] = k-th limb of every row
};
```

Each `LimbColumn` is a slice at **fixed significance spanning all rows**: row
`i`'s word in `limbs[k]` covers value bits `2^(exponent + 64k)` through
`2^(exponent + 64k + 63)`. They are ordered least-significant first, and the
topmost one is special — it carries the two's-complement sign bit.

Naming caveat: `LimbColumn` sits next to `Column` in `ColumnBlockMatrix` and a
`LimbColumn` is *not* a column of the matrix. Worth a comment at the
declaration; naming conventions may be revisited.

### Why this layout

The carry dependency runs along limbs within one entry, which is inherently
sequential. Rows are independent. Limb-major puts those on different axes: walk
limb positions sequentially for the carry chain, process 8 rows per `vpaddq`.
No cross-lane carry propagation, no shuffles.

Two structural wins matter more than the add kernel:

- **Widening becomes an append.** `ensure_limbs` currently re-strides every
  entry, O(N·L). With one allocation per limb position, existing arrays do not
  move at all: push back a new one and sign-fill it. The fill is exactly
  `_mm512_srai_epi64(top, 63)` — arithmetic shift right by 63 yields all-ones or
  all-zeros per lane. O(N·L) restride becomes O(N) append.
- **Rescale becomes a streaming funnel shift.**
  `new[k] = (old[k-w] << b) | (old[k-w-1] >> (64-b))` is `_mm512_shldi_epi64`
  (VBMI2, present on Zen 4). Sequential reads, no gather.

The awkward part is addend placement: each row's mantissa lands at its own bit
offset, spanning limb positions `off` and `off+1`, so each position needs masked
selects. Compute all offsets, start at the horizontal minimum, and exit via
`_kortestz` once every lane's carry has died.

Cost accepted: scalar `add(i, j, v)` touches L cache lines instead of one or
two. Bulk `add_matrix` gets much faster; single-element updates get slower.

## Allocation: 64-byte aligned, deliberately skewed

Each `LimbColumn` gets its own `aligned_alloc`, overallocated and offset:

```cpp
struct LimbColumn {
    void*          alloc;  // aligned_alloc(64, rows*8 + 512), kept for free()
    std::uint64_t* data;   // alloc + ((index & 7) * 64), still 64B-aligned
};
```

The offset must stay a multiple of 64 so alignment survives.

**Alignment** is worth ~17%: with cache-set conflicts controlled, an 8-byte
misalignment cost 17.3% and a 32-byte misalignment only 2.0%. The datapath is
256-bit, so a 64-byte access runs as two 32-byte halves; 32-byte alignment keeps
each half inside one line. Align to 64 anyway, it is free. Plain
`std::vector<uint64_t>` measured 48-byte alignment here, so it is not enough.

**Skew** is not about throughput, it is about variance. A limb array is
`rows * 8` bytes, and round row counts make those mutually 4KB-congruent —
`rows = 512` is exactly 4096. Running the identical benchmark twice in one
process, differing only in heap state:

```
first call:   no skew 28.8 GB/s    skewed 91.1 GB/s
second call:  no skew 79.6 GB/s    skewed 89.6 GB/s
```

Unskewed performance is a lottery drawn by the allocator, 2.8x between tickets.
Skewed is ~90 both times. 448 bytes on a 512KB array (0.09%) buys a
deterministic outcome.

`K=8` (offsets `(k & 7) * 64`) sufficed everywhere measured. Treat it as a
tunable, not a derived constant — see the caveats.

## Limb width: keep 64-bit storage, tune the radix

**Measured conclusion: limb width is not the lever; deferring carries is.**
Storage stays `uint64_t` on both targets, and the one parameter worth exposing
is the *radix* — how many of those 64 bits are used before headroom begins.

RTX 3060, 2^20 rows x 256 accumulations, ~512-bit accumulator held in
registers, reported as accumulator bits updated per second:

| scheme | limbs | storage / 512 bits | throughput |
| --- | --- | --- | --- |
| canonical, 64-bit limbs | 8 | 64 B (100%) | 29.1 Tbit/s |
| canonical, 32-bit limbs | 16 | 64 B (100%) | 27.0 Tbit/s |
| **carry-save, 52-bit radix** | 10 | 80 B (81%) | **82.2 Tbit/s** |
| carry-save, 32-bit radix | 16 | 128 B (50%) | 49.4 Tbit/s |

Two results, both against the direction this document previously took:

- **32-bit limbs do not help on a 32-bit machine.** They are 8% slower than
  64-bit canonically and 40% slower under carry-save. The 32-bit ALU implements
  a 64-bit limb at exactly proportional cost, so there is nothing to reclaim by
  matching the hardware width.
- **A reduced radix beats a narrow one.** 52-in-64 wins because it carries more
  usable bits per register and per byte moved (81% density against 50%), while
  2^11 deferred additions is ample for a batch. Normalising every few thousand
  accumulations is cheap amortised.

Caveat on the magnitudes: the benchmark adds a full-width value, whereas the
real kernel adds a 53-bit mantissa landing in two or three limbs with carries
that usually die immediately. That makes canonical look worse here than it will
be in practice, so treat 2.8x as an upper bound on the carry-save advantage.
The ordering should hold — carry-save over canonical, 52 over 32 — because a
53-bit mantissa spans two limbs at radix 52 and up to three at radix 32.

### Why width is neutral on both targets

The two architectures charge for a limb in proportion to its width, which is
why the measurement above comes out flat.

On the GPU the ALU is 32-bit, so one 64-bit add is already two chained 32-bit
adds — visible directly in the SASS:

```
IADD3   R10, P1, R2, R4, RZ           // low 32 bits, carry out to predicate P1
IADD3.X R0,  P1, R3, R5, RZ, P1, !PT  // high 32 bits, carry in from P1
```

A 64-bit limb is therefore honest rather than free: it costs exactly twice a
32-bit limb and covers exactly twice the bits. Nothing is reclaimed by matching
the hardware width.

On AVX-512 the same cancellation appears as lane arithmetic. Pushing `R` rows
through a `W`-bit column costs `(W/64) * (R/8) = RW/512` instructions with
64-bit limbs and `(W/32) * (R/16) = RW/512` with 32-bit. Doubling the lanes
also doubles the limb positions needed to span the same width, and the two
cancel — while doubling each row's serial carry chain.

### Why carry-save wins

Headroom, not width. A radix narrower than the storage leaves spare bits per
lane, so many values can be accumulated before overflow is possible: 2^11 at
radix 52, 2^31 at radix 32. The inner loop becomes a bare add — no carry-out
test, no mask op, and no dependency between limb positions — with carries
propagated once, lazily, before readback.

The density term is what then favours 52 over 32: fewer limbs for the same
accumulator width means fewer registers, fewer bytes moved, and a 53-bit
mantissa landing in two limbs instead of three.

A redundant (non-canonical) representation stays exact — every add is still
exact, just stored redundantly — but `to_double`, `to_exact_decimal` and
`is_zero` must normalize first, since a nonzero encoding can denote zero.

### Historical note

`__uint128_t` **is** supported in device code, verified on CUDA 12.9 / sm_86:
the exact carry shape `limbs.cpp` uses compiles, runs, and produces correct
carry-out at the `2^64-1 + 2^64-1 + 1` boundary, lowering to `add.cc.s64` /
`addc.cc.s64`. An earlier draft claimed the opposite and made it the
load-bearing argument for a 32-bit radix. Measurement then removed the
conclusion too: 32-bit limbs are slower than 64-bit on both targets.

### What is limb-width-dependent, and what is not

Limb-width-dependent, must be parameterized:

- `limb_t`, and `p[n-1] >> 63` in `is_negative`
- every `__uint128_t` carry/borrow/product intermediate (portable to device
  code as-is, but it should track the radix rather than be hard-coded)
- `__builtin_clzll` / `__builtin_ctzll` where the operand is a limb
- the decimal chunk `10^19` / 19 digits, and the `5^27` stride in
  `to_exact_decimal` — both sized to the largest power fitting in a limb

The measured conclusion is that storage stays 64-bit everywhere, so most of
this list is dormant. It still matters: a reduced radix means limbs are no
longer full 64-bit values, so anything reading a limb as a whole word — the
sign test, `bit_length`, the decimal conversion — must work in radix bits
rather than storage bits.

**Fixed at 64 bits regardless**, because these describe IEEE-754 doubles rather
than the accumulator: everything in `decompose` (52, 1075, `0x7FF`, the sign
shift), `DoubleParts::mantissa`, `extract_u64`'s return type, and the rounding
path in `to_double`. Keeping this boundary clean is most of the portability
work. `60711df` fixed the two places that had already blurred it.

### The addend span, and another reason 52 is a good radix

`add_shifted` takes the addend as two limbs (`lo`, `hi`), sized
`ceil((53 + radix - 1) / radix)` in general. A 53-bit mantissa at intra-limb
offset `o` occupies bits `o .. o+52`:

- radix 64: `o <= 63`, top bit 115, spans **2 limbs**
- radix 52: `o <= 51`, top bit 103, spans **2 limbs**
- radix 32: `o <= 31`, top bit 83, spans **3 limbs**

So a 52-bit radix keeps the existing two-limb addend and needs no structural
change here, while a 32-bit radix would force a third. 52 is well matched to a
53-bit significand.

### De-risking without a GPU

Template on the radix and instantiate **more than one on the CPU**, then assert
they produce bit-identical `to_exact_decimal` output over the existing test
corpus. Exact results must not depend on the radix, so this is a strict
equality test, and it makes the representation question testable today —
including carry-save against canonical, which is where a normalization bug
would otherwise hide until it corrupted a result silently.

This argues for templates over a `CBFP_LIMB_BITS` macro: a macro forces one
configuration per build, so two instantiations can never be compared in a
single test binary.

### Other CUDA notes

- Dynamic growth mid-kernel is impractical. `reserve_for` already exists as the
  pre-sizing escape hatch and is what a GPU path should require.
- Many small `cudaMalloc`s are expensive, which pushes against one allocation
  per limb position. The append-without-touching property still argues for it;
  a pool or arena is the likely compromise.
- The skew concern maps to partition camping rather than L1 set conflicts, but
  the medicine is the same: avoid power-of-two strides between limb arrays.
- Coalescing is fine either way: thread `i` handling row `i` reads `limbs[k][i]`
  across a warp as 32 consecutive words — 256 bytes at 64-bit storage, 128 at
  32-bit. Both are fully coalesced with no wasted bytes, which is part of why
  the width comparison came out flat.

## Implemented

**Bit-manipulation `decompose`** (`8565474`). Reads the IEEE-754 fields
directly; normal and subnormal collapse to `e = max(biased, 1) - 1075`. Took
~28% off per-element cost end to end (12.3 → 8.8 ns narrow, 14.4 → 10.6 wide).

In isolation it is ~9x faster than the `frexp`/`ldexp` version, but that
overstates the win: `frexp`/`ldexp` were already inlined as builtins inside
`accumulate`, so the microbenchmark measured an out-of-line call the real code
never made.

## Remaining work, in order

1. **Column-at-a-time bulk path.** `accumulate` is called per element and
   re-does the bounds check, `isfinite`, column lookup and `ensure_limbs` test
   every time. The layout change buys nothing until this exists. A vectorized
   decompose also gives a branch-free finite check: one `kortestz` per 16
   elements instead of an `isfinite` each.
2. **Limb-major layout + SIMD kernels**, per above.
3. **Threading across columns** — columns are already independent, so this is
   less work than intrinsics and composes with them.

A 16-wide AVX-512 `decompose` prototype exists and verifies against the scalar
version on ~2M inputs (specials, random bit patterns including subnormals and
non-finite, random finite values across exponents -1080..1020). It runs at
0.275 ns/elem versus 0.689 for the scalar bit-twiddle, so ~2.5x — worth having
once step 1 makes 16-wide output consumable, not before.

## Measurement caveats

Two benchmarks in this investigation gave confidently wrong answers before
being caught. Both are worth remembering when re-measuring.

- An alignment benchmark reported *misaligned as 4x faster*. It had allocated
  each array separately, landing them 4KB-congruent, so it was measuring
  cache-set conflicts with alignment held constant.
- A skew benchmark produced non-monotonic results across limb counts that
  contradicted the L1-set model behind them (predicted: more limbs, more need
  for skew; measured: no benefit at 16, ~5% cost at 32). Reproducible to ~3%,
  so not noise — the model was simply wrong. It was discarded.
- The claim that nvcc rejects `__uint128_t` in device code was asserted from
  recollection, not measured, and is false on CUDA 12.9. It had been made the
  load-bearing argument for a 32-bit radix. Benchmarking then took the
  conclusion as well: 32-bit limbs are slower than 64-bit on both targets. Two
  rounds of plausible reasoning, both overturned the moment either was
  measured — check the toolkit against the toolkit, and the hardware against
  the hardware.

The synthetic carry walk used throughout has none of the real work: no
decompose, no per-lane variable shifts, no masked selects. Re-measure the skew
constant against the actual kernel before treating any of these numbers as
settled.

## Current baseline

4096x64, 32 batches, after the `decompose` change:

| case | ns/elem | column width |
| --- | --- | --- |
| narrow (1 binade) | 8.8 | 192 bits |
| wide (120 binades) | 10.6 | 512 bits |
| wide, pre-reserved | 10.2 | 512 bits |
| rescale every batch | 30.5 | 2048 bits |
| widen every batch | 10.9 | 1024 bits |

Marginal cost of an extra limb is ~0.42 ns (~2 cycles), so at 3 limbs the limb
arithmetic is only ~10% of runtime. That is why the per-element overhead in
step 1 comes before the SIMD work in step 2.
