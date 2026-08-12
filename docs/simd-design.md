# SIMD and portability design record

**Status: design, not implemented.** The only pieces of this that have landed
are the bit-manipulation `decompose` (`8565474`) and the `kLimbBits` cleanup
(`60711df`). Everything below the "Implemented" section is a decision record,
not a description of the code.

Targets both AVX-512 on the CPU and a future CUDA implementation. The CUDA
constraint is the more restrictive one and should drive the abstraction.

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

## Limb width: separate storage from radix

CUDA wants 32-bit limbs. The forcing constraint is not throughput: **nvcc has
no `__uint128_t` in device code**, and every carry path in `limbs.cpp` goes
through one today. A 32-bit radix makes the double-width intermediate a plain
`uint64_t`, which device code handles natively.

For AVX-512, lane count alone is *not* an argument for narrow limbs. Pushing
`R` rows through a `W`-bit column costs `(W/64) * (R/8) = RW/512` instructions
with 64-bit limbs and `(W/32) * (R/16) = RW/512` with 32-bit. Doubling the
lanes also doubles the limb positions needed to span the same value width, and
the two cancel exactly. It also doubles each row's serial carry chain.

The real argument for a narrow radix is **headroom**. A 32-bit radix stored in
64-bit lanes leaves 32 spare bits per lane, so ~2^31 values can be accumulated
before overflow is possible. The inner loop becomes a bare `vpaddq` with no
carry-out test, no mask op, and no dependency between limb positions; carries
are propagated once, lazily, before readback.

So parameterize two constants, not one:

| storage | radix | deferred adds | density | target |
| --- | --- | --- | --- | --- |
| 64 | 64 | 1 (canonical) | 100% | today's CPU code |
| 32 | 32 | 1 (canonical) | 100% | CUDA canonical |
| 64 | 32 | 2^31 | 50% | AVX-512 carry-save |
| 64 | 52 | 2^11 | 81% | reduced radix, better density |

A redundant (non-canonical) representation stays exact — every add is still
exact, just stored redundantly — but `to_double`, `to_exact_decimal` and
`is_zero` must normalize first, since a nonzero encoding can denote zero.

### What is limb-width-dependent, and what is not

Limb-width-dependent, must be parameterized:

- `limb_t`, and `p[n-1] >> 63` in `is_negative`
- every `__uint128_t` carry/borrow/product intermediate
- `__builtin_clzll` / `__builtin_ctzll` where the operand is a limb
- the decimal chunk `10^19` / 19 digits (32-bit limbs want `10^9` / 9)
- the `5^27` stride in `to_exact_decimal` (32-bit limbs want `5^13`)

**Fixed at 64 bits regardless**, because these describe IEEE-754 doubles rather
than the accumulator: everything in `decompose` (52, 1075, `0x7FF`, the sign
shift), `DoubleParts::mantissa`, `extract_u64`'s return type, and the rounding
path in `to_double`. Keeping this boundary clean is most of the portability
work. `60711df` fixed the two places that had already blurred it.

### One structural change beyond a typedef

`add_shifted` takes the addend as two limbs (`lo`, `hi`). That holds only for a
64-bit radix: a 53-bit mantissa offset by up to 63 bits spans 116 bits, which
is 2 limbs. At a 32-bit radix, 53 bits offset by up to 31 spans 84 bits, which
is **3 limbs**. The addend should become a small fixed-size array sized
`ceil((53 + radix - 1) / radix)`.

### De-risking without a GPU

Template on the limb parameters and instantiate **both widths on the CPU**,
then assert the two produce bit-identical `to_exact_decimal` output over the
existing test corpus. Exact results must not depend on limb width, so this is a
strict equality test, and it makes the whole portability question testable
today. If the code is limb-width-clean, the CUDA port is a backend swap; if it
is not, that shows up now rather than after the kernels are written.

This argues for templates over a `CBFP_LIMB_BITS` macro: a macro forces one
width per build, so the two instantiations can never be compared in a single
test binary.

### Other CUDA notes

- Dynamic growth mid-kernel is impractical. `reserve_for` already exists as the
  pre-sizing escape hatch and is what a GPU path should require.
- Many small `cudaMalloc`s are expensive, which pushes against one allocation
  per limb position. The append-without-touching property still argues for it;
  a pool or arena is the likely compromise.
- The skew concern maps to partition camping rather than L1 set conflicts, but
  the medicine is the same: avoid power-of-two strides between limb arrays.
- With a 32-bit radix, thread `i` handling row `i` reads `limbs[k][i]` across a
  warp as 32 consecutive 32-bit words — 128 bytes, perfectly coalesced.

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
