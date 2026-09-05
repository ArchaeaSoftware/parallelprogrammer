# SIMD and portability design record

**Status: largely implemented.** The limb-major layout, the skewed aligned
allocation and the AVX-512 survey and accumulate kernels have all landed, and
the measurements below are from the real kernels unless a section says
otherwise. The radix question is now closed too: carry-save at
52 was implemented and measured, and 64 stays. Threading, producer-supplied
surveys on both targets, symmetric storage and multi-batch folding have all
landed since. This
header previously read "design, not implemented", which stopped being true at
`fdc8b75`.

Targets AVX-512 on the CPU and CUDA on the device, both implemented. CUDA
turned out to impose no special constraint on the representation — measurement
found the two targets want the same thing — so the abstraction is driven by the
radix choice below rather than by either machine's word size. This paragraph
read "a future CUDA implementation" until that stopped being true.

Measurements are from a Ryzen 7 7700X (Zen 4: AVX512F/BW/DQ/VL, VBMI2,
VPOPCNTDQ; 256-bit datapath, so AVX-512 ops are double-pumped) and an
RTX 3060 (sm_86, CUDA 12.9).

## Decisions of record

| decision | status |
| --- | --- |
| **64-bit limb storage (`uint64_t`) on every target** | **settled** |
| Limb-major layout: one array per limb position | settled |
| 64-byte aligned allocation, skewed per limb array | settled |
| `LimbColumn` as the type name | settled, naming may be revisited |
| Bit-manipulation `decompose` | landed (`8565474`) |
| Flat free-function kernels; runtime dispatch via target attributes | settled |
| No templates in the public API; radix templated in kernels only | settled |
| Radix: 52-bit carry-save vs 64-bit canonical | **closed: stays at 64** |
| Whether limb arrays get separate allocations on CUDA | measured neutral, kept separate |
| Survey supplied by the producer; accumulation matrix pre-sized from it | landed, both targets |
| Symmetric matrices store one triangle, LAPACK `uplo` | landed, both targets |
| Readback returns the residual alongside the rounded double | landed |
| Multi-batch folding, K matrices per pass | landed, CPU |

Storage width is closed: 32-bit limbs were measured slower than 64-bit on both
targets, so `uint64_t` is the plan of record everywhere and the limb type does
not need to be a template parameter.

The **radix** was a separate question, open for most of this document's life
and now closed — on the opposite target from the one it spent its time on, and
against the change. The reasoning below is kept because the path from "2.8x on
the GPU" to "measured, and not worth doing on either" is the useful part.

An earlier draft had it that carry-save was a GPU win (2.8x, below) and that
the CPU was simply unmeasured. Bounding it on both targets says the reverse.
The bound is the cheap version of the experiment: take the real kernel and
delete the carry chain outright — two limb positions written unconditionally,
no carry-out test, no break, no dependency between them. That is
arithmetically wrong and is not an implementation; it is the *shape*
carry-save would have, so whatever it saves is the most carry-save could ever
save.

| target | what deleting the carry chain buys |
| --- | --- |
| CUDA, 2-limb columns, 4096–65536 rows | −2.4% to −0.3% |
| CUDA, 8-limb columns, spread 0–400 | −5.2% to +2.8% |
| AVX-512, 1–4 limbs | **25–36%** |
| AVX-512, 8 limbs, spread 200–400 | **45–52%** |

**On the GPU the carry chain is free.** Removing it entirely is worth nothing
outside noise, and the sign is not even consistent. That kernel is bound by
memory, so the arithmetic hides behind it.

**On the CPU it is a third to a half of the kernel.** AVX-512 runs at ~2.3 IPC
against a double-pumped 256-bit datapath, so the six extra 512-bit ops a
carrying limb position costs — the complement, the carry seed, two `cmplt`
compares, the second add, the `pending` test — are all paid at full price.

Two things the bound does not charge, so the real figure is lower than 52%:
at radix 52 a block's offsets are `shift/52` rather than `shift/64`, so the
*range* a block spans grows ~23%, and that range is exactly what drives the
AVX-512 loop's trip count — it bites hardest in the divergent case where the
bound looks best. Normalization before readback is not charged either.

So the CPU path stays canonical at radix 64 for now, which is what the code
does, but the reason has changed again -- see below, where it was implemented
and measured rather than bounded.

### Radix 52: implemented, measured, not pursued

The bound above was collected by deleting the carry chain, which is not an
implementation. A real radix-52 carry-save AVX-512 accumulate was then written
and checked against an independent integer reference before being timed, and
it charges three things the bound never did:

- **The offset is a division.** `shift/52` where radix 64 has `shift >> 6`.
  There is no SIMD integer divide, so it is a magic multiply plus a shift, and
  a second multiply to recover `bit`. Radix 64 pays none of this.
- **The same range needs more limbs**: 1 -> 2, 2 -> 3, 4 -> 5, 8 -> 10. That
  is 1.25-2x the accumulation matrix memory, which costs nothing while it fits in
  cache and a great deal when it does not.
- **Normalization**, which turned out to be negligible: under 0.002 ns/elem
  amortized over the 2048 addends the 11 bits of headroom allow.

ns/elem, single-threaded, radix 52 including amortized normalization:

| case | r64 | r52 | | past L3, r64 | r52 | |
| --- | --- | --- | --- | --- | --- | --- |
| 1 limb, spread 0 | 0.729 | 0.758 | **0.96x** | 0.769 | 0.945 | **0.81x** |
| 2 limbs, spread 40 | 1.102 | 0.826 | 1.33x | 1.195 | 0.894 | 1.34x |
| 4 limbs, spread 150 | 1.871 | 1.208 | 1.55x | 2.434 | 1.680 | 1.45x |
| 8 limbs, spread 400 | 3.229 | 1.914 | 1.69x | 4.650 | 3.928 | 1.18x |

**Carry-save is not a global win, and the bound's shape was misleading.** It
loses outright on one-limb columns -- the narrowest and most common case, and
the one the CPU's best single-thread rate is measured on -- because there is
no carry chain to remove there and the width doubles anyway. And the widest
case, where the bound looked best at 45-52%, is where the extra limbs cost
most once the accumulation matrix outgrows L3: 1.69x collapses to 1.18x.

Where it pays is the middle, 2-4 limbs, at 1.33-1.55x, and that holds up past
L3. So the useful form of this is not a global radix but a **per-column**
choice on width, which the container is already shaped for since every column
carries its own exponent and limb count. That doubles the kernel surface and
every readback path, which is the cost to weigh against 1.3-1.5x on some
columns and a loss on others.

Two asymmetries in the comparison, in both directions: the radix-52 prototype
does not compute the detector flags the shipped radix-64 kernel does (worth
0.8-4.8%), and it does not interleave two row blocks the way that kernel does
(worth ~7-10% at 2-4 limbs). Treat 1.3-1.5x as the shape of the answer, not a
figure to three digits.

**Decision: the radix stays at 64 and this line of work is closed.** Per-column
selection is the only form that measures well, and it buys 1.3-1.5x on middling
columns at the price of a second accumulate kernel, a second normalization
path, a second readback, and a width-dependent branch in the container -- for
something that is a loss on the narrowest columns and nearly gone on the widest
once past L3. The kernels stay templated on the radix, so the hook remains if
the column-width distribution of a real workload ever argues differently, but
nothing further is planned. Reopening this should mean new evidence about
column widths, not a fresh look at the same numbers.

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

Naming caveat: `LimbColumn` sits next to `Column` in `AccumulationMatrix` and a
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

Cost accepted: a scalar element-at-a-time add touches L cache lines instead of
one or two, and charges the whole column for one addend. Bulk `add_matrix`
gets much faster; single-element updates get slower. That entry point was
later removed for exactly those reasons: a caller who wants one cell builds
the column, which is all it ever did.

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

**Skew** keeps limb arrays off mutual 4KB congruence. A limb array is
`rows * 8` bytes, so round row counts make them congruent — `rows = 512` is
exactly 4096 — and one 8-row block touches every live limb array at the same
row offset, which then lands every one of them in the same L1 set against an
8-way cache.

Measured against the real accumulate kernel, with the arrays carved from one
arena at a controlled stride so the allocator is not a variable, 4096 rows and
exponents spread across the full column width so every limb array is live:

| limbs | congruent | skewed | skew buys |
| --- | --- | --- | --- |
| 2 | 1.158 | 1.157 | nothing |
| 4 | 2.06 | 1.92 | 6.8% |
| 8 | 3.61 | 3.42 | 5.3% |
| 16 | 7.16 | 6.54 | 8.7% |

Two qualifications, both of which matter more than the headline:

- **The effect is 4KB congruence specifically, not the offset.** A stride
  delta of 4096 measures identically to a delta of 0 — still congruent — while
  64, 128, 256 and 512 are indistinguishable from each other. `kSkewStep = 64`
  is not a tuned constant, it is the cheapest way to not be a multiple of 4096.
- **It only appears when the column is wide *and* its exponents diverge.**
  With the exponents clustered, a block touches ~2 limb arrays whatever the
  column's width, and the skew measures as noise. That is the same regime that
  makes the AVX-512 accumulate lose to scalar, so the skew pays off mainly
  where the vector kernel should not be running anyway.

An earlier draft of this section claimed 2.8x, from a bandwidth-bound
synthetic carry walk (28.8 vs 91.1 GB/s). That does not reproduce against the
real kernel, which is issue-bound at ~2.3 IPC rather than waiting on L1. The
skew is worth keeping at 448 bytes on a 512KB array (0.09%); it is not worth
2.8x.

## Limb width: 64-bit storage everywhere (settled)

**Measured conclusion: limb width is not the lever; deferring carries is.**
Storage is `uint64_t` on both targets — the plan of record — and the only
parameter worth exposing is the *radix*, meaning how many of those 64 bits are
used before headroom begins.

RTX 3060, 2^20 rows x 256 accumulations, ~512-bit accumulation matrix held in
registers, reported as accumulation matrix bits updated per second:

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
  2^11 deferred additions is ample for a batch. Normalizing every few thousand
  accumulations is cheap amortized.

**This benchmark held the accumulation matrix in registers, and that is what makes it
misleading.** With no memory traffic the carry chain is the only cost, so
removing it is the whole game. The real kernel keeps the accumulation matrix in global
memory, where traffic dominates and the chain is free — measured, and recorded
under the radix decision above. Treat the table below as a statement about
register-resident accumulation, not about this workload.

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
accumulation matrix width means fewer registers, fewer bytes moved, and a 53-bit
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
than the accumulation matrix: everything in `decompose` (52, 1075, `0x7FF`, the sign
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

## Differentiating scalar, AVX-512 and CUDA

Keep the code flat. The kernel surface is a handful of free functions over raw
pointers — accumulate a column, shift a column left, sign-fill a new limb
array, normalize (if carry-save), bulk readback — taking
`std::uint64_t* const* limbs, std::size_t nlimbs, std::size_t rows, ...`, which
is what limb-major with separate allocations already hands you. No virtuals, no
CRTP, no policy types.

**CPU: one kernel per variant, selected once at runtime.**

```cpp
// kernels_scalar.cpp   compiled at baseline
// kernels_avx512.cpp   compiled with -mavx512f -mavx512dq -mavx512vbmi2
static const auto add_limbs =
    __builtin_cpu_supports("avx512f") ? add_limbs_avx512 : add_limbs_scalar;
```

Two ways to keep AVX-512 codegen out of a baseline binary, both verified on
this toolchain:

- **Per-TU flags.** The AVX-512 kernels get their own translation unit with
  its own flags and no attributes. Better codegen: inside that TU the compiler
  may auto-vectorize everything, not only the hand-written intrinsics, and
  inline freely. **Hazard:** `-mavx512f` licenses AVX-512 anywhere in that TU,
  including static initializers or inlined header templates. Keep it to pure
  leaf kernels with no static init, reached only after the CPU check.
- **Per-function `__attribute__((target(...)))`.** Everything stays in one TU
  compiled at baseline; `zmm` codegen stays confined to the attributed function
  (measured: 0 occurrences in the scalar kernel, 8 in the AVX-512 one). Safer
  by construction and needs no per-file build rules, at the cost of a more
  constrained optimizer.

Per-TU flags are the default choice; the attribute is the fallback where a
kernel must sit next to code that runs unconditionally.

Dispatch is one indirect call per *column*, not per element.

### CRTP was considered and rejected

CRTP would replace the function pointer with a compile-time-resolved call.
Measured at column granularity (4096 rows, 200k dispatches):

```
function pointer (runtime dispatch)    366.6 ns/dispatch
CRTP (compile-time, no indirection)    368.4 ns/dispatch
direct call (upper bound)              367.0 ns/dispatch
```

CRTP is 0.5% *slower* — that is, indistinguishable. There is no indirection
cost to remove when one dispatch covers 4096 elements; the hot loop lives
entirely inside a single kernel call. CRTP earns its keep when the call sits in
the hot loop, which is not the case here.

It also does not eliminate the need for per-TU flags or target attributes:
those come from how the compiler is invoked, not from the pattern. Plain free
functions in a separately-compiled TU get exactly the same codegen, verified.

The cost is real, though. CRTP makes the backend a *type*, and since AVX-512
availability is a runtime property the branch still has to exist — now
selecting between two types, which forces `AccumulationMatrix` to become a
template or to be type-erased behind a vtable, the very thing CRTP was meant to
avoid. The whole class body would also be instantiated per backend when only
the handful of kernels differ.

**CUDA is a different seam, at a different level.** Swapping a CPU kernel
leaves the memory and the object identical, so it is a true backend swap.
CUDA's data lives in device memory, so putting it behind the same function
pointer would disguise host/device transfers as ordinary calls. The CUDA path
wants its own container with device-resident limb arrays; what it shares with
the CPU is the *algorithm* — when to rescale, when to widen, how the exponent
moves — not the memory. Unifying all three behind one dispatch table would
produce exactly the leaky abstraction this design is trying to avoid.

## Templates: only in the kernels, never in the API

Settling storage at 64 bits removed the case for templating the public class —
there is no type left to vary. Templating `AccumulationMatrix` would force
header-only or explicit instantiation and make the most-read code less flat for
no benefit.

The one legitimate parameter is the radix, and it belongs on the kernel bodies
only. Those are already isolated behind the flat function interface, so they
can be templated and explicitly instantiated at radix 52 and 64 without a
template ever appearing in a public header. If the radix question closes at 64,
the parameter is deleted and nothing else changes.

### De-risking without a GPU

Template on the radix and instantiate **more than one on the CPU**, then assert
they produce bit-identical `to_exact_decimal` output over the existing test
corpus. Exact results must not depend on the radix, so this is a strict
equality test, and it makes the representation question testable today —
including carry-save against canonical, which is where a normalization bug
would otherwise hide until it corrupted a result silently.

This argues for templates over a `TRUESUM_LIMB_BITS` macro: a macro forces one
configuration per build, so two instantiations can never be compared in a
single test binary.

### Other CUDA notes

- Dynamic growth mid-kernel is impractical, and two ways around it are closed
  on this toolkit rather than merely inadvisable — see "The CUDA path" below,
  which also carries the plan for growing a column *between* launches.
- Many small `cudaMalloc`s are expensive, which pushed against one allocation
  per limb position. **Resolved: the stream-ordered pool is the compromise.**
  `cudaMallocAsync` makes per-limb-position allocation affordable, so the
  device keeps the append-without-touching property rather than packing a
  column into one strided arena. A first cut did use a single arena per column
  with limb `k` at `base + k*pitch`, which is smaller and simpler to address
  but makes widening copy the whole column — surrendering the one structural
  win limb-major exists for.

  The extra indirection this costs the inner loop — a `limb_t* const*` load per
  limb position instead of a multiply-add — measured as nothing: 4.394 against
  4.402 Gelem/s at 4096x64, and every other shape within 2% with no systematic
  direction. A warp's threads all work different rows of the *same* limb
  column, so the load is a broadcast of one value, and the whole pointer array
  is `nlimbs * 8` bytes — 128 for a 16-limb column — hence L1-resident for the
  life of the kernel.
- **The skew does not carry over, and the device path deliberately omits it.**
  An earlier draft assumed the concern maps to partition camping and that the
  medicine was the same. It does not: the CPU needs the skew because eight
  lanes touch `nlimbs` arrays at the same row offset in succession and collide
  in an 8-way L1 set. A warp instead reads 32 consecutive rows of *one* limb
  column as a single coalesced transaction, and visits limb positions
  sequentially within a thread, so there is no equivalent collision — and
  Ampere hashes physical addresses across channels precisely so software need
  not. What the device does need is 256-byte alignment for coalescing, which
  the allocator gives for free and which a mis-sized skew would break. This was
  a decision not to add unmeasured complexity, not a measurement.
- Coalescing is fine either way: thread `i` handling row `i` reads `limbs[k][i]`
  across a warp as 32 consecutive words — 256 bytes at 64-bit storage, 128 at
  32-bit. Both are fully coalesced with no wasted bytes, which is part of why
  the width comparison came out flat.

## The CUDA path

Implemented and cross-checked against the CPU container by comparing stored
limbs, which is stronger than comparing readback: identical limbs at an
identical column exponent imply identical everything downstream, and unlike a
decimal comparison it cannot be passed by two implementations wrong in the same
way. Numbers below are an RTX 3060.

### Pre-sizing is forced, not chosen

A single accumulation pass can widen a column by many limbs at once, because
`max_addend_bits` is a *range* from the column exponent to the highest bit
reached rather than an increment. Measured on the CPU container: a fresh column
given one pass spanning 2^-1074 to 2^1023 goes from 1 limb to 33; a column
settled at 1.0 that then sees 2^1000 goes to 16; a single `add(2^900)` goes to
15. So there is no cheap slack to pre-allocate — sizing for the worst case is
2112 bits an entry against the 128-192 real columns use, which is the flat
allocation this structure exists to avoid.

Two escapes were measured and both are closed:

- **Device-heap memory is unreachable from the host.** A kernel can read a
  pointer from device-side `malloc`, but `cudaMemcpy` off it returns
  `invalid argument`, so device-allocated limb arrays would break readback.
- **CDP parent-child synchronization is gone.** Device-side
  `cudaDeviceSynchronize` does not merely warn on CUDA 12.9, it does not
  compile. A kernel cannot launch a child to grow a column and wait for it.

### The host-input path: read in place, do not copy

Input always arrives from host memory, and is read exactly once — the survey
runs on the host, so the device never needs a second look at it. That single
fact is what licenses the whole design: the kernel reads the host buffer
directly over PCIe as it works, rather than a transfer copying it to device
memory first and the kernel reading it from there.

Measured at 65536x64 from a pinned host buffer, copying 33.6 MB and then
reading it from device memory takes 1768 us; reading it in place takes 1257.
The transfer alone costs 1258, so the accumulation matrix arithmetic hides entirely
behind the bus instead of queueing after it.

| | 4096x64 | 16384x64 | 65536x64 |
| --- | --- | --- | --- |
| staged: the library copies the caller's buffer | 2.85 | 2.71 | 1.45 |
| borrowed from `acquire_input`, read in place | **3.13** | **3.27** | **3.27** |

The first row is now history: **`add_matrix_col_major` requires page-locked,
device-mapped input and rejects anything else.** No staging copy remains to
fall back on, so the second row is what every caller gets. Either the
accumulation matrix owns the buffer (`acquire_input`) or the caller does
(cudaHostAlloc/cudaHostRegister with the mapped flag); either way the call
returns an `InputRead`, a CUDA event recorded after the launch, and the rule is
the same for everyone: do not write the buffer until it clears.

That uniformity matters as much as the speed. The staged form had the library
making a copy the caller could have avoided by knowing to ask, and the version
before it inferred safety from the pointer's memory type -- one function, two
lifetime contracts, nothing in the signature to tell them apart, 63% of entries
corrupted. One rule, stated in the return type, replaces both. At 65536x64:

| | us a batch |
| --- | --- |
| staged, caller's own buffer | 3203 |
| in place, waiting on the handle each batch | 2024 |
| in place, two buffers in rotation | **1585** |

2.02x, and 1585 us is within reach of the ~1760 us the kernel spends on the
bus, so the host has stopped being the limiter. The event costs 0.16 us to
create and destroy, which is why one per submission is affordable rather than
needing a pool.

3.27 Gelem/s is 26.1 GB/s of doubles against the 26.7 this link achieves, so
the path sits at 97% of PCIe — and flat, where the copying version falls off
as batches outgrow whatever was hiding the copy.

Neither figure counts the caller generating its data, which is not the
library's cost. Getting that wrong produced two contradictory revisions of this
section: first a benchmark that filled one buffer and resubmitted it, which
measured the right thing by accident while demonstrating the unsafe reuse
pattern; then a "correction" that timed a loop copying pre-generated test data
into the buffer each batch, concluded the path was producer-bound rather than
bus-bound, and was measuring the harness. What is timed now is the library:
for a borrowed buffer, the survey and the kernel; for a staged one, those plus
the copy the library itself performs.

3.27 Gelem/s is 26.1 GB/s of doubles against the 26.7 this link achieves, so
**PCIe is the floor for this path**: eight bytes an element must cross it, and
no kernel change beats ~3.3 Gelem/s while input starts on the host. If that
ever stops being true it will be through GPUDirect from a network adapter, not
a faster staging path.

Ordinary host memory still works and is staged through a mapped buffer, which
is why the staged row falls off as batches grow: that host copy becomes the
limit. `acquire_input` lends out a buffer the device can read directly, which
`add_matrix_col_major` recognizes and reads in place.

#### What the survey is actually hidden behind, measured through the container

The 1257 us above is a bare prototype kernel. Timed through the container at
the same shape, so the per-column descriptor handling and the host survey are
in it too:

| | us a batch |
| --- | --- |
| host survey standalone, two-pass | 1376 |
| the survey as the library runs it, plus the launch | 705 |
| staging copy of 33.6 MB | 1331 |
| the kernel streaming the mapped buffer over PCIe | ~1760 |
| ordinary-buffer path, end to end | 3179 |

Three things follow, and the first corrects the framing above.

**There is no transfer for the survey to hide behind.** Read-in-place removed
it; the only `cudaMemcpyAsync` left on this path carries the 1 KB of
`DeviceSurvey`. What hides the survey is the *previous batch's kernel*, via the
two slots -- ~705 us of survey against ~1760 us of kernel, so it is free with
room to spare. The conclusion in this section survives; the mechanism it gives
does not.

**Only on the borrowed path.** From an ordinary buffer the staging copy alone
(1331 us) plus the survey already exceed the kernel, and end to end the batch
costs 3179 us -- host-bound, with the device waiting. The staged/borrowed table
above shows the same thing as throughput; this is where it comes from.

**The library's survey is cheaper than a standalone one** -- 705 against 1376 --
because `validate_host_survey` passes the column's current exponent as the
floor, so the significand pass is skipped once a column is settled.

**Fusing the copy with the survey was tried properly and does not pay.** The
staged path walks 33.6 MB twice, to copy and then to survey. A real fusion
shares the *loads*: eight values arrive in a vector register on their way to the
destination and the survey's first pass reads them there, so that pass costs no
memory traffic at all. It was written that way -- streaming stores, alignment
prologue for triangular slices, `sfence` -- and measured at each step, against a
mapped destination:

| | mapped dest |
| --- | --- |
| interleaved per column, no shared loads | ~0.92x |
| shared loads, ordinary stores | 1.01x |
| + non-temporal stores | 1.03x |
| + two independent accumulation matrix sets | **1.09x** |

Two things surfaced getting there. The copy is not the weak part: an AVX-512
loop with streaming stores beats `memcpy` outright, 1136 us against 1172. And
the single-accumulation matrix fused version was slower than running both operations
*separately* -- 2758 against 1136 + 1340 -- because the survey's min and max are
loop-carried, and one set of them serializes iterations the copy could
otherwise overlap. Two sets recover most of that, the same reason the
accumulate kernel interleaves two row blocks.

So the fusion does work, at 1.09x of host copy+survey. It is **invisible end to
end**: 2952-3167 us a batch against a 2941-2997 baseline, because the batch is
not purely host-bound and the kernel's PCIe time absorbs the saving. Reverted
on that basis -- a third kernel variant, non-temporal stores with an alignment
prologue and a memory-ordering fence, for a component nobody measures.

The staging copy was therefore never going to be made cheaper, and has since
been removed instead: the host path requires pinned input and reads it where it
lies.

Resolved, and it was nothing. 1257 us for the bare prototype against ~1760
"through the container" looked like real overhead and was an arithmetic
artifact: the 1760 was never measured, it was 2467 minus 705, subtracting a
host-side timing taken while a previous kernel was in flight from a serialized
total on a different path. Measured directly at the same shape:

| | us |
| --- | --- |
| standalone prototype | 1256.4 |
| library kernel, pre-sized | 1261.4 |
| library adaptive, total 1917.7 less host survey 658.4 | 1259.2 |

Within 0.4%. The per-column descriptor handling named as the candidate costs
nothing measurable. The lesson is the same one the per-block cost taught: a
figure obtained by subtracting two measurements is only as good as their having
been taken under the same conditions, and it is not a measurement.

The pinning measurement that justified the earlier design still stands and
still matters, because the staged path uses it: `cudaMemcpyAsync` from pageable
memory returns only after 3.92 ms of a 3.96 ms transfer, leaving no window to
survey in, while pinned returns in 0.01 ms of 2.51 and is 58% faster besides.

Validation stays fail-fast and costs nothing. The host survey finishes before
the accumulate is launched, so a batch that does not fit is rejected without
having touched the accumulation matrix — which a device-side survey can only promise
by draining the pipeline to ask.

**Why threading that survey buys nothing, stated carefully because the obvious
summary is wrong.** The survey reads every element of a batch, and looks like
exposed serial host work now that reading in place has removed the transfer it
used to hide behind. Threading it changes end-to-end throughput not at all —
but not because it resists threading:

| survey, ms per batch | 1 thread | 2 | 4 | 8 | kernel streams it in |
| --- | --- | --- | --- | --- | --- |
| 4096x64 | 0.055 | 0.031 | 0.029 | 0.055 | 0.079 |
| 16384x64 | 0.160 | 0.078 | 0.055 | 0.069 | 0.314 |
| 65536x64 | 0.729 | 0.356 | 0.188 | 0.150 | 1.257 |

It threads well, 4.9x at the largest shape. It is that the single-threaded
survey already fits inside the kernel at every shape, and two slots in rotation
put it there: batch N+1 is surveyed while batch N streams. Making hidden work
faster moves nothing, and end to end it measured flat on this path, +24% on the
staged path at 16384x64 and -8% at 65536x64.

So the conclusion is conditional. The survey is free *because PCIe is slow*.
The margin is thinnest at 4096x64 — 0.055 against 0.079 — and a slower host, or
a column wide enough to need the survey's second pass, could cross it. And if
input ever reaches the device faster than the bus, the survey becomes the
critical path and the 4.9x is there to collect. It is absent from the code
because today it would buy nothing, not because it does not work.

**Ping-pong buffers, and why the alternative was a trap.** Reading in place
means the kernel is still streaming a buffer after the call that submitted it
returned. The first version chose the in-place path by probing the pointer's
memory type, which gave one function two lifetime contracts: ordinary memory
was copied and could be reused at once, mapped memory could not, and nothing
said so. Switching allocator for speed silently corrupted 63% of entries.

`acquire_input` lends out a buffer the accumulation matrix owns and blocks only until
the device has finished with the one being handed back, so the host fills one
while the device streams the other. Ownership rather than memory type decides
whether a buffer is read in place, so a pointer the caller passes is theirs
again on return whatever kind of memory it is.

### The grid rule: cap the whole grid, not its x extent

Both kernels grid-stride, so any smaller grid is correct; it only gives each
thread more rows. Capping the *total* at 1024 blocks is what makes the survey's
reduction pay for itself — uncapped, 65536 rows over 64 columns launches 16384
blocks, each folding one value per thread and then paying a full eight-step
shared-memory tree to reduce it.

| | survey kernel | device path, end to end |
| --- | --- | --- |
| 16384 rows | | 5.543 -> 6.172 Gelem/s |
| 65536 rows | 195 -> 113 us | 5.853 -> 6.715 |
| 262144 rows | 746 -> 423 us | 5.931 -> 6.840 |

The accumulate is not a reduction and was expected to be indifferent; measured
through the host path, which uses the AVX-512 survey and so isolates it, it is
unchanged or slightly better.

### Getting results back to the host

Three rules, all measured, for the small per-column results both kernels
produce:

- **Never atomics over PCIe.** Reducing into mapped host memory costs 673 us
  against 2.59 in device memory — 260x. (They are *correct*; the cost is the
  whole objection.)
- **Mapped host memory wins for small results and falls off a cliff.** Writing
  1, 16 or 64 values costs ~4.8 us against ~11.9 for device memory plus a
  copy back; at 256 it is 30.2 against 13.2, and at 1024, 121.9 against 13.0.
  It is transaction count that matters, not volume — 1024 slots is 16 KB, under
  a microsecond of bandwidth.
- **So stage in device memory, elect one block to deliver.** Every block
  reduces into device memory with ordinary atomics, then `__threadfence()` and
  one `atomicAdd` on a ticket; the block drawing `gridDim - 1` copies the
  finished per-column results to mapped memory and re-arms the sentinels for
  the next launch, so there is no H2D initialization either. The election costs
  ~3.4 us and the PCIe write itself 0.64. Do **not** mark the device
  accumulation matrix `volatile` — copied from the NVIDIA sample, it cost 5.71 us by
  defeating caching on every atomic, and turned a 1.12x win into a 0.90x loss.

The election serializes on one counter, so it is O(blocks); it wins below about
2000 blocks and loses above, which the grid cap keeps it under.

### The survey belongs to whoever produced the matrix (landed)

Every path here pays for a survey — a full read of the batch to learn each
column's lowest true-ulp exponent and highest bit — and every path pays for it
*again* on top of the read the accumulate already does. On the host path that
is 0.729 ms at 65536x64, currently hidden behind the bus. On the device path it
is a kernel plus a round-trip whose verdict the host must see before the
accumulate may launch, which is 19.9 us of fixed cost a batch and 79% of a
1024-row one.

But whoever computed the matrix already touched every element, with the values
in registers. Producing the same four numbers per column there costs a
comparison or two and no extra memory traffic at all. And the result is
minuscule beside the payload it describes:

| | 65536x64 |
| --- | --- |
| matrix elements | 33.6 MB |
| per-column survey | 768 B |
| ratio | 1 : 43690 |

`Survey` is 12 bytes: two `int` exponents and two `bool` flags. The exponents
were `long long` at first, which made the struct 24 bytes rather than the 16 a
count of its fields suggests -- alignment padding, and two of the four fields
being one byte each. A double's true-ulp exponent lives in [-1074, 1023] and
its top in [-1073, 1024], so `int` carries six orders of magnitude more range
than the format can produce, and a `static_assert` now pins the size.

So it can simply be carried alongside the matrix. What that buys is out of
proportion to its size:

- **No survey, on either path.** The read disappears rather than being
  hidden or threaded.
- **No device round-trip.** The verdict is known before anything is launched,
  so the survey kernel, its drain and its 19.9 us go away — and with them the
  reason the device path ever needed to ask the host a question mid-batch.
- **Exact pre-sizing, before any processing begins.** `reserve_for` currently
  infers a column's exponent and width from one representative batch and a
  count. Given the survey for every batch up front, an accumulation matrix can be
  allocated with precisely the limb columns the whole stream will need, and
  then never rescale or widen at all — which is the expensive half of adaptive
  sizing, and the half that must happen between launches.

The cost is a trust boundary, and the failure modes are not symmetric.

**Over-reporting is harmless.** A column wider or lower than it needed to be
holds the right answer, in more memory, marginally slower.

**Under-reporting is silent, not an error**, which is the part worth designing
for. Both directions were demonstrated as bugs earlier in this project:

- *Exponent too high.* `shift = e - exponent` goes negative, the unsigned
  conversion makes `off` enormous, and the accumulate loop never executes —
  the value is dropped. This is the rescale bug in `eae167d`, where the shift
  itself landed correctly and the addend that triggered it vanished, giving 8
  where 9 was right.
- *Width too small.* `require_fit` derives `max_addend_bits` **from the
  survey**, so a survey that under-reports fools the check meant to catch it,
  and the sum wraps. This is `749707f`: twenty batches of 2^60 into a one-limb
  column gave 4611686018427387904 for 23058430092136939520.

Neither raises anything, and the reason is the decision the whole structure
rests on — the width is *derived* so that no entry can overflow, which is what
lets the carry chain run without a per-limb overflow test and be vectorized at
all. Take the derivation from an untrusted source and that guarantee leaves
with it.

Verifying the metadata up front is not an option: that is the survey. But
detecting the contradiction afterwards costs almost nothing, because the
accumulate already decomposes every element and the facts are in registers —
`shift < 0`, `off` beyond `nlimbs`, a non-finite significand, a carry out of
the top limb. Any of those means reality disagreed with what it was told. (The
first three shipped with the detector; the fourth, as a sign test on the top
limb only, came later and is what finally catches the `749707f` case above —
every addend fit, the sum did not.)
Reported through the mapped channel the occupancy figure already uses and
checked at the next readback, that turns a silent wrong answer into a loud one
without putting a second read on the critical path.

**The interface.** The metadata arrives as a vector of per-column surveys per
matrix — `vector<vector<Survey>>`, outer indexed by matrix and inner by column
— given to the constructor. The outer size *is* the count, so nothing separate
has to be passed to derive the headroom, and preallocation is offered only when
both are available. A caller who does not know how many matrices are coming, or
what is in them, gets exactly the code paths that exist today: survey the
pending batch, grow the limb columns from running statistics.

The constructor reduces it to one reservation per column — the minimum of the
minima for the exponent, the maximum of the maxima for the top, and the count
from the outer size — which is `reserve_for`'s arithmetic with the extents
supplied rather than inferred. Note what that reduction implies: the
accumulation matrix never needs to know *which* matrix it is being handed, because any
matrix inside the aggregate extents fits the reservation. Submitting them in
any order therefore falls out rather than having to be arranged, which suits a
structure whose exactness already makes accumulation order-independent.

Three consequences beyond the survey disappearing:

- **No rescaling and no widening, ever.** Both exist to react to a batch that
  did not fit; nothing will not fit. That is the expensive half of adaptive
  sizing and the half that has to happen between launches.
- **No running bookkeeping.** `require_fit` maintains `max_addend_bits` and
  `add_count` solely to derive a width that is now given, so the whole path
  shortens rather than merely speeding up.
- **The device path stops asking the host anything mid-batch**, which is what
  its 19.9 us of fixed cost is: a survey kernel and the round-trip carrying its
  verdict. That is 79% of a 1024-row batch, so small matrices gain most.

Two constraints the interface has to state rather than imply. Submissions must
not exceed the declared count, because the width bound holds for that many and
no more; fewer is fine, and merely over-allocated. And *any order* does not
mean *concurrently* — two batches read-modify-write the same limbs of the same
rows, so they still serialize, whatever order they arrive in.

`Survey` presently lives in `src/kernels.hpp` and would have to become public
vocabulary, shared between whoever produces a matrix and both accumulation matrices.
That argues for a small header of its own rather than exposing the kernel
surface around it.

**If the producer cannot supply it**, the fallback is to keep the survey on the
device but stop asking the host about it. Survey and accumulate launch back to
back, and the accumulate reads the verdict first: if the batch fits — the
common case, and the only one when `reserve_for` has done its job — it proceeds
with no host involvement. If it does not, the accumulate does nothing, leaving
the accumulation matrix untouched, and the host discovers it at the next check and
performs the rescale or widen it needs before resubmitting. That is optimiztic
execution rather than deferred validation: a rejected batch is skipped, not
half-applied, so there is no corrupted state to explain. The survey can also be
folded into the copy where one is happening anyway, so the data lands in device
memory and its extent is known from the same pass.

### Measured occupancy in place of a derived width (considered, dropped)

`fit_column` sizes a column from

    bits = max_addend_bits + ceil_log2(add_count + 1) + 1

which is conservative twice over: it assumes every addend was as large as the
largest ever seen, and that all of them accumulated in the same direction.
`max_addend_bits` is a high-water mark that never recedes, so after
`1e300 + 1 - 1e300` it stands near 1000 bits while the value occupies one.

The accumulate can report what the accumulation matrix actually holds, nearly free — it
has already decomposed every value, and the carry loop already knows the
highest limb it disturbed. Measured, the per-element statistics cost nothing
(-0.2% in the divergent case, where the kernel does most work per element);
the whole cost is the fixed election, 10% of a 4096-row batch and 0.3% of a
262144-row one.

Combined with the host survey of the *pending* batch, that answers exactly the
question the code cannot currently ask:

    fits = max(occupancy, pending_max_top - exponent) + 1 + 1 <= nlimbs * 64

This needs no deferred validation. The host survey inspects the pending input
before the accumulate runs, so non-finite values and out-of-range exponents are
still rejected before anything is written. An earlier version of this plan
assumed fusing the survey into the accumulate and accepting a corrupted
accumulation matrix on rejection; that was unnecessary, and the two statistics answer
different questions — the survey describes the input, the accumulate describes
the accumulation matrix.

The structural consequence is what makes it worth the election: it is what lets
the device widen at all. Widening is an error there today because the fit test
is conservative and made before the fact. With exact occupancy the host knows,
between launches, whether the limb arrays it has will hold the batch it is
about to send — and appending limb arrays disturbs nothing already allocated,
which is the property per-limb-position allocation exists to preserve.


## Symmetric matrices: store one triangle

`AccumulationMatrix(n, Uplo)` and `CudaAccumulationMatrix(n, Uplo)` keep only
`i >= j` (Lower) or `i <= j` (Upper). Column `j` then holds `n-j` entries or
`j+1`, and either way its stored rows stay contiguous — which is the whole
reason the triangle fits this layout at all. A warp still reads 32 consecutive
rows of one limb column as a single coalesced access; eight SIMD lanes still
walk one allocation. Only the column's *length* changes, and that was already
a per-column quantity in `LimbColumn`.

The memory halves, exactly: 268.4 to 134.3 MB of device allocation at
n = 4096. But the reason to store one copy rather than mirror two is not the
memory.

**Every column carries its own exponent and width.** In full storage `A(i,j)`
and `A(j,i)` therefore hold equal values in *different limbs*, and symmetry is
a property the data has to keep earning rather than one the structure
guarantees. Stored once it cannot drift, and the tests assert
`entry_limbs(i,j) == entry_limbs(j,i)` — an assertion full storage could not
pass.

Accumulation reads only the stored slice of an input matrix, so input traffic
halves with the storage. The input is taken to be symmetric and that is not
checked, because checking means reading the half this exists to avoid reading.

Threading needed a different split. A triangular column's work runs from `n`
entries down to 1, so an equal split by column *count* hands one worker most of
the matrix; `partition_columns` splits on cumulative stored entries instead,
which stays contiguous because the count is monotonic.

| n | threads | full | triangular | |
| --- | --- | --- | --- | --- |
| 512 | 1 | 0.223 | 0.451 | 2.02x |
| 2048 | 1 | 0.174 | 0.377 | 2.17x |
| 2048 | 8 | 0.501 | 0.843 | 1.68x |
| 1024 | GPU | 2.442 | 3.716 | 1.52x |
| 4096 | GPU | 1.261 | 1.817 | 1.44x |

Slightly *over* 2x single-threaded on the CPU rather than exactly 2x: at
n = 2048 full storage is 34.8 MB against 32 MiB of L3 while the triangle is
18.1 MB, so halving the footprint also buys back cache residency.

The GPU's 1.4-1.5x above is short of the CPU's 2x, and **the kernel is not
where it goes**. Timed on the device path with the accumulation matrix pre-sized, so
that nothing but the kernel is in the measurement, n = 2048:

| shape | stored | ns/stored |
| --- | --- | --- |
| full | 4194304 | 0.1230 |
| uniform `n/2 x n` | 2097152 | 0.1272 |
| Upper | 2098176 | 0.1273 |
| Lower | 2098176 | 0.1271 |

Every shape costs the same per stored element, within 3.5%. The triangle is
already getting ~1.94x out of the kernel; the shortfall in the 1.4-1.5x figure
is host-side work in `add_matrix_col_major`, which stages and surveys before it
launches, and not the device at all.

Two wrong explanations were recorded here before that was measured -- grid load
imbalance, then a ~470 ns per-block fixed cost -- and both are corrected in the
caveats below.

## Readback: the rounded double and its residual

`to_double(i, j, &residual)` also returns what the rounding discarded: exactly
`(stored value - returned double)`, itself correctly rounded.

The subtraction happens in the accumulation matrix's own fixed point, not in floating
point — forming `exact - hi` in doubles is precisely the cancellation this
container exists to avoid. The rounding step returns the exact reconstruction
of what it produced, `m * 2^scale` with `scale >= exp`, so `m << (scale - exp)`
lands back in the column's fixed point and the difference is one exact limb
subtraction using the same 128-bit-at-an-offset decomposition the accumulate
kernel uses.

Verified two ways: the suite re-accumulates each entry, subtracts the returned
double and requires what remains to equal the reported residual exactly; and
400 random multi-term sums were compared against exact rational arithmetic
outside this codebase, matching bit for bit including the sign of zero.

What it is worth, and what it is not. The pair carries ~106 bits against the
double's ~53 — measured worst relative error 2^-107 against 2^-53 over those
400 sums. But `lo` is a `double`, so it is a fixed +53 bits and *not* "the rest
of the value": columns here are routinely 128-192 bits wide, and everything
below the pair stays in the limbs. And the pair has to be kept unevaluated —
`hi + lo` in double arithmetic returned `hi` unchanged in 399 of those 400
cases, because `|lo| <= ulp(hi)/2` by construction.

## Multi-batch: folding K matrices into one pass

`add_matrices_col_major(b, count)` applies `count` matrices to each column in a
single pass. Identical results to a loop over `add_matrix_col_major` — exact
accumulation does not care about order or grouping — but one at a time each
batch reads and writes every limb it touches, while folded the limbs are read
once, all K addends applied in registers, and written once. Traffic per element
goes from ~40 bytes to `8 + 32/K`.

Column-major only. Folding row-major input would mean staging K columns per
worker, and writing and re-reading that staging is the traffic this exists to
avoid.

| 64 columns, Gelem/s, median of 3 | seq | K=2 | K=4 | K=8 |
| --- | --- | --- | --- | --- |
| 65536 rows, 8 threads | 1.14 | 1.84 | 2.65 | 3.01 |
| any shape, 1 thread | 0.57 | 0.69 | 0.59 | 0.49 |

**2.6x at the DRAM-bound shape, and K is a judgement rather than "as large as
possible".** Single-threaded only K=2 pays: there is no DRAM pressure to
relieve, and K input columns are K concurrent streams instead of one. Inside L3
on eight threads the run-to-run spread swamps the difference — sequential alone
varies 1.75 to 2.67 across runs — so nothing is claimed there.

Folding does **not** require pre-sizing, and the folded matrices need not be
all of them or come first. Without surveys the fold surveys its K columns and
fits the column once before adding any of them, so it settles the column itself
rather than requiring it settled; folded and single adds interleave in any
order, including a fold that rescales a column earlier single adds populated.

Separating the two levers at 65536x64 on eight threads:

| | adaptive | pre-sized |
| --- | --- | --- |
| one matrix at a time | 1.184 | 1.195 |
| folded, K=8 | 2.148 | 3.068 |

Folding alone is 1.81x; **pre-sizing alone is 1.01x**. The uplift is
accumulation matrix traffic, not allocation — both figures in the headline table were
already pre-sized, so neither reallocated at all. Pre-sizing then adds 1.43x on
top of folding, because an adaptive fold still reads its K columns to survey
them, and once folding has removed the dominant traffic those reads are what is
left to remove.

Columns wider than eight limbs fall back to one batch at a time. Not a
concession: the single-batch kernel stops at the first dead carry, while the
fold must write back every limb it loaded, so past that width folding would
move *more* memory.

## Implemented

**Bit-manipulation `decompose`** (`8565474`). Reads the IEEE-754 fields
directly; normal and denormal collapse to `e = max(biased, 1) - 1075`. Took
~28% off per-element cost end to end (12.3 → 8.8 ns narrow, 14.4 → 10.6 wide).

In isolation it is ~9x faster than the `frexp`/`ldexp` version, but that
overstates the win: `frexp`/`ldexp` were already inlined as builtins inside
`accumulate`, so the microbenchmark measured an out-of-line call the real code
never made.

**Peeled tails in the AVX-512 kernels.** Every loop used to compute a tail
mask per 8-row block and thread it through the body, though it is all-ones for
every block but the last. Peeling the final partial block out and passing a
constant mask to the rest folds the masking away entirely: the survey's pass-1
went 0.144 to 0.123 ns/elem (**14.7%**) and the accumulate 1.159 to 1.073
(**7.4%**) at four limbs.

Note what this was *not*: jagged row counts were never the cost. 4093 rows and
4096 rows measured identically both before and after, because the ragged block
is one iteration in 512. The cost was the per-block tax paid by every block for
the possibility of a tail. The survey takes the larger share because its pass-1
body is only ~8 instructions, so three extra ones are a large fraction of it.

This also removes the case for requiring row counts to be a multiple of 8, or
for padding the caller's column: with the mask out of the loop there is
essentially nothing left to buy, and a zero-padding contract on
`add_matrix_col_major` would fail silently on the last column.

### Threading across columns

`set_threads(n)` partitions columns across workers. Columns are independent in
storage and each is touched by exactly one worker, so the hot path takes no
locks and results are bit-identical to a serial run — same exponent, same
width, same limbs, asserted at 2, 4 and 8 threads and clean under
ThreadSanitizer, which is the check that matters and which ASan and UBSan do
not make.

The only state the columns shared was the staging buffer the row-major path
gathers each strided column through; that is now one per worker.

Workers are parked between calls rather than created per call. Spawning per
batch measured 2.73 Gelem/s at 4096x64 against 4.64 for a pool, because a
batch is only about a hundred microseconds and thread creation is a real
fraction of that.

## What limits each target

Both are now bound by moving bytes rather than by arithmetic, but at different
points, which is the whole shape of the comparison.

| | 4096x64 | 16384x64 | 65536x64 |
| --- | --- | --- | --- |
| CPU, 1 thread | 0.80 | 0.77 | 0.67 |
| CPU, 8 threads | **4.31** | **4.65** | 1.11 |
| GPU, host input | 3.11 | 3.26 | **3.23** |
| GPU, input already resident | 4.40 | 6.07 | 6.65 |

**The CPU is issue-bound while its working set fits in cache, and scales
almost linearly there** — 2.32 IPC at a 1% miss rate, 6.0x over eight cores.
At 65536x64 the input alone is 33.6 MB against 32 MiB of L3 and eight cores
queue on DRAM at about 58 GB/s of the 83 available, so threading buys 1.7x
instead of 6.

**The GPU is bandwidth-bound everywhere.** With input resident the kernel
moves about 48 bytes an element at 330 GB/s, which is 100% of what the card
streams in a plain copy loop. With input on the host it is PCIe-bound at 97%.

So the crossover is cache capacity, not compute: below it the CPU wins, above
it the GPU does, and the ratio at the top is just the bandwidth ratio.

## Remaining work, in order

Ping-pong input buffers, multi-batch folding on the CPU, and taking the survey
out of the library have all landed since the last revision of this list, along
with symmetric storage and readback residuals, which were not on it. Carry-save
at radix 52 came off it by being measured rather than by being done.

1. **Getting the input into device memory in the first place.** With the
   staging copy gone the host path is as good as it gets, and PCIe is what is
   left. Measured pre-sized, per batch, from the same data:

   | rows x 64 | staged (since removed) | pinned, in place | device-resident |
   | --- | --- | --- | --- |
   | 4096 | 93 us | 88 | 37 |
   | 16384 | 353 | 322 | 135 |
   | 65536 | 2278 | 1263 | **516** |

   Device-resident is 2.35-2.45x the best host path, consistently across
   shapes, and 516 us at 65536x64 is 0.123 ns an element -- exactly the
   kernel-only figure measured above. So device residency *is* the pure kernel
   cost, and every microsecond above it is the bus.

   That is the argument for a NIC writing into device memory over GPUDirect, or
   anything else that skips host memory: it is worth more than any host-side
   tuning left, because there is no host side left.

   `survey_matrix_col_major_device` exists for the same reason. A producer
   whose data never touches host memory could not previously describe it: the
   survey entry points were host-only, and the device survey kernel was private
   to `reserve_for_device`, which consumed its result and discarded it. Now the
   survey can be taken where the data is and sent ahead -- 768 bytes against
   33.6 MB, so it arrives long before the matrix does and the accumulation matrix is
   sized before the first byte of bulk data lands.

   Its `out` is a **device** pointer, and that is the interesting part of the
   signature. Device memory keeps the survey where a caller feeding it to
   something else on the device wants it, with no copy at all; a caller wanting
   it on the host passes a mapped pinned pointer and thereby says so. Returning
   it in host memory would have made that choice for everyone and hidden a copy
   inside a function whose whole point is that 768 bytes need not travel the
   same road as 33.6 MB.

   Writing the caller's buffer directly also removed the internal
   `DeviceSurvey`. It existed because `Survey`'s two bools are not
   atomic-friendly -- but only the two extents need atomics, and they are ints
   at offsets 0 and 4 with the struct's 12-byte size keeping every entry
   4-aligned. The flags need none: every block that sets one sets it true, so
   concurrent stores of the same value race benignly. The twin bought nothing
   but a conversion at every boundary.

   Still missing, and known: neither side has a survey entry point that covers
   only a symmetric accumulation matrix's *stored triangle*. A caller pre-sizing a
   triangular accumulation matrix has to walk `column_rows(j)` values from
   `column_first_row(j)` itself. The test suite hand-rolls exactly that helper,
   which is usually the sign it belongs in the library.
   `add_matrix_col_major_device` already accepts such input; what is missing is
   a way for the producer to deliver into a buffer the accumulation matrix will read --
   the device-side counterpart of `acquire_input`.

   The block reductions were since cut back to what actually needs one. Only
   the extents and the occupancy are reduced through a shared-memory tree;
   `any` is implied by whether the minimum is still its sentinel, and the
   non-finite and detector flags are written straight to their staging by
   whichever threads raise one -- which in a clean batch is none, so the
   ordinary case performs no atomic and needs no array. Shared memory halves,
   4096 to 2048 bytes for the survey and 2049 to 1025 for the accumulate, and
   throughput is unchanged to within noise.

   Occupancy does *not* change on this card, which is worth saying because it
   was the reason to expect a gain. Both kernels are register-bound: the
   survey at 22 registers is limited to 11 resident blocks and the accumulate
   at 46 to 5, while shared memory allowed 25 and 49 before the change. Halving
   it moves those to 50 and 99, and the binding constraint never moves. The
   simplification stands on being simpler; the headroom would only matter on a
   part with less shared memory a SM, or if the registers came down first.

   The per-block fixed cost, by contrast, is not worth chasing: measured at
   ~12.7 ns a block against ~0.12 ns an element, which is 2.5% of a block at
   4096 rows a column. It is 62% at 64 rows, but the whole launch there is
   20.8 us. Attributing it further: of that ~12 ns, the shared-memory
   reduction tree is ~6.3, the `__threadfence` and ticket ~2.8, the elected
   block's copy ~2.0, and the per-column atomics ~0.9.

2. ~~Multi-batch folding on the device.~~ **Landed.**
   `add_matrices_col_major_device` takes up to 16 device-resident matrices and
   applies all of them in one pass. Every matrix's addend lands in the same
   limbs of the same row, back to back, so L1 absorbs the repeats and DRAM sees
   one read and one writeback however many are folded:

   | K | us a batch | ns/elem | bytes/elem | implied GB/s | vs K=1 |
   | --- | --- | --- | --- | --- | --- |
   | 1 | 516 | 0.1230 | 40.0 | 325.2 | 1.00x |
   | 2 | 627 | 0.0747 | 24.0 | 321.3 | 1.65x |
   | 4 | 838 | 0.0499 | 16.0 | 320.3 | 2.46x |
   | 8 | 1326 | 0.0395 | 12.0 | 303.6 | **3.11x** |

   Implied bandwidth never leaves 303-325 GB/s against the card's 326.5, so the
   kernel stays bandwidth-bound throughout and the whole gain is moving fewer
   bytes. Predicted 1.67 / 2.50 / 3.33 from `8 + 32/K`, measured 1.65 / 2.46 /
   3.11. 25.3 Gelem/s at K=8, against 8.13 at K=1 and 6.65 for the previous
   best on this card.

   **Blocked over the input set, not held in registers**, and that was the
   decision worth getting right. A register array has to be indexed at compile
   time or it spills to local memory, so registers would mean templating the
   kernel on the limb count -- which columns of one matrix do not share, since
   each carries its own width, so a single launch could not serve them. Letting
   L1 do the reuse sidesteps the question: the working set is 256 threads x
   nlimbs x 8 bytes, 4 KB at two limbs against 128 KB of L1. Measure the cheap
   version before paying for the templated one.

   The shortfall at K=8 -- 303.6 GB/s against 325 -- is the eight input streams
   beginning to cost. The same effect was worth 2.75x on the CPU and is worth
   about 7% here, which is what a machine built to hide memory latency buys.

   **When folding is worth anything, stated plainly, because it usually is
   not.** Sixteen resident matrices is not a workload; cycling two or three
   buffers is. And at PCIe rates the accumulation matrix is not the bottleneck to begin
   with, so making it faster changes nothing:

   | 65536x64 | us a batch |
   | --- | --- |
   | producer filling a device buffer, H2D at 26.7 GB/s | 1255 |
   | accumulation matrix alone, K=1 | 516 |
   | two buffers rotating, producer on its own stream | 1301 |
   | two buffers rotating, producer on the null stream | 1776 |

   The accumulation matrix already has 2.4x of headroom under a PCIe-fed producer, and
   a rotation hides it to within 4% of the fill. Folding would take 516 us to
   166 and the batch would still cost 1256. It earns its keep only where input
   arrives faster than about 65 GB/s -- GPUDirect from a NIC, NVLink, or data
   generated on the device -- which is the same conclusion the GPUDirect item
   above reaches from the other direction.

   **The null-stream result is a usage trap worth naming.** This accumulation matrix's
   stream is a blocking stream, so it implicitly synchronizes with the legacy
   null stream. A producer using plain `cudaMemcpy` serializes against the
   accumulate and gets none of the overlap -- 1776 against 1301, which is the
   whole difference between a pipeline and a queue. The header says so at the
   entry point.

~~`reserve_column` issues a `cudaMemsetAsync` per limb position.~~ **Done.**
The zeroing is deferred and then done in one launch, which is 14x cheaper than
the calls it replaces -- 512 arrays took 827 us of `cudaMemsetAsync` against 59
us for a single kernel. Reserving a whole pre-sized accumulation matrix is where the
count is largest, so that is where it pays:

| shape | limb arrays | before | after | |
| --- | --- | --- | --- | --- |
| 65536x64, 2 limbs | 128 | 935 us | 563 | 1.66x |
| 65536x64, 6 limbs | 384 | 1575 | 567 | **2.78x** |
| 4096x256 | 512 | 2083 | 1175 | 1.77x |
| 16384x1024 | 2048 | 7435 | 3835 | 1.94x |

Worth knowing what is left: the allocations, not the zeroing. At 512 arrays
`cudaMallocAsync` costs 944 us against the zeroing's 59, so allocation is now
the whole of it. One slab carved into slices would fix that and would give up
the property the layout exists for -- widening appends without moving anything
-- so it is not obviously worth having.

~~Measured occupancy could relax the derived width bound.~~ **Dropped, not
done.** It is worth under a limb where cancellation has kept a column small,
against a width bound that has to stay conservative for the carry chain to
run without a per-limb overflow test at all -- which is the property the whole
accumulate path rests on. Trading that for a fraction of a limb is not a trade.

That leaves `column_occupancy` with no consumer inside the library; relaxing the
bound was what it was for. It stays as a diagnostic, because a caller can use it
to see whether the surveys they supplied were over-conservative, and because
removing it would not buy much: the kernel epilogue is 12.7 ns a block in total,
occupancy is roughly 4 ns of that, and at realistic column lengths that is under
1% of a block. The detector flags share the same reduction and elected block and
are load-bearing, so the machinery does not go away either way.

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
- An 18% speedup was attributed to interleaving two row blocks' dependency
  chains. Rebuilding with the unrolling removed but the same prepare/apply
  structure produced identical end-to-end throughput (1.111-1.124 against
  1.123-1.125 Gelem/s), so the gain came from the restructure and the
  attribution was wrong. In isolation the interleaving is worth ~7% at two
  limbs and ~10% at four, and nothing at one; 4x matches 2x within noise and
  8x regresses. Measure the thing you are claiming, not a proxy for it.
- The claim that nvcc rejects `__uint128_t` in device code was asserted from
  recollection, not measured, and is false on CUDA 12.9. It had been made the
  load-bearing argument for a 32-bit radix. Benchmarking then took the
  conclusion as well: 32-bit limbs are slower than 64-bit on both targets. Two
  rounds of plausible reasoning, both overturned the moment either was
  measured — check the toolkit against the toolkit, and the hardware against
  the hardware.

- Carry-save was recorded as a 2.8x GPU win, from a benchmark holding the
  accumulation matrix in registers and adding a full-width value. Bounding it against
  the real kernels says the GPU gains nothing at all — it is bandwidth-bound,
  so the arithmetic is free — while the CPU, dismissed here as unmeasured,
  gains 25-52%. The same claim was wrong about the size of the effect *and*
  about which machine it was on. A microbenchmark that changes where the
  accumulation matrix lives is not measuring the same question.
- The skew was claimed at 2.8x from a synthetic carry walk. Re-measured
  against the real kernel it is 5-9%, and only when the column is wide and its
  exponents diverge. The walk was bandwidth-bound; the real kernel is
  issue-bound, so the one number did not transfer to the other. A first attempt
  at re-measuring it let `aligned_alloc` place the arrays and produced
  non-reproducible swings in both directions — the allocator has to be taken
  out of the experiment before the effect is visible at all.

- The multi-batch fold was first written scalar, on the strength of the
  traffic argument alone, and lost everywhere -- a flat 3x deficit at every
  shape, 0.19 against 0.57 Gelem/s single-threaded. The argument was sound and
  the implementation gave up eight lanes to collect on it. Two measurements
  said what to do rather than guessing: thread scaling showed the mechanism
  working (the fold scaled 7.0x across eight cores where sequential managed
  2.4x), and running K=8 over eight distinct buffers against the same buffer
  eight times -- identical arithmetic, identical accumulation matrix traffic -- gave
  1.000 against 2.748 Gelem/s, so most of the remaining loss was the input
  *stream count*, not the fold. A traffic model that counts bytes and not
  streams will mispredict this.
- The fold's uplift was initially attributed to avoiding reallocation. It is
  not: both paths in the headline measurement were already pre-sized, so
  neither reallocated, and separating the levers gives 1.81x for folding alone
  against 1.01x for pre-sizing alone. When two changes ship together, measure
  the 2x2 before crediting either.

- The triangle's shortfall on the GPU got two wrong explanations before it got
  a right one, and the second was worse than the first because it came with a
  number. First: grid load imbalance, one column per y index. Refuted by a
  uniform `n/2 x n` matrix -- no triangle, no imbalance -- costing the same per
  stored element. Second: a ~470 ns fixed cost per block, obtained by sweeping
  rows at a fixed column count and reading off the intercept. That sweep timed
  `add_matrix_col_major`, which stages a host copy and runs a full host-side
  survey before it launches anything, so what was divided by the block count
  was mostly host work that has no block in it at all. Re-run on the device
  path with the accumulation matrix pre-sized, the intercept is **12.7 ns**, and an
  independent mimic of just the epilogue agrees at ~12. Then the real answer
  fell out: every shape costs the same per stored element, so the kernel was
  never the problem and the loss is host-side.

  The lesson is narrow and worth stating plainly: an intercept is only a fixed
  *per-block* cost if the only thing that varies with the sweep is per-block
  work. Timing an API call instead of a kernel puts host work in the intercept,
  where it looks exactly like a device constant.

The synthetic carry walk used throughout has none of the real work: no
decompose, no per-lane variable shifts, no masked selects. Treat anything still
resting on it as unsettled.

## Current baseline

Ryzen 7 7700X (8 cores, 32 MiB L3, AVX-512 double-pumped) and an RTX 3060
(12 GB, 192-bit, 328.8 GB/s measured in a copy loop, 26.7 GB/s over PCIe).
Gelem/s, exact accumulations per second, mixed signs.

| | 4096x64 | 16384x64 | 65536x64 |
| --- | --- | --- | --- |
| CPU, 1 thread | 0.80 | 0.77 | 0.67 |
| CPU, 8 threads | 4.31 | 4.65 | 1.11 |
| CPU, 8 threads, folded K=8 | — | — | **3.01** |
| GPU, host input | 3.11 | 3.26 | 3.23 |
| GPU, input resident | 4.40 | 6.07 | 6.65 |

The 65536x64 row is the one that moved. It was the DRAM-bound case and the
worst number in the table; folding eight matrices per pass takes it from 1.11
to 3.01, past the GPU's host-input rate. Nothing else in the table changed,
because nothing else was bound by accumulation matrix traffic.

Pre-sizing does not appear here because it is worth ~1% on the single-matrix
path. Where it shows is the device's small-batch fixed cost, which it removes
outright: 56.5 to 8.9 us for a 1024-row batch read from device memory, since
there is no longer a survey kernel or a round trip to wait for.

The CPU reaches 1.22 Gelem/s single-threaded on one-limb columns, which is the
narrowest case and the one the pre-threading figure of 1.122 was measured on.

Per-element costs behind those numbers, single-threaded CPU:

| case | ns/elem |
| --- | --- |
| 1-limb columns | 0.82 |
| 2-limb columns | 1.27 |
| 8 limbs, exponents spread 400 binades | 5.2 |

and on the device, 19.88 us fixed per batch plus 0.1456 ns/elem marginal,
the marginal rate being 6.87 Gelem/s or 330 GB/s — the card's full streaming
bandwidth.
