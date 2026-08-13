# SIMD and portability design record

**Status: largely implemented.** The limb-major layout, the skewed aligned
allocation and the AVX-512 survey and accumulate kernels have all landed, and
the measurements below are from the real kernels unless a section says
otherwise. What remains open is the radix (see the table — the CPU, not the
GPU, is where it would pay) and threading. This
header previously read "design, not implemented", which stopped being true at
`fdc8b75`.

Targets both AVX-512 on the CPU and a future CUDA implementation. CUDA turned
out to impose no special constraint on the representation — measurement found
the two targets want the same thing — so the abstraction is driven by the radix
choice below rather than by either machine's word size.

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
| Radix: 52-bit carry-save vs 64-bit canonical | **open, but on the CPU** |
| Whether limb arrays get separate allocations on CUDA | open |

Storage width is closed: 32-bit limbs were measured slower than 64-bit on both
targets, so `uint64_t` is the plan of record everywhere and the limb type does
not need to be a template parameter.

The **radix** is a separate question and remains open — but on the opposite
target from the one this document spent its time on.

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
does, but the reason has changed: it is no longer "the measurement has not
been run" but "the measurement says this is where the win is, and collecting
it means implementing carry-save properly.

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

**This benchmark held the accumulator in registers, and that is what makes it
misleading.** With no memory traffic the carry chain is the only cost, so
removing it is the whole game. The real kernel keeps the accumulator in global
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
selecting between two types, which forces `ColumnBlockMatrix` to become a
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
there is no type left to vary. Templating `ColumnBlockMatrix` would force
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

This argues for templates over a `CBFP_LIMB_BITS` macro: a macro forces one
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
The transfer alone costs 1258, so the accumulator arithmetic hides entirely
behind the bus instead of queueing after it.

| | 4096x64 | 16384x64 | 65536x64 |
| --- | --- | --- | --- |
| staged, buffer refilled each batch | 2.74 | 2.36 | 1.46 |
| borrowed, refilled each batch | 2.29 | 2.24 | 1.55 |

**These are not bus-bound, and an earlier revision of this section said they
were.** That claim came from a benchmark that filled one mapped buffer once and
resubmitted it, so the producer cost nothing; it read 3.11 / 3.26 / 3.23 and
was quoted as 97% of PCIe. It was also, not coincidentally, the unsafe reuse
pattern -- the same benchmark shape that corrupts 63% of entries once the
buffer is actually refilled.

Measured against the producer that a real caller has to run:

| | 33.6 MB batch |
| --- | --- |
| host writing the batch, any buffer kind | 1.6-1.7 ms (~20 GB/s) |
| host surveying it | ~0.5 ms |
| PCIe streaming it | 1.26 ms |

Writing a batch costs more than transferring it. The path is therefore bound by
the host producing and surveying the data, not by the bus, and reading in place
buys only the staging copy back -- which is real but is one of three comparable
costs rather than the whole of it. Zero-copy measures at 0.84-1.06x against
staging, not 2.2x.

What it does buy unconditionally is the device-side input buffer: 33.6 MB at
this shape, returned to the accumulator, which is what wants it.

Ordinary host memory still works and is staged through a mapped buffer, which
is why the middle column above improves least: that host copy then becomes the
limit. `allocate_input` hands back memory the device can read directly, and
`add_matrix_col_major` recognises it and skips staging.

The pinning measurement that justified the earlier design still stands and
still matters, because the staged path uses it: `cudaMemcpyAsync` from pageable
memory returns only after 3.92 ms of a 3.96 ms transfer, leaving no window to
survey in, while pinned returns in 0.01 ms of 2.51 and is 58% faster besides.

Validation stays fail-fast and costs nothing. The host survey finishes before
the accumulate is launched, so a batch that does not fit is rejected without
having touched the accumulator — which a device-side survey can only promise
by draining the pipeline to ask.

**Ping-pong buffers, and why the alternative was a trap.** Reading in place
means the kernel is still streaming a buffer after the call that submitted it
returned. The first version chose the in-place path by probing the pointer's
memory type, which gave one function two lifetime contracts: ordinary memory
was copied and could be reused at once, mapped memory could not, and nothing
said so. Switching allocator for speed silently corrupted 63% of entries.

`acquire_input` lends out a buffer the accumulator owns and blocks only until
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
  the next launch, so there is no H2D initialisation either. The election costs
  ~3.4 us and the PCIe write itself 0.64. Do **not** mark the device
  accumulator `volatile` — copied from the NVIDIA sample, it cost 5.71 us by
  defeating caching on every atomic, and turned a 1.12x win into a 0.90x loss.

The election serializes on one counter, so it is O(blocks); it wins below about
2000 blocks and loses above, which the grid cap keeps it under.

### Design intent: measured occupancy in place of a derived width

`fit_column` sizes a column from

    bits = max_addend_bits + ceil_log2(add_count + 1) + 1

which is conservative twice over: it assumes every addend was as large as the
largest ever seen, and that all of them accumulated in the same direction.
`max_addend_bits` is a high-water mark that never recedes, so after
`1e300 + 1 - 1e300` it stands near 1000 bits while the value occupies one.

The accumulate can report what the accumulator actually holds, nearly free — it
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
accumulator on rejection; that was unnecessary, and the two statistics answer
different questions — the survey describes the input, the accumulate describes
the accumulator.

The structural consequence is what makes it worth the election: it is what lets
the device widen at all. Widening is an error there today because the fit test
is conservative and made before the fact. With exact occupancy the host knows,
between launches, whether the limb arrays it has will hold the batch it is
about to send — and appending limb arrays disturbs nothing already allocated,
which is the property per-limb-position allocation exists to preserve.


## Implemented

**Bit-manipulation `decompose`** (`8565474`). Reads the IEEE-754 fields
directly; normal and subnormal collapse to `e = max(biased, 1) - 1075`. Took
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

Widening, rescaling, `reserve_for` and threading have all landed since the
last revision of this list.

1. **Ping-pong input buffers**, per "The host-input path" above. A correctness
   gap rather than an optimization: a caller cannot currently know when it is
   safe to refill a buffer the device is reading.
2. **Multi-batch accumulation, on the CPU.** Folding K batches into one
   accumulator read-modify-write takes traffic from ~40 bytes an element to
   8 + 32/K. Measured on the GPU as an upper bound — 8.25 to 27.2 Gelem/s at
   K=8, with implied bandwidth flat at ~330 GB/s throughout, which is what
   proves the kernel bandwidth-bound. The lever belongs on the **CPU**, whose
   65536x64 case is DRAM-bound at 1.11 Gelem/s; on the GPU's host-input path
   the bus carries the input regardless, so batching cannot reduce what
   crosses it.
3. **Carry-save at radix 52 on the CPU.** The largest measured win for the
   cache-resident sizes, where the CPU is issue-bound and every instruction
   removed is time removed: deleting the carry chain outright is worth 25-36%
   at 1-4 limbs and 45-52% at eight with divergent exponents. Expect less --
   at radix 52 a block's offsets span ~23% further, and that range drives the
   AVX-512 loop's trip count. De-risk as prescribed above: instantiate both
   radices and assert bit-identical `to_exact_decimal` over the corpus.
4. **The device path's 19.9 us fixed cost.** Only affects input that is
   already resident, and only small batches -- 79% of a 1024-row batch, 3% of
   a 65536-row one. The survey's verdict has to come back before the
   accumulate may launch, which drains the pipeline; surveying batch N+1
   while batch N accumulates would hide it, since the two do not depend on
   each other.

Smaller, known: `reserve_column` issues a `cudaMemsetAsync` per limb position,
`cols * nlimbs` of them, one-time at reserve rather than per batch. And
measured occupancy could relax the derived width bound where cancellation has
kept a column small, which is worth under a limb in the ordinary case and is
the lowest-value item here.

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
  accumulator in registers and adding a full-width value. Bounding it against
  the real kernels says the GPU gains nothing at all — it is bandwidth-bound,
  so the arithmetic is free — while the CPU, dismissed here as unmeasured,
  gains 25-52%. The same claim was wrong about the size of the effect *and*
  about which machine it was on. A microbenchmark that changes where the
  accumulator lives is not measuring the same question.
- The skew was claimed at 2.8x from a synthetic carry walk. Re-measured
  against the real kernel it is 5-9%, and only when the column is wide and its
  exponents diverge. The walk was bandwidth-bound; the real kernel is
  issue-bound, so the one number did not transfer to the other. A first attempt
  at re-measuring it let `aligned_alloc` place the arrays and produced
  non-reproducible swings in both directions — the allocator has to be taken
  out of the experiment before the effect is visible at all.

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
| GPU, host input | 3.11 | 3.26 | 3.23 |
| GPU, input resident | 4.40 | 6.07 | 6.65 |

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
