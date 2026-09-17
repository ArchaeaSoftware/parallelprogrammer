# Loading an identity index vector without touching memory

Materializing the identity permutation `[0, 1, 2, …, N-1]` into an AVX2/AVX-512
register using **only immediates encoded in the instruction stream** — no
`vmovdqa [rip+const]` load from `.rodata`.

## Why load an identity vector?

A lane-index ramp `[0, 1, 2, …, N-1]` is one of the most reused constants in
vector code, because it is the vector that answers "which lane am I?" — and a
surprising amount of SIMD work is really per-lane bookkeeping. It shows up as:

- **Gather/scatter indices.** `vpgatherdd` and friends take a vector of indices;
  the identity is the base from which strided or computed access patterns are
  built (`base + idx*stride`, `idx*stride + offset`).
- **Permutation control.** For the register-to-register lookup instructions
  `vpermd`/`vpermps`, the identity is the do-nothing permutation you start from
  and perturb to express reversals, rotations, or compaction.
- **Tail masking.** Comparing the ramp against a remaining-element count
  (`idx < n`) produces the predicate mask for the final partial vector of a loop,
  the standard way to handle counts that are not a multiple of the vector width.
- **Iota / sequence generation.** Any arithmetic sequence is the ramp scaled and
  shifted; positions, coordinates, and prefix-sum scaffolding all begin here.

Because it recurs so often — frequently once per loop iteration or once per call
in a hot kernel — computing an identity vector cheaply and *unobtrusively* actually
matters.

The obvious lowering of `_mm256_set_epi32(7,6,5,4,3,2,1,0)` is a single 32-byte
constant load and, as we shall see, compilers are prone to emitting that instruction.
It is the fastest option when L1 is warm and uncontested — but it consumes a load-port
slot and a cache line, which is exactly why we are here to examine alternatives.

**Environment:** AMD Ryzen 7 7700X (Zen 4) · g++ 9.4.0 · `-O2 -march=native` ·
Linux. Byte counts are from `objdump -d`; cycle figures are TSC ticks/iteration.

## Negative Results

When I first started exploring this question months ago, I wanted to formulate
the solutions in terms of *immediate loads*, keeping the constants away from the
load-port slots. Immediates are encoded directly into the instruction stream
consumed by the CPU. I could think of a few different ways to start with an
immediate and generate eight (8) 32-bit values in a vector register. The problem
is that even the x86-64 instruction set does not provide for 16- or 32-byte
immediates; one way to approach the problem is to start with a 32-bit immediate
whose 8-bit lanes contain our numbers, then promote to 16- and then 32-bit:

```
__attribute((noinline))
__m256i get_indices2()
{
    __m128i m8 = _mm_cvtsi32_si128( 0x03020100 );
    __m128i m16 = _mm_unpacklo_epi8( m8, _mm_setzero_si128() );
    __m128i m16hi = _mm_add_epi16( m16, _mm_set1_epi16( 4 ) );
    __m128i m32 = _mm_unpacklo_epi16( m16, _mm_setzero_si128() );
    __m128i m32hi = _mm_add_epi32( m32, _mm_set1_epi32( 4 ) );

    return _mm256_inserti128_si256( _mm256_castsi128_si256( m32 ), m32hi, 1 );
}
```

x86-64 does have a 64-bit immediate load instruction, so we *could* start by loading
a single GPR with `0x0706050403020100`. But AVX2 doesn't have ergonomic ways to 
promote the lanes in that half-register to 16- and then 32-bit precision. In particular,
`vpunpck*` instructions operate within each 128-bit lane of the 256-bit AVX2 registers!

For whatever reason, I tabled the question for a good long time, but [this tweet](https://x.com/pskocik/status/2100607837198053884?s=20) inspired 
me to resurrect the code base and sic Claude on it.

And, though I feel like I have almost an encyclopedic memory of the
AVX2 instructions (though it often feels like writing a legal brief - constantly
hitting the books to double-check case law, dot your i's and cross your t's,
because there are so many non-orthogonal features of the instruction set - like the 
half-register updates of `vpunpck*`), Claude's first reaction was to use an AVX2 
instruction I had never heard of.

---

## The key instruction: `vpmovzxbd reg, reg`

`vpmovzxbd` is usually thought of as a load-and-zero-extend. Its **register-to-
register** form (AVX2, `ymm, xmm`) zero-extends 8 bytes in the low half of an
XMM into 8 dwords across a YMM in a single uop — collapsing an entire
byte→word→dword unpack ladder into one instruction.

Coupled with the 64-bit move-immediate instruction, the identity vector load can
be accomplished in just three instructions:

### AVX2 — 8 lanes `[0..7]`

```asm
movabs $0x0706050403020100, %rax   ; 8 identity bytes, immediate in the stream
vmovq  %rax, %xmm0                 ; GPR -> XMM, no memory
vpmovzxbd %xmm0, %ymm0             ; 8 bytes -> 8 dwords, reg-to-reg
```

20 bytes of payload, zero memory references. The 8-byte `movabs` immediate is the
minimal way to inject 8 distinct values; there is no instruction that places a
64-bit immediate directly into a vector register, so the GPR round-trip is
needed.

With that, Claude immediately obsoleted every intrinsics- and inline assembly-
based solution I'd been benchmarking against each other.

### AVX-512 — 16 lanes `[0..15]`

`vpmovzxbd zmm, xmm` fans 16 source bytes into 16 dwords, but `movabs` injects
only 8 bytes at a time, so the 16-byte source must be assembled first:

```asm
movabs $0x0706050403020100, %rax
movabs $0x0f0e0d0c0b0a0908, %rdx
vmovq   %rax, %xmm0
vpinsrq $1, %rdx, %xmm0, %xmm0     ; xmm0 = bytes 0..15
vpmovzxbd %xmm0, %zmm0             ; 16 bytes -> 16 dwords
```

Because `[0..15]` is an affine ramp (`high 8 = low 8 + 8`), the second immediate
can be *derived* rather than injected. Several such variants were tried — see the
table. **Deriving does not beat two `movabs` on total code size**: the broadcast/
add/mask machinery costs more bytes than the 10 bytes saved. It does reduce the
*constant data* injected (8 bytes vs 16) at equal speed, if that is the metric
you care about.

---

## Byte counts

Payload = materialization instructions only (excludes the `endbr64` CET prologue
and `ret`). Full = including them.

```
+----------------------------------------+---------+---------+------+--------+--------------------+
| Variant                                | ISA     | Payload | Full | Instrs | Memory             |
+----------------------------------------+---------+---------+------+--------+--------------------+
| get_indices1  (compiler default)       | AVX2    |    10 B | 15 B |      3 | reads 32 B .rodata |
| get_indices11 (movabs+vmovq+vpmovzxbd) | AVX2    |    20 B | 25 B |      5 | none               |
| 512_intrin    (compiler default)       | AVX-512 |    10 B | 15 B |      3 | reads 16 B .rodata |
| 512_asm       (two movabs + vpinsrq)   | AVX-512 |    37 B | 42 B |      7 | none               |
| 512_addb      (vpaddb derive)          | AVX-512 |    40 B | 45 B |      9 | none               |
| 512_halves    (vinserti32x8)           | AVX-512 |    48 B | 53 B |     10 | none               |
| 512_bcastadd  (masked +8)              | AVX-512 |    52 B | 57 B |     11 | none               |
+----------------------------------------+---------+---------+------+--------+--------------------+
```

- **Smallest immediate-only, per ISA:** `get_indices11` (20 B) and `512_asm` (37 B).
- **Least constant data at equal speed (AVX-512):** `512_addb` — one 8-byte
  immediate plus a `0x08` broadcast, 3 bytes larger than `512_asm`.
- **Avoid:** `512_bcastadd` — biggest *and* slowest (the `kmov` + masked-add
  chain adds latency).

### The intrinsic trap

The compiler defeats the intent of the *intrinsic* forms. Both
`_mm256_set_epi32(...)` and `_mm512_cvtepu8_epi32(<constant>)` are folded to a
`.rodata` load (`vmovdqa`/`vpmovzxbd [rip], …`). **Only the inline-asm forms
guarantee the immediate stays in the instruction stream.** If you want a load,
let the compiler emit one; if you want immediates, you must write the asm.

---

## Benchmark 1 — warm & uncontested (`loadidx.cpp`)

Each function called through a pointer in a tight loop; dominated by call +
`volatile`-store overhead, so this measures little beyond "nothing is
pathological." Clocks/iteration:

```
+-----------------------------------------+-------------+
| Variant                                 | clocks/iter |
+-----------------------------------------+-------------+
| AVX2    get_indices1  (memory)          |        6.55 |
| AVX2    get_indices11 (immediate)       |        6.55 |
| AVX2    unpack-ladder variants 2-7      |        ~7.4 |
| AVX-512 512_intrin (memory)             |        4.13 |
| AVX-512 512_asm / 512_addb / 512_halves |        6.55 |
| AVX-512 512_bcastadd                    |        7.37 |
+-----------------------------------------+-------------+
```

The memory-load `512_intrin` looks best here — **because the constant sits in
warm, idle L1.** That advantage is an artifact of the microbenchmark, which is
the whole point of the next one.

## Benchmark 2 — under cache & load-port pressure (`contention.cpp`)

Identity vectors live in real instruction streams where the load ports and L1 are
already busy. This benchmark models that:

- **Forces re-materialization** of *both* variants every iteration (`asm
  volatile`), since a loop-invariant constant load would otherwise be hoisted out
  and cost nothing — and so would the immediate sequence.
- **Identical surrounding work:** 8 independent strided background loads per
  iteration into 8 accumulators (throughput-bound, not latency-bound), with the
  identity folded into one accumulator so its uops overlap.
- **Sweeps the working-set size** so the background traffic escalates from
  L1-resident to DRAM.

Only difference between columns: how the identity vector is produced.

```
+---------------+-------+----------+-----------+------------+-----------+
| Working set   |  none | mem-load | immediate | imm vs mem | winner    |
+---------------+-------+----------+-----------+------------+-----------+
| L1 (16 KB)    |  7.05 |     7.51 |      7.92 |      -5.4% | mem-load  |
| L1 (32 KB)    |  7.09 |     7.51 |      7.92 |      -5.5% | mem-load  |
| L2 (256 KB)   |  7.15 |     7.57 |      8.00 |      -5.7% | mem-load  |
| L2 (1 MB)     |  8.56 |     8.91 |      8.98 |      -0.7% | mem-load  |
| L3 (8 MB)     |  8.25 |     9.29 |      8.23 |     +11.4% | immediate |
| L3 (32 MB)    | 17.79 |    14.73 |     13.86 |      +5.8% | immediate |
| DRAM (128 MB) | 21.99 |    21.81 |     21.81 |      +0.0% | tie       |
+---------------+-------+----------+-----------+------------+-----------+
```

Crossover is stable across runs.

### Reading the result

- **Working set ≤ L2:** the constant load hits warm L1 and the load ports have
  slack, so the 1-uop load beats the immediate's 3 uops by ~5 %.
- **Working set L3-resident:** the background traffic now contends for L1 fill
  bandwidth and cache capacity. The memory load's extra cache line and load slot
  become a real cost; the immediate sequence touches neither and **wins by up to
  11 %.**
- **Pure DRAM footprint:** both are memory-latency bound; the identity's cost
  disappears into the stalls.

So the "fastest" way to load an identity vector depends entirely on what the rest
of the code is doing to the memory system. When the cache has better things to
do, immediates win.

Note that results were gathered on a rather old AMD Ryzen CPU (their first AVX-512
implementation), so YMMV on more recent microarchitectures.

---

## Building

```sh
# Microbenchmark + correctness (needs asmjit for the JIT variant)
g++ -O2 -march=native -I<asmjit>/src loadidx.cpp -lasmjit -o a.out && ./a.out

# Contention benchmark
g++ -O2 -march=native contention.cpp -o contention && ./contention [iterations]
```
