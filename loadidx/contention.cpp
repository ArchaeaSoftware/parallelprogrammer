// contention.cpp
//
// The warm-L1 microbenchmark in loadidx.cpp says the compiler's memory-load
// (vmovdqa [rip], reg) is as fast as anything: the constant sits in L1 and the
// load port is idle. But identity vectors live in *real* instruction streams,
// where the load ports and L1 are already doing useful work. This benchmark
// contrives that pressure and asks: does materializing from immediates (no load
// port, no cache line) beat loading the constant once the machine is busy?
//
// It sweeps both vector widths:
//   AVX2    256-bit, 8-lane  [0..7]   mem constant = 32 bytes
//   AVX-512 512-bit, 16-lane [0..15]  mem constant = 64 bytes (a full line)
//
// Two honesty requirements:
//   1. Force re-materialization every iteration. A loop-invariant constant load
//      is hoisted out of the loop by any decent compiler and costs nothing --
//      but so would the immediate sequence. We use asm volatile for BOTH so the
//      comparison is per-iteration cost, as it would be in un-hoistable code.
//   2. Identical surrounding work. Both loops do the same background memory
//      traffic (BG independent loads from a working set of size S) and the same
//      arithmetic; the ONLY difference is how the identity vector is produced.
//
// Sweep S across L1 / L2 / L3 / DRAM to see where the difference appears.

#include <x86intrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// Real in-memory constants for the "load" variant to read each iteration.
alignas(64) static const uint32_t g_ident256[8]  = { 0,1,2,3,4,5,6,7 };
alignas(64) static const uint32_t g_ident512[16] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };

// Number of independent background loads per iteration (load-port pressure).
// Each goes to its own accumulator so the loop is throughput- not latency-bound.
#ifndef BG
#define BG 8
#endif

enum Mode { MODE_MEM = 0, MODE_IMM = 1, MODE_NONE = 2 };

// ---- Per-width vector traits -----------------------------------------------
template<int W> struct V;

template<> struct V<256> {
    typedef __m256i t;
    static const int LANES = 8;                 // dwords per vector
    static inline t zero()          { return _mm256_setzero_si256(); }
    static inline t add(t a, t b)   { return _mm256_add_epi32(a, b); }
    static inline t loadu(const uint32_t* p) {
        t d; asm volatile("vmovdqu %1, %0" : "=x"(d) : "m"(*p)); return d;
    }
    static inline uint32_t drain(t v) {
        uint32_t s[8]; _mm256_storeu_si256((__m256i*)s, v); return s[0] ^ s[7];
    }
};

template<> struct V<512> {
    typedef __m512i t;
    static const int LANES = 16;
    static inline t zero()          { return _mm512_setzero_si512(); }
    static inline t add(t a, t b)   { return _mm512_add_epi32(a, b); }
    static inline t loadu(const uint32_t* p) {
        t d; asm volatile("vmovdqu64 %1, %0" : "=v"(d) : "m"(*p)); return d;
    }
    static inline uint32_t drain(t v) {
        uint32_t s[16]; _mm512_storeu_si512((void*)s, v); return s[0] ^ s[15];
    }
};

// ---- Identity producers (asm volatile => re-emitted every iteration) -------
template<int W, int MODE> static inline typename V<W>::t identity();

// AVX2, 8-lane [0..7]
template<> inline __m256i identity<256, MODE_MEM>() {
    __m256i idx;                                    // 32-byte load, one L1 line
    asm volatile("vmovdqa %1, %0" : "=x"(idx) : "m"(g_ident256[0]));
    return idx;
}
template<> inline __m256i identity<256, MODE_IMM>() {
    __m128i m8; __m256i idx;                        // movabs + vmovq + vpmovzxbd
    asm volatile("movabs $0x0706050403020100, %%rax\n\t"
                 "vmovq %%rax, %0" : "=x"(m8) : : "rax");
    asm volatile("vpmovzxbd %1, %0" : "=x"(idx) : "x"(m8));
    return idx;
}
template<> inline __m256i identity<256, MODE_NONE>() { return _mm256_setzero_si256(); }

// AVX-512, 16-lane [0..15]
template<> inline __m512i identity<512, MODE_MEM>() {
    __m512i idx;                                    // 64-byte load, a full line
    asm volatile("vmovdqa64 %1, %0" : "=v"(idx) : "m"(g_ident512[0]));
    return idx;
}
template<> inline __m512i identity<512, MODE_IMM>() {
    __m128i m8; __m512i idx;                        // two movabs + vpinsrq + vpmovzxbd
    asm volatile("movabs $0x0706050403020100, %%rax\n\t"
                 "movabs $0x0f0e0d0c0b0a0908, %%rdx\n\t"
                 "vmovq   %%rax, %0\n\t"
                 "vpinsrq $1, %%rdx, %0, %0" : "=x"(m8) : : "rax", "rdx");
    asm volatile("vpmovzxbd %1, %0" : "=v"(idx) : "v"(m8));
    return idx;
}
template<> inline __m512i identity<512, MODE_NONE>() { return _mm512_setzero_si512(); }

// ---- The measured kernel ----------------------------------------------------
template<int W, int MODE>
__attribute__((noinline))
static double run(const uint32_t* buf, size_t words, size_t N)
{
    typedef V<W> VT;
    typedef typename VT::t vec;
    const int L = VT::LANES;
    const size_t mask = words - 1;              // words is a power of two

    // BG independent accumulators: no cross-iteration dependency, so the loop
    // is bound by load-port THROUGHPUT, not add latency. The identity vector is
    // folded into one accumulator -- its uops overlap the load-bound loop.
    vec acc[BG];
    for (int b = 0; b < BG; b++) acc[b] = VT::zero();
    size_t p = 0;

    uint64_t start = __rdtsc();
    for (size_t i = 0; i < N; i++) {
        vec d[BG];
        #pragma GCC unroll 16
        for (int b = 0; b < BG; b++)
            d[b] = VT::loadu(&buf[(p + (size_t)b * L) & mask]);
        #pragma GCC unroll 16
        for (int b = 0; b < BG; b++)
            acc[b] = VT::add(acc[b], d[b]);

        p = (p + (size_t)L * BG) & mask;

        // One extra "produce the identity vector" per iteration. In MEM mode
        // this is an extra load (load-port pressure + an extra cache line); in
        // IMM mode it is ALU/shuffle uops that overlap the load-bound loop.
        acc[0] = VT::add(acc[0], identity<W, MODE>());
    }
    uint64_t et = __rdtsc() - start;

    volatile uint32_t consume = VT::drain(acc[0]);  // keep the work alive
    for (int b = 1; b < BG; b++) consume ^= VT::drain(acc[b]);
    (void)consume;

    return (double)et / N;
}

struct Level { const char* name; size_t bytes; };

template<int W>
static void sweep(const char* label, const uint32_t* buf,
                  const Level* levels, int nl, size_t N)
{
    printf("\n== %s ==\n", label);
    printf("%-13s  %8s  %8s  %8s   %s\n",
           "working set", "none", "mem-load", "immediate", "imm vs mem");
    printf("-------------------------------------------------------------------\n");
    for (int i = 0; i < nl; i++) {
        size_t words = levels[i].bytes / sizeof(uint32_t);
        run<W, MODE_NONE>(buf, words, N / 10);      // warm/settle
        double t_none = run<W, MODE_NONE>(buf, words, N);
        double t_mem  = run<W, MODE_MEM >(buf, words, N);
        double t_imm  = run<W, MODE_IMM >(buf, words, N);
        double delta  = (t_mem - t_imm) / t_mem * 100.0;
        printf("%-13s  %8.2f  %8.2f  %8.2f   %+6.1f%%  %s\n",
               levels[i].name, t_none, t_mem, t_imm, delta,
               t_imm < t_mem ? "immediate wins" : "mem-load wins");
    }
}

int main(int argc, char** argv)
{
    size_t N = (argc > 1) ? strtoull(argv[1], 0, 0) : 20000000ull;

    // Working-set sizes chosen relative to Zen 4: 32K L1D, 1M L2, 32M L3.
    Level levels[] = {
        { "L1  (16 KB) ",  16ull << 10 },
        { "L1  (32 KB) ",  32ull << 10 },
        { "L2  (256 KB)", 256ull << 10 },
        { "L2  (1 MB)  ",   1ull << 20 },
        { "L3  (8 MB)  ",   8ull << 20 },
        { "L3  (32 MB) ",  32ull << 20 },
        { "DRAM(128 MB)", 128ull << 20 },
    };
    const int nl = (int)(sizeof(levels) / sizeof(levels[0]));

    // Buffer big enough for the largest level, plus a pad so a full LANES-wide
    // load starting at any masked index stays in bounds.
    size_t maxwords = levels[nl - 1].bytes / sizeof(uint32_t);
    size_t alloc = maxwords + 32;
    uint32_t* buf = (uint32_t*)aligned_alloc(64, ((alloc * 4 + 63) & ~size_t(63)));
    for (size_t i = 0; i < alloc; i++) buf[i] = (uint32_t)i;

    printf("BG loads/iter = %d, N = %zu, clocks/iter (lower is better)\n", BG, N);

    sweep<256>("AVX2 (256-bit, 8-lane [0..7], 32B constant)",  buf, levels, nl, N);
    sweep<512>("AVX-512 (512-bit, 16-lane [0..15], 64B constant)", buf, levels, nl, N);

    free(buf);
    return 0;
}
