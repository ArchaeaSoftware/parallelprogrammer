// contention.cpp
//
// The warm-L1 microbenchmark in loadidx.cpp says the compiler's memory-load
// (vmovdqa [rip], ymm) is as fast as anything: the constant sits in L1 and the
// load port is idle. But identity vectors live in *real* instruction streams,
// where the load ports and L1 are already doing useful work. This benchmark
// contrives that pressure and asks: does materializing from immediates (no load
// port, no cache line) beat loading the constant once the machine is busy?
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
#include <string.h>

// A real 32-byte in-memory constant for the "load" variant to read each iter.
alignas(64) static const uint32_t g_ident[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

// Number of independent background loads per iteration (load-port pressure).
// Each goes to its own accumulator so the loop is throughput- not latency-bound.
#ifndef BG
#define BG 8
#endif

enum Mode { MODE_MEM = 0, MODE_IMM = 1, MODE_NONE = 2 };

// Produce the identity vector. asm volatile prevents hoisting/CSE, so it is
// re-emitted every iteration exactly as un-hoistable real code would force.
template<int MODE>
static inline __m256i identity()
{
    __m256i idx;
    if (MODE == MODE_MEM) {
        // vmovdqa from a genuine memory constant -- consumes a load port + an
        // L1 line, exactly like the compiler's _mm256_set_epi32 lowering.
        asm volatile("vmovdqa %1, %0" : "=x"(idx) : "m"(g_ident[0]));
    } else if (MODE == MODE_IMM) {
        // movabs + vmovq + vpmovzxbd -- no memory, no load port.
        __m128i m8;
        asm volatile("movabs $0x0706050403020100, %%rax\n\t"
                     "vmovq %%rax, %0" : "=x"(m8) : : "rax");
        asm volatile("vpmovzxbd %1, %0" : "=x"(idx) : "x"(m8));
    } else {
        idx = _mm256_setzero_si256(); // baseline: background work only
    }
    return idx;
}

template<int MODE>
__attribute__((noinline))
static double run(const uint32_t* buf, size_t words, size_t N)
{
    const size_t mask = words - 1;          // words is a power of two
    // BG independent accumulators: no cross-iteration dependency, so the loop
    // is bound by load-port THROUGHPUT, not add latency. The identity vector is
    // folded into one accumulator -- its uops overlap the load-bound loop.
    __m256i acc[BG];
    for (int b = 0; b < BG; b++) acc[b] = _mm256_setzero_si256();
    size_t p = 0;

    uint64_t start = __rdtsc();
    for (size_t i = 0; i < N; i++) {
        __m256i d[BG];
        #pragma GCC unroll 16
        for (int b = 0; b < BG; b++) {
            asm volatile("vmovdqu %1, %0"
                         : "=x"(d[b])
                         : "m"(buf[(p + (size_t)b * 8) & mask]));
        }
        #pragma GCC unroll 16
        for (int b = 0; b < BG; b++)
            acc[b] = _mm256_add_epi32(acc[b], d[b]);

        p = (p + 8 * BG) & mask;

        // One extra "produce the identity vector" per iteration. In MEM mode
        // this is a 9th load (extra load-port pressure + an extra L1 line); in
        // IMM mode it is ALU/shuffle uops that overlap the load-bound loop.
        acc[0] = _mm256_add_epi32(acc[0], identity<MODE>());
    }
    uint64_t et = __rdtsc() - start;

    // Reduce and keep alive so nothing is eliminated.
    __m256i s = _mm256_setzero_si256();
    for (int b = 0; b < BG; b++) s = _mm256_add_epi32(s, acc[b]);
    uint32_t sink[8];
    _mm256_storeu_si256((__m256i*)sink, s);
    volatile uint32_t consume = sink[0] ^ sink[7];
    (void)consume;

    return (double)et / N;
}

struct Level { const char* name; size_t bytes; };

int main(int argc, char** argv)
{
    size_t N = (argc > 1) ? strtoull(argv[1], 0, 0) : 20000000ull;

    // Working-set sizes chosen relative to Zen 4: 32K L1D, 1M L2, 32M L3/iter.
    Level levels[] = {
        { "L1  (16 KB) ",  16ull << 10 },
        { "L1  (32 KB) ",  32ull << 10 },
        { "L2  (256 KB)", 256ull << 10 },
        { "L2  (1 MB)  ",   1ull << 20 },
        { "L3  (8 MB)  ",   8ull << 20 },
        { "L3  (32 MB) ",  32ull << 20 },
        { "DRAM(128 MB)", 128ull << 20 },
    };

    // One buffer big enough for the largest level.
    size_t maxwords = (levels[sizeof(levels)/sizeof(levels[0]) - 1].bytes) / sizeof(uint32_t);
    uint32_t* buf = (uint32_t*)aligned_alloc(64, maxwords * sizeof(uint32_t));
    for (size_t i = 0; i < maxwords; i++) buf[i] = (uint32_t)i;

    printf("BG loads/iter = %d, N = %zu, clocks/iter (lower is better)\n", BG, N);
    printf("%-13s  %8s  %8s  %8s   %s\n",
           "working set", "none", "mem-load", "immediate", "imm vs mem");
    printf("-------------------------------------------------------------------\n");

    for (auto& L : levels) {
        size_t words = L.bytes / sizeof(uint32_t);
        // warm/settle
        run<MODE_NONE>(buf, words, N/10);
        double t_none = run<MODE_NONE>(buf, words, N);
        double t_mem  = run<MODE_MEM >(buf, words, N);
        double t_imm  = run<MODE_IMM >(buf, words, N);
        double delta  = (t_mem - t_imm) / t_mem * 100.0;
        printf("%-13s  %8.2f  %8.2f  %8.2f   %+6.1f%%  %s\n",
               L.name, t_none, t_mem, t_imm, delta,
               t_imm < t_mem ? "immediate wins" : "mem-load wins");
    }

    free(buf);
    return 0;
}
