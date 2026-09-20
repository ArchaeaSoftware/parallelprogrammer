// bench_remote.cpp
//
// Portable version of contention.cpp: same experiment (memory-load vs immediate
// identity vector, under background load-port/cache pressure, swept across the
// working-set size), but the cache-level LABELS are detected from this machine's
// own /sys cache sizes so they are honest on any CPU. Sweep points are fixed
// powers of two (the index masking needs a power-of-two word count); each is
// labeled L1/L2/L3/DRAM from the detected sizes.
//
// Build:  g++ -O2 -march=native bench_remote.cpp -o bench_remote
// Run:    ./bench_remote [iterations]

#include <x86intrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

alignas(64) static const uint32_t g_ident256[8]  = { 0,1,2,3,4,5,6,7 };
alignas(64) static const uint32_t g_ident512[16] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };

#ifndef BG
#define BG 8
#endif

enum Mode { MODE_MEM = 0, MODE_IMM = 1, MODE_NONE = 2 };

template<int W> struct V;
template<> struct V<256> {
    typedef __m256i t;
    static const int LANES = 8;
    static inline t zero()        { return _mm256_setzero_si256(); }
    static inline t add(t a, t b) { return _mm256_add_epi32(a, b); }
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
    static inline t zero()        { return _mm512_setzero_si512(); }
    static inline t add(t a, t b) { return _mm512_add_epi32(a, b); }
    static inline t loadu(const uint32_t* p) {
        t d; asm volatile("vmovdqu64 %1, %0" : "=v"(d) : "m"(*p)); return d;
    }
    static inline uint32_t drain(t v) {
        uint32_t s[16]; _mm512_storeu_si512((void*)s, v); return s[0] ^ s[15];
    }
};

template<int W, int MODE> static inline typename V<W>::t identity();
template<> inline __m256i identity<256, MODE_MEM>() {
    __m256i idx; asm volatile("vmovdqa %1, %0" : "=x"(idx) : "m"(g_ident256[0])); return idx;
}
template<> inline __m256i identity<256, MODE_IMM>() {
    __m128i m8; __m256i idx;
    asm volatile("movabs $0x0706050403020100, %%rax\n\t"
                 "vmovq %%rax, %0" : "=x"(m8) : : "rax");
    asm volatile("vpmovzxbd %1, %0" : "=x"(idx) : "x"(m8));
    return idx;
}
template<> inline __m256i identity<256, MODE_NONE>() { return _mm256_setzero_si256(); }
template<> inline __m512i identity<512, MODE_MEM>() {
    __m512i idx; asm volatile("vmovdqa64 %1, %0" : "=v"(idx) : "m"(g_ident512[0])); return idx;
}
template<> inline __m512i identity<512, MODE_IMM>() {
    __m128i m8; __m512i idx;
    asm volatile("movabs $0x0706050403020100, %%rax\n\t"
                 "movabs $0x0f0e0d0c0b0a0908, %%rdx\n\t"
                 "vmovq   %%rax, %0\n\t"
                 "vpinsrq $1, %%rdx, %0, %0" : "=x"(m8) : : "rax", "rdx");
    asm volatile("vpmovzxbd %1, %0" : "=v"(idx) : "v"(m8));
    return idx;
}
template<> inline __m512i identity<512, MODE_NONE>() { return _mm512_setzero_si512(); }

template<int W, int MODE>
__attribute__((noinline))
static double run(const uint32_t* buf, size_t words, size_t N)
{
    typedef V<W> VT;
    typedef typename VT::t vec;
    const int L = VT::LANES;
    const size_t mask = words - 1;
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
        acc[0] = VT::add(acc[0], identity<W, MODE>());
    }
    uint64_t et = __rdtsc() - start;
    volatile uint32_t consume = VT::drain(acc[0]);
    for (int b = 1; b < BG; b++) consume ^= VT::drain(acc[b]);
    (void)consume;
    return (double)et / N;
}

// ---- Cache detection from sysfs --------------------------------------------
static size_t read_cache_size(const char* path) {
    FILE* f = fopen(path, "r"); if (!f) return 0;
    char buf[64] = {0};
    if (!fgets(buf, sizeof buf, f)) { fclose(f); return 0; }
    fclose(f);
    size_t v = strtoull(buf, 0, 10);
    const char* p = buf; while (*p >= '0' && *p <= '9') p++;
    if (*p == 'K' || *p == 'k') v <<= 10; else if (*p == 'M' || *p == 'm') v <<= 20;
    return v;
}
static void detect_caches(size_t& l1, size_t& l2, size_t& l3) {
    l1 = l2 = l3 = 0;
    for (int i = 0; i < 12; i++) {
        char lp[160], tp[160], sp[160];
        snprintf(lp, sizeof lp, "/sys/devices/system/cpu/cpu0/cache/index%d/level", i);
        snprintf(tp, sizeof tp, "/sys/devices/system/cpu/cpu0/cache/index%d/type", i);
        snprintf(sp, sizeof sp, "/sys/devices/system/cpu/cpu0/cache/index%d/size", i);
        FILE* f = fopen(lp, "r"); if (!f) continue;
        int lvl = 0; if (fscanf(f, "%d", &lvl) != 1) lvl = 0; fclose(f);
        char type[32] = {0};
        f = fopen(tp, "r"); if (f) { if (fscanf(f, "%31s", type) != 1) type[0] = 0; fclose(f); }
        size_t sz = read_cache_size(sp);
        if (lvl == 1 && type[0] == 'D') l1 = sz;
        else if (lvl == 2) l2 = sz;
        else if (lvl == 3) l3 = sz;
    }
}
static void cpu_model(char* out, size_t n) {
    out[0] = 0;
    FILE* f = fopen("/proc/cpuinfo", "r"); if (!f) return;
    char line[256];
    while (fgets(line, sizeof line, f)) {
        if (strncmp(line, "model name", 10) == 0) {
            char* c = strchr(line, ':');
            if (c) { c += 2; c[strcspn(c, "\n")] = 0; strncpy(out, c, n - 1); out[n - 1] = 0; }
            break;
        }
    }
    fclose(f);
}
static void fmt_size(size_t bytes, char* out, size_t n) {
    if (bytes >= (1u << 20)) snprintf(out, n, "%zu MB", bytes >> 20);
    else                     snprintf(out, n, "%zu KB", bytes >> 10);
}

struct Level { size_t bytes; char label[32]; };

template<int W>
static void sweep(const char* wlabel, const uint32_t* buf,
                  const Level* levels, int nl, size_t N)
{
    printf("\n== %s ==\n", wlabel);
    printf("%-13s  %8s  %8s  %8s   %s\n",
           "working set", "none", "mem-load", "immediate", "imm vs mem");
    printf("-----------------------------------------------------------------\n");
    for (int i = 0; i < nl; i++) {
        size_t words = levels[i].bytes / sizeof(uint32_t);
        run<W, MODE_NONE>(buf, words, N / 10);
        double t_none = run<W, MODE_NONE>(buf, words, N);
        double t_mem  = run<W, MODE_MEM >(buf, words, N);
        double t_imm  = run<W, MODE_IMM >(buf, words, N);
        double delta  = (t_mem - t_imm) / t_mem * 100.0;
        printf("%-13s  %8.2f  %8.2f  %8.2f   %+6.1f%%  %s\n",
               levels[i].label, t_none, t_mem, t_imm, delta,
               t_imm < t_mem ? "immediate" : "mem-load");
    }
}

int main(int argc, char** argv)
{
    size_t N = (argc > 1) ? strtoull(argv[1], 0, 0) : 40000000ull;

    size_t l1, l2, l3; detect_caches(l1, l2, l3);
    char model[256]; cpu_model(model, sizeof model);
    char b1[24], b2[24], b3[24];
    fmt_size(l1, b1, sizeof b1); fmt_size(l2, b2, sizeof b2); fmt_size(l3, b3, sizeof b3);

    printf("CPU: %s\n", model);
    printf("Caches: L1d=%s  L2=%s  L3=%s\n", b1, b2, b3);
    printf("BG loads/iter = %d, N = %zu, clocks/iter (lower is better)\n", BG, N);

    // Fixed power-of-two working sets; labeled from the detected cache sizes.
    static const size_t sizes[] = {
        16ull<<10, 32ull<<10, 128ull<<10, 512ull<<10,
        1ull<<20, 2ull<<20, 4ull<<20, 8ull<<20, 32ull<<20, 128ull<<20,
    };
    const int nl = (int)(sizeof(sizes)/sizeof(sizes[0]));
    Level levels[nl];
    for (int i = 0; i < nl; i++) {
        levels[i].bytes = sizes[i];
        const char* lvl = "DRAM";
        if (l1 && sizes[i] <= l1) lvl = "L1";
        else if (l2 && sizes[i] <= l2) lvl = "L2";
        else if (l3 && sizes[i] <= l3) lvl = "L3";
        char sz[24]; fmt_size(sizes[i], sz, sizeof sz);
        snprintf(levels[i].label, sizeof levels[i].label, "%-4s %s", lvl, sz);
    }

    size_t maxwords = sizes[nl - 1] / sizeof(uint32_t);
    size_t alloc = maxwords + 32;
    uint32_t* buf = (uint32_t*)aligned_alloc(64, ((alloc * 4 + 63) & ~size_t(63)));
    for (size_t i = 0; i < alloc; i++) buf[i] = (uint32_t)i;

    sweep<256>("AVX2 (256-bit, 8-lane [0..7], 32B constant)",  buf, levels, nl, N);
    sweep<512>("AVX-512 (512-bit, 16-lane [0..15], 64B constant)", buf, levels, nl, N);

    free(buf);
    return 0;
}
