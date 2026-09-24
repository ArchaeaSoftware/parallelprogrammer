// The CPU kernels against each other, per element, across column widths and
// exponent spreads. The kernel is chosen once per process, so run it twice --
// plain, and with TRUESUM_KERNEL=scalar -- and compare.
//
// Columns are pre-sized from surveys so only the accumulate is timed. "spread"
// is the range of true-ulp exponents in binades, uniform over the column; the
// clustered case sets the same width with one outlier per column while the
// bulk sits within 16 binades.
//
// This was written to evaluate a per-block dispatch in the AVX-512 kernel,
// which measured a 1.9x regression on clustered columns and is not in the
// library; docs/simd-design.md records what it measured and why it was
// rejected, along with the TRUESUM_SPAN knob that went with it. What the
// driver still shows is the crossover between the two kernels, which is why
// TRUESUM_KERNEL=scalar exists.
//
// Usage: bench_kernels     (built when TRUESUM_BUILD_BENCH is ON)
//
#include <truesum/accumulation_matrix.hpp>
#include <truesum/survey.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <vector>
using namespace truesum;
using Clock = std::chrono::steady_clock;

static std::uint64_t s = 88172645463325252ULL;
static std::uint64_t rnd() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }

// Odd 53-bit mantissas at exponents uniform over `spread` binades. If `outliers`,
// the width is set by a handful of values at the extremes while the bulk sits
// in a 16-binade cluster: wide column, clustered exponents.
static void fill(std::vector<double> &m, int spread, bool outliers)
{
    for (std::size_t i = 0; i < m.size(); ++i) {
        const double mant = static_cast<double>((rnd() >> 11) | 1);  // odd, < 2^53
        int e;
        if (outliers) e = (i % 4096 == 0) ? spread : static_cast<int>(rnd() % 17);
        else e = spread ? static_cast<int>(rnd() % (spread + 1)) : 0;
        const double v = std::ldexp(mant, e - 52 - spread / 2);  // centered, stays finite
        m[i] = (rnd() & 1) ? -v : v;
    }
}

static double ns_per_elem(std::size_t rows, std::size_t cols, int spread,
                          bool outliers, int subs, bool batch, std::size_t *limbs,
                          unsigned threads = 1)
{
    const int k = batch ? 8 : 1;
    std::vector<std::vector<double>> bufs(k, std::vector<double>(rows * cols));
    std::vector<std::vector<Survey>> surveys;
    std::vector<const double *> ptrs;
    for (auto &b : bufs) {
        fill(b, spread, outliers);
        std::vector<Survey> sv(cols);
        survey_matrix_col_major(sv.data(), b.data(), rows, cols);
        for (int i = 0; i < subs; ++i) surveys.push_back(sv);
        ptrs.push_back(b.data());
    }
    const std::vector<double> &m = bufs[0];
    double best = 1e30;
    for (int rep = 0; rep < 5; ++rep) {
        AccumulationMatrix acc(rows, cols, surveys);
        acc.set_threads(threads);
        const auto t0 = Clock::now();
        for (int i = 0; i < subs; ++i) {
            if (batch) acc.add_matrices_col_major(ptrs.data(), k);
            else acc.add_matrix_col_major(m.data());
        }
        const double t = std::chrono::duration<double, std::nano>(Clock::now() - t0).count();
        if (rep > 0 && t < best) best = t;    // rep 0 warms the pages
        *limbs = acc.column_limbs(0);
    }
    return best / (static_cast<double>(rows) * cols * subs * k);
}

static double survey_ns(std::size_t rows, std::size_t cols, int spread)
{
    std::vector<double> m(rows * cols);
    fill(m, spread, false);
    std::vector<Survey> sv(cols);
    double best = 1e30;
    for (int rep = 0; rep < 7; ++rep) {
        const auto t0 = Clock::now();
        for (int i = 0; i < 20; ++i) survey_matrix_col_major(sv.data(), m.data(), rows, cols);
        const double t = std::chrono::duration<double, std::nano>(Clock::now() - t0).count();
        if (rep > 0 && t < best) best = t;
    }
    return best / (static_cast<double>(rows) * cols * 20);
}

int main()
{
    AccumulationMatrix probe(8, 1);
    std::printf("kernel: %s\n", probe.describe().find("avx512") != std::string::npos ? "avx512" : "scalar");
    struct Case { const char *name; std::size_t rows, cols; int spread; bool out, batch; int subs; unsigned thr; };
    const Case cases[] = {
        {"spread 650",                    4096, 16, 650, false, false,  30, 1},
        {"spread 950",                    4096, 16, 950, false, false,  20, 1},
        {"spread 1900",                   4096, 16,1900, false, false,  15, 1},
        {"spread 0",                      4096, 16,   0, false, false, 200, 1},
        {"spread 40",                     4096, 16,  40, false, false, 200, 1},
        {"spread 150",                    4096, 16, 150, false, false, 100, 1},
        {"spread 400",                    4096, 16, 400, false, false,  60, 1},
        {"8 limbs, clustered + outlier",  4096, 16, 400, true,  false,  60, 1},
        {"batched K=8 distinct, spr 40",  4096, 16,  40, false, true,   25, 1},
        {"65536x64 spr 40, 1 thread",    65536, 64,  40, false, false,   6, 1},
        {"65536x64 spr 40, 6 threads",   65536, 64,  40, false, false,  10, 6},
    };
    for (const Case &c : cases) {
        std::size_t limbs = 0;
        const double ns = ns_per_elem(c.rows, c.cols, c.spread, c.out, c.subs, c.batch, &limbs, c.thr);
        std::printf("%-30s %2zu limbs  %.3f ns/elem\n", c.name, limbs, ns);
    }
    std::printf("%-30s %8s  %.3f ns/elem\n", "survey, spread 40", "", survey_ns(4096, 16, 40));
    std::printf("%-30s %8s  %.3f ns/elem\n", "survey, spread 40, DRAM", "", survey_ns(65536, 64, 40));
    return 0;
}
