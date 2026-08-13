// Accumulates a batch of double matrices exactly and contrasts the result
// with naive double summation.
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "cbfp/column_accumulator.hpp"

namespace {

constexpr std::size_t kRows = 4;
constexpr std::size_t kCols = 4;
constexpr int kBatches = 20000;
constexpr long kMaxUlps = 1000000;

// How many representable doubles separate a from b.
long
ulp_distance(double a, double b)
{
    long ulps = 0;
    double walk = a;
    for (; ulps < kMaxUlps && walk != b; ++ulps) {
        walk = std::nextafter(walk, b);
    }
    return ulps;
}

}  // namespace

int
main()
{
    // Column 0: values near 1. Column 1: tiny. Column 2: huge. Column 3:
    // mixed magnitudes that cancel, which is where naive summation falls
    // apart.
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);

    std::vector<std::vector<double>> batches;
    batches.reserve(kBatches);
    for (int b = 0; b < kBatches; ++b) {
        std::vector<double> m(kRows * kCols);
        for (std::size_t i = 0; i < kRows; ++i) {
            m[i * kCols + 0] = unit(rng);
            m[i * kCols + 1] = std::ldexp(unit(rng), -80);
            m[i * kCols + 2] = std::ldexp(unit(rng), 80);
            m[i * kCols + 3] = std::ldexp(unit(rng), (b % 2) ? 60 : -60);
        }
        batches.push_back(std::move(m));
    }

    cbfp::ColumnBlockMatrix exact(kRows, kCols);
    exact.reserve_for(batches[0].data(), kBatches);

    std::vector<double> naive(kRows * kCols, 0.0);
    for (const auto& m : batches) {
        exact.add_matrix(m.data());
        for (std::size_t k = 0; k < kRows * kCols; ++k) naive[k] += m[k];
    }

    std::printf("%s\n", exact.describe().c_str());

    std::printf("%-4s %-24s %-24s %s\n", "cell", "exact (rounded)", "naive sum",
                "ulp error of naive");
    for (std::size_t i = 0; i < kRows; ++i) {
        for (std::size_t j = 0; j < kCols; ++j) {
            const double e = exact.to_double(i, j);
            const double n = naive[i * kCols + j];

            char cell[16];
            std::snprintf(cell, sizeof cell, "%zu,%zu", i, j);
            std::printf("%-4s %-24.17g %-24.17g %ld\n", cell, e, n,
                        ulp_distance(n, e));
        }
    }

    std::printf("\nexact value of cell 0,3:\n  %s\n",
                exact.to_exact_decimal(0, 3).c_str());
    return 0;
}
