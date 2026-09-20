// Shared timing helpers for the benchmark drivers. Header-only so each
// driver stays a single translation unit.
#pragma once

#include <chrono>
#include <cstdio>
#include <vector>

namespace bench {

using Clock = std::chrono::steady_clock;

inline double
ms(Clock::time_point a, Clock::time_point b)
{
    return std::chrono::duration<double, std::milli>(b - a).count();
}

// Median of a small sample: the drivers report medians of three so a
// scheduling hiccup on one run does not become the number in the paper.
inline double
median3(double a, double b, double c)
{
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return b;
}

// Values for a column-major matrix spanning a moderate dynamic range: the
// case the library is designed for and the one Table 1 of the paper reports.
inline void
fill_moderate(std::vector<double> &m, unsigned long long seed)
{
    unsigned long long s = seed;
    for (auto &x : m) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        const double u = static_cast<double>(s >> 11) * (1.0 / 9007199254740992.0);
        // magnitudes across ~2^-8 .. 2^8, both signs
        const int e = static_cast<int>((s >> 3) % 17) - 8;
        double v = (u + 0.5);
        for (int k = 0; k < e; ++k) v *= 2.0;
        for (int k = 0; k > e; --k) v *= 0.5;
        x = (s & 1) ? -v : v;
    }
}

}  // namespace bench
