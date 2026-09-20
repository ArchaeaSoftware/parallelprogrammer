// Cost of the two structural adjustments on the CPU, per submission, as a
// column grows: a widen appends one limb-column per step (each addend is
// 2^64 times the last), a rescale lowers the exponent by 64 bits per step
// (each addend is 2^-64 times the last). Produces Table 2 of the paper.
//
// Usage: bench_adjust_cpu [rows] [cols] [steps]   (default 65536 16 12)
#include <truesum/accumulation_matrix.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "bench_common.hpp"

using namespace truesum;

namespace {

double
submit(AccumulationMatrix &m, std::vector<double> &buf, double v)
{
    for (auto &x : buf) x = v;
    const auto t0 = bench::Clock::now();
    m.add_matrix_col_major(buf.data());
    return bench::ms(t0, bench::Clock::now());
}

}  // namespace

int
main(int argc, char **argv)
{
    const std::size_t rows = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 65536;
    const std::size_t cols = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 16;
    const int steps = argc > 3 ? std::atoi(argv[3]) : 12;
    std::vector<double> buf(rows * cols);

    // Baseline: the same matrix again, so nothing is due but survey + add.
    AccumulationMatrix base(rows, cols);
    submit(base, buf, 1.0);
    double b0 = 0;
    for (int k = 0; k < 5; ++k) b0 += submit(base, buf, 1.0);

    AccumulationMatrix w(rows, cols), r(rows, cols);
    submit(w, buf, 1.0);
    submit(r, buf, 1.0);

    std::printf("%zu x %zu, one CPU thread\n", rows, cols);
    std::printf("submission with nothing to adjust: %.1f ms\n\n", b0 / 5);
    std::printf("%5s | %6s %9s | %6s %10s\n", "step", "limbs", "widen ms",
                "limbs", "rescale ms");
    for (int k = 1; k <= steps; ++k) {
        const double tw = submit(w, buf, std::ldexp(1.0, 64 * k));
        const std::size_t lw = w.column_limbs(0);
        const double tr = submit(r, buf, std::ldexp(1.0, -64 * k));
        const std::size_t lr = r.column_limbs(0);
        std::printf("%5d | %6zu %9.1f | %6zu %10.1f\n", k, lw, tw, lr, tr);
    }
    return 0;
}
