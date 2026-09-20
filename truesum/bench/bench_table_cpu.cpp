// Replicates the CPU rows of the throughput table: 1 thread, 8 threads, and
// 8 threads with eight matrices batched per pass, for 4096x64, 16384x64 and
// 65536x64 column-major inputs. Reports billions of input values accumulated
// exactly per second, medians of three runs of ten submissions each.
//
// Usage: bench_table_cpu [threads]     (default: 8)
#include <truesum/accumulation_matrix.hpp>
#include <truesum/survey.hpp>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "bench_common.hpp"

using namespace truesum;

namespace {

constexpr std::size_t kCols = 64;
constexpr int kSubmissions = 10;
constexpr std::size_t kBatch = 8;

double
gelem_per_s(std::size_t elems, double ms)
{
    return static_cast<double>(elems) / (ms * 1e-3) / 1e9;
}

// One matrix at a time, `threads` workers, on a pre-sized accumulation matrix:
// the surveys of every submission are supplied at construction, so no
// submission surveys, rescales or widens. This is the configuration the
// throughput table reports.
double
run_sequential(std::size_t rows, unsigned threads)
{
    std::vector<double> m(rows * kCols);
    bench::fill_moderate(m, 1234);
    std::vector<Survey> sv(kCols);
    survey_matrix_col_major(sv.data(), m.data(), rows, kCols);
    const std::vector<std::vector<Survey>> surveys(kSubmissions, sv);

    double best[4];
    for (int rep = 0; rep < 4; ++rep) {
        AccumulationMatrix acc(rows, kCols, surveys);
        acc.set_threads(threads);
        const auto t0 = bench::Clock::now();
        for (int k = 0; k < kSubmissions; ++k) acc.add_matrix_col_major(m.data());
        best[rep] = bench::ms(t0, bench::Clock::now());
    }
    const double t = bench::median3(best[1], best[2], best[3]);
    return gelem_per_s(rows * kCols * kSubmissions, t);
}

// Eight distinct matrices submitted in one pass over the accumulation matrix.
double
run_batched(std::size_t rows, unsigned threads)
{
    std::vector<std::vector<double>> ms(kBatch, std::vector<double>(rows * kCols));
    std::vector<const double *> ptrs(kBatch);
    for (std::size_t b = 0; b < kBatch; ++b) {
        bench::fill_moderate(ms[b], 777 + b);
        ptrs[b] = ms[b].data();
    }
    std::vector<std::vector<Survey>> surveys;
    surveys.reserve(kSubmissions * kBatch);
    for (std::size_t b = 0; b < kBatch; ++b) {
        std::vector<Survey> sv(kCols);
        survey_matrix_col_major(sv.data(), ms[b].data(), rows, kCols);
        for (int k = 0; k < kSubmissions; ++k) surveys.push_back(sv);
    }

    double best[4];
    for (int rep = 0; rep < 4; ++rep) {
        AccumulationMatrix acc(rows, kCols, surveys);
        acc.set_threads(threads);
        const auto t0 = bench::Clock::now();
        for (int k = 0; k < kSubmissions; ++k) {
            acc.add_matrices_col_major(ptrs.data(), kBatch);
        }
        best[rep] = bench::ms(t0, bench::Clock::now());
    }
    const double t = bench::median3(best[1], best[2], best[3]);
    return gelem_per_s(rows * kCols * kSubmissions * kBatch, t);
}

}  // namespace

int
main(int argc, char **argv)
{
    const unsigned threads = argc > 1 ? static_cast<unsigned>(std::atoi(argv[1])) : 8;
    const std::size_t shapes[] = {4096, 16384, 65536};

    std::printf("CPU exact accumulation, Gelem/s, %zu columns, medians of 3\n", kCols);
    std::printf("%-30s", "rows");
    for (std::size_t r : shapes) std::printf(" %10zu", r);
    std::printf("\n");

    {
        std::vector<double> probe(shapes[0] * kCols);
        bench::fill_moderate(probe, 1234);
        AccumulationMatrix a(shapes[0], kCols);
        a.reserve_for(probe.data(), kSubmissions);
        std::printf("column width after sizing: %zu limbs\n", a.column_limbs(0));
    }
    std::printf("%-30s", "1 thread");
    for (std::size_t r : shapes) std::printf(" %10.2f", run_sequential(r, 1));
    std::printf("\n");

    std::printf("%u threads%*s", threads, threads >= 10 ? 19 : 20, "");
    for (std::size_t r : shapes) std::printf(" %10.2f", run_sequential(r, threads));
    std::printf("\n");

    std::printf("%u threads, batched K=%zu%*s", threads, kBatch, threads >= 10 ? 6 : 7, "");
    for (std::size_t r : shapes) std::printf(" %10.2f", run_batched(r, threads));
    std::printf("\n");
    return 0;
}
