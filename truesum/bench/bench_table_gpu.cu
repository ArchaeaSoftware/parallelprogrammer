// Replicates the GPU rows of the throughput table: input read from mapped
// pinned host memory, input resident in device memory, and eight resident
// matrices batched per pass, for 4096x64, 16384x64 and 65536x64 column-major
// inputs. Reports billions of input values accumulated exactly per second,
// medians of three runs of ten submissions each.
#include <truesum/cuda_accumulation_matrix.hpp>
#include <truesum/survey.hpp>

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "bench_common.hpp"

using namespace truesum;

namespace {

constexpr std::size_t kCols = 64;
constexpr int kSubmissions = 10;
constexpr std::size_t kBatch = 8;

#define CUDA_CHECK(x)                                                        \
    do {                                                                     \
        cudaError_t err_ = (x);                                              \
        if (cudaSuccess != err_) {                                           \
            std::fprintf(stderr, "%s failed: %s\n", #x,                      \
                         cudaGetErrorString(err_));                          \
            std::exit(1);                                                    \
        }                                                                    \
    } while (0)

double
gelem_per_s(std::size_t elems, double ms)
{
    return static_cast<double>(elems) / (ms * 1e-3) / 1e9;
}

// Input in mapped pinned host memory, streamed over the bus by the kernel.
double
run_host_input(std::size_t rows)
{
    const std::size_t n = rows * kCols;
    std::vector<double> m(n);
    bench::fill_moderate(m, 1234);
    double *host = nullptr;
    CUDA_CHECK(cudaHostAlloc(&host, n * sizeof(double), cudaHostAllocMapped));
    for (std::size_t i = 0; i < n; ++i) host[i] = m[i];

    std::vector<Survey> sv(kCols);
    survey_matrix_col_major(sv.data(), m.data(), rows, kCols);
    const std::vector<std::vector<Survey>> surveys(kSubmissions, sv);

    double best[4];
    for (int rep = 0; rep < 4; ++rep) {
        CudaAccumulationMatrix acc(rows, kCols, surveys);
        acc.synchronize();
        const auto t0 = bench::Clock::now();
        for (int k = 0; k < kSubmissions; ++k) acc.add_matrix_col_major(host);
        acc.synchronize();
        best[rep] = bench::ms(t0, bench::Clock::now());
    }
    CUDA_CHECK(cudaFreeHost(host));
    return gelem_per_s(n * kSubmissions, bench::median3(best[1], best[2], best[3]));
}

// Input already resident in device memory, one matrix per submission.
double
run_resident(std::size_t rows)
{
    const std::size_t n = rows * kCols;
    std::vector<double> m(n);
    bench::fill_moderate(m, 1234);
    double *dev = nullptr;
    CUDA_CHECK(cudaMalloc(&dev, n * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(dev, m.data(), n * sizeof(double), cudaMemcpyHostToDevice));

    std::vector<Survey> sv(kCols);
    survey_matrix_col_major(sv.data(), m.data(), rows, kCols);
    const std::vector<std::vector<Survey>> surveys(kSubmissions, sv);

    double best[4];
    for (int rep = 0; rep < 4; ++rep) {
        CudaAccumulationMatrix acc(rows, kCols, surveys);
        acc.synchronize();
        const auto t0 = bench::Clock::now();
        for (int k = 0; k < kSubmissions; ++k) acc.add_matrix_col_major_device(dev);
        acc.synchronize();
        best[rep] = bench::ms(t0, bench::Clock::now());
    }
    CUDA_CHECK(cudaFree(dev));
    return gelem_per_s(n * kSubmissions, bench::median3(best[1], best[2], best[3]));
}

// Eight distinct resident matrices per pass over the accumulation matrix.
double
run_resident_batched(std::size_t rows)
{
    const std::size_t n = rows * kCols;
    std::vector<double *> dev(kBatch, nullptr);
    std::vector<double> m(n);
    std::vector<std::vector<Survey>> surveys;
    surveys.reserve(kSubmissions * kBatch);
    for (std::size_t b = 0; b < kBatch; ++b) {
        bench::fill_moderate(m, 777 + b);
        CUDA_CHECK(cudaMalloc(&dev[b], n * sizeof(double)));
        CUDA_CHECK(cudaMemcpy(dev[b], m.data(), n * sizeof(double), cudaMemcpyHostToDevice));
        std::vector<Survey> sv(kCols);
        survey_matrix_col_major(sv.data(), m.data(), rows, kCols);
        for (int k = 0; k < kSubmissions; ++k) surveys.push_back(sv);
    }
    std::vector<const double *> ptrs(dev.begin(), dev.end());

    double best[4];
    for (int rep = 0; rep < 4; ++rep) {
        CudaAccumulationMatrix acc(rows, kCols, surveys);
        acc.synchronize();
        const auto t0 = bench::Clock::now();
        for (int k = 0; k < kSubmissions; ++k) {
            acc.add_matrices_col_major_device(ptrs.data(), kBatch);
        }
        acc.synchronize();
        best[rep] = bench::ms(t0, bench::Clock::now());
    }
    for (double *p : dev) CUDA_CHECK(cudaFree(p));
    return gelem_per_s(n * kSubmissions * kBatch, bench::median3(best[1], best[2], best[3]));
}

}  // namespace

int
main()
{
    const std::size_t shapes[] = {4096, 16384, 65536};
    cudaDeviceProp prop{};
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::printf("GPU exact accumulation on %s, Gelem/s, %zu columns, medians of 3\n",
                prop.name, kCols);
    std::printf("%-30s", "rows");
    for (std::size_t r : shapes) std::printf(" %10zu", r);
    std::printf("\n");
    std::printf("%-30s", "host input (mapped pinned)");
    for (std::size_t r : shapes) std::printf(" %10.2f", run_host_input(r));
    std::printf("\n");
    std::printf("%-30s", "input resident");
    for (std::size_t r : shapes) std::printf(" %10.2f", run_resident(r));
    std::printf("\n");
    std::printf("%-30s", "resident, batched K=8");
    for (std::size_t r : shapes) std::printf(" %10.2f", run_resident_batched(r));
    std::printf("\n");
    return 0;
}
