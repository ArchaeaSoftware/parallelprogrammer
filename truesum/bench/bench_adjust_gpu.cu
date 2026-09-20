// The device-side counterpart of bench_adjust_cpu: per-submission cost of a
// widen and a rescale as a column grows, plus the cost of the stream-ordered
// allocation on its own, which the paper's CUDA discussion cites.
//
// Usage: bench_adjust_gpu [rows] [cols] [steps]   (default 65536 16 12)
#include <truesum/cuda_accumulation_matrix.hpp>

#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "bench_common.hpp"

using namespace truesum;

namespace {

#define CUDA_CHECK(x)                                                        \
    do {                                                                     \
        cudaError_t err_ = (x);                                              \
        if (cudaSuccess != err_) {                                           \
            std::fprintf(stderr, "%s failed: %s\n", #x,                      \
                         cudaGetErrorString(err_));                          \
            std::exit(1);                                                    \
        }                                                                    \
    } while (0)

// cudaMallocAsync of one limb-column's worth of bytes, cold and then warm.
void
allocation_cost(std::size_t bytes)
{
    cudaStream_t s = nullptr;
    CUDA_CHECK(cudaStreamCreate(&s));
    std::vector<void *> p;
    double first = 0, rest = 0;
    for (int k = 0; k < 16; ++k) {
        void *q = nullptr;
        const auto t0 = bench::Clock::now();
        CUDA_CHECK(cudaMallocAsync(&q, bytes, s));
        CUDA_CHECK(cudaStreamSynchronize(s));
        const double t = bench::ms(t0, bench::Clock::now());
        if (0 == k) first = t; else rest += t;
        p.push_back(q);
    }
    for (void *q : p) CUDA_CHECK(cudaFreeAsync(q, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
    CUDA_CHECK(cudaStreamDestroy(s));
    std::printf("cudaMallocAsync of %zu KB: first %.3f ms (pool growth), "
                "then %.4f ms each\n", bytes >> 10, first, rest / 15);
}

}  // namespace

int
main(int argc, char **argv)
{
    const std::size_t rows = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 65536;
    const std::size_t cols = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 16;
    const int steps = argc > 3 ? std::atoi(argv[3]) : 12;
    const std::size_t n = rows * cols;

    cudaDeviceProp prop{};
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    std::printf("%zu x %zu on %s\n", rows, cols, prop.name);
    allocation_cost(rows * sizeof(double));

    double *host = nullptr;
    CUDA_CHECK(cudaHostAlloc(&host, n * sizeof(double), cudaHostAllocMapped));
    for (std::size_t i = 0; i < n; ++i) host[i] = 1.0;

    auto prime = [&](CudaAccumulationMatrix &m) {
        for (std::size_t j = 0; j < cols; ++j) m.reserve_column(j, 0, 64);
        m.add_matrix_col_major(host);
        m.synchronize();
    };
    auto submit = [&](CudaAccumulationMatrix &m, double v) {
        for (std::size_t i = 0; i < n; ++i) host[i] = v;
        const auto t0 = bench::Clock::now();
        m.add_matrix_col_major(host);
        m.synchronize();
        return bench::ms(t0, bench::Clock::now());
    };

    CudaAccumulationMatrix base(rows, cols), w(rows, cols), r(rows, cols);
    prime(base); prime(w); prime(r);
    double b0 = 0;
    for (int k = 0; k < 5; ++k) b0 += submit(base, 1.0);
    std::printf("submission with nothing to adjust: %.2f ms\n\n", b0 / 5);

    std::printf("%5s | %6s %9s | %6s %10s\n", "step", "limbs", "widen ms",
                "limbs", "rescale ms");
    for (int k = 1; k <= steps; ++k) {
        const double tw = submit(w, std::ldexp(1.0, 64 * k));
        const double tr = submit(r, std::ldexp(1.0, -64 * k));
        std::printf("%5d | %6zu %9.2f | %6zu %10.2f\n", k, w.column_limbs(0),
                    tw, r.column_limbs(0), tr);
    }
    CUDA_CHECK(cudaFreeHost(host));
    return 0;
}
