/*
 *
 * bench.cpp
 *
 * Benchmark: fused operator chains against separate passes over memory.
 *
 * Build with: g++ -std=c++17 -O2 -mavx512f -mavx512dq -mavx512bw -mavx512vl -mfma bench.cpp
 *
 * Copyright (C) 2026 by Nicholas Wilt.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this 
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this
 *    list of conditions and the following disclaimer in the documentation and/or 
 *    other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "ops.hpp"
#include "fusion.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>

static double now_s() {
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static float* alloc_f(std::size_t n) {
    void* p = nullptr;
    if (posix_memalign(&p, 64, n * sizeof(float)) != 0) { std::perror("alloc"); std::exit(1); }
    return static_cast<float*>(p);
}

// Runtime-sourced constants: opaque to the optimiser, hoistable out of loops.
static volatile float g_seed = 1.0f;

static double checksum(const float* p, std::size_t n) {
    double s = 0;
    for (std::size_t i = 0; i < n; i += 997) s += p[i];
    return s;
}

template <int D>
struct Chain {
    Scale   s;
    FmaBias f;
    Clamp   c;
    Poly<D> p;
    Relu    r;
    Chain(float seed, const float* coeffs)
        : s(0.5f + seed * 0.001f),
          f(1.25f + seed * 0.001f, 0.125f),
          c(-8.0f, 8.0f),
          p(coeffs),
          r(0.0f) {}
};

template <std::size_t U, int D>
static void measure(std::size_t n, const float* in, float* out, float* tmp,
                    const float* coeffs, bool verify, const char* label) {
    Chain<D> ch(g_seed, coeffs);

    auto time_it = [&](auto&& fn) {
        // Warm up: first touch of the output pages, and let the core reach
        // a steady clock before anything is recorded.
        for (int w = 0; w < 3; ++w) fn();

        // Adaptive repetition: aim for ~40 ms of work per trial.
        int reps = 1;
        double dt = 0;
        for (;;) {
            double t0 = now_s();
            for (int r = 0; r < reps; ++r) fn();
            dt = now_s() - t0;
            if (dt > 0.04 || reps >= (1 << 22)) break;
            reps *= 2;
        }
        // Best of several trials: the minimum is the least contaminated by
        // scheduling, interrupts, and frequency excursions.
        double best = dt / reps;
        for (int trial = 1; trial < 5; ++trial) {
            double t0 = now_s();
            for (int r = 0; r < reps; ++r) fn();
            double d = (now_s() - t0) / reps;
            if (d < best) best = d;
        }
        return best;
    };

    double t_sep = time_it([&] {
        avx512_separate_passes(in, out, tmp, n, ch.s, ch.f, ch.c, ch.p, ch.r);
    });
    double sum_sep = checksum(out, n);

    double t_vec = time_it([&] {
        fused_avx512_vecmajor<U>(in, out, n, ch.s, ch.f, ch.c, ch.p, ch.r);
    });
    double sum_vec = checksum(out, n);

    double t_op = time_it([&] {
        fused_avx512_transform<U>(in, out, n, ch.s, ch.f, ch.c, ch.p, ch.r);
    });
    double sum_op = checksum(out, n);

    double t_op2 = time_it([&] {
        fused_avx512_opmajor_byval<U>(in, out, n, ch.s, ch.f, ch.c, ch.p, ch.r);
    });
    double sum_op2 = checksum(out, n);

    if (verify) {
        double tol = 1e-6 * std::fabs(sum_sep);
        if (std::fabs(sum_sep - sum_vec) > tol ||
            std::fabs(sum_sep - sum_op)  > tol ||
            std::fabs(sum_sep - sum_op2) > tol) {
            std::printf("  !! MISMATCH sep=%.6f vec=%.6f op=%.6f op2=%.6f\n",
                        sum_sep, sum_vec, sum_op, sum_op2);
        }
    }

    double gib = double(n) * sizeof(float) / (1 << 30);
    std::printf("%-8s %10.1f KB  %8.3f %8.3f %8.3f %8.3f  %6.2fx %6.2fx %6.2fx\n",
                label, double(n) * sizeof(float) / 1024.0,
                t_sep * 1e3, t_vec * 1e3, t_op * 1e3, t_op2 * 1e3,
                t_sep / t_vec, t_sep / t_op, t_sep / t_op2);
    (void)gib; (void)sum_op2;
}

int main(int argc, char** argv) {
    const std::size_t NMAX = 1u << 24;   // 64 MB of floats
    float* in  = alloc_f(NMAX);
    float* out = alloc_f(NMAX);
    float* tmp = alloc_f(NMAX);

    for (std::size_t i = 0; i < NMAX; ++i)
        in[i] = 0.25f + float((i * 2654435761u) >> 20) * (1.0f / 4096.0f);

    float coeffs[16];
    for (int i = 0; i < 16; ++i) coeffs[i] = 0.1f * float(i + 1);

    std::printf("chain: scale -> fma_bias -> clamp -> poly<D> -> relu\n");
    std::printf("times in ms per call; speedups are separate/fused\n\n");

    std::printf("=== size sweep, U=4, D=5 ===\n");
    std::printf("%-8s %13s  %8s %8s %8s %8s  %6s %6s %6s\n",
                "impl", "bytes", "separate", "vec/loop", "vec/pack", "opmaj/val",
                "vec", "pack", "op/val");
    for (std::size_t n : {std::size_t(1) << 12, std::size_t(1) << 16,
                          std::size_t(1) << 20, std::size_t(1) << 24})
        measure<4, 5>(n, in, out, tmp, coeffs, true, "sweep");

    std::printf("\n=== non-power-of-two sizes, U=4, D=5 ===\n");
    for (std::size_t n : {std::size_t(1000), std::size_t(4095),
                          std::size_t(65537), std::size_t(1000003)})
        measure<4, 5>(n, in, out, tmp, coeffs, true, "odd");

    std::printf("\n=== unroll sweep at n=2^20 (4 MB, L3), D=5 ===\n");
    measure<1, 5>(1u << 20, in, out, tmp, coeffs, true, "U=1");
    measure<2, 5>(1u << 20, in, out, tmp, coeffs, true, "U=2");
    measure<4, 5>(1u << 20, in, out, tmp, coeffs, true, "U=4");
    measure<8, 5>(1u << 20, in, out, tmp, coeffs, true, "U=8");

    std::printf("\n=== polynomial degree sweep at n=2^20, U=4 ===\n");
    measure<4, 3>(1u << 20, in, out, tmp, coeffs, true, "D=3");
    measure<4, 5>(1u << 20, in, out, tmp, coeffs, true, "D=5");
    measure<4, 7>(1u << 20, in, out, tmp, coeffs, true, "D=7");
    measure<4, 9>(1u << 20, in, out, tmp, coeffs, true, "D=9");
    measure<4, 11>(1u << 20, in, out, tmp, coeffs, true, "D=11");

    std::free(in); std::free(out); std::free(tmp);
    (void)argc; (void)argv;
    return 0;
}
