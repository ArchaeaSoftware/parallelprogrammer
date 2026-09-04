/*
 *
 * bench_unrolled_binary_search.cpp
 *
 * Benchmarks for unrolled_binary_search.hpp: the shipped search against
 * std::lower_bound and a plain loop, and the four ways of writing the
 * unrolled descent that the accompanying article compares.
 *
 * Build with: g++ -std=c++17 -O2 -march=native bench_unrolled_binary_search.cpp
 *
 * Copyright (C) 2026 by Archaea Software, LLC.
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

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

#include "unrolled_binary_search.hpp"

#if defined(__clang__)
#  define BENCH_CC "clang " __clang_version__
#elif defined(__GNUC__)
#  define BENCH_CC "g++ " __VERSION__
#else
#  define BENCH_CC "unknown compiler"
#endif

// ---------------------------------------------------------------------------
// Query set: half the values are present, half fall between two elements, and
// the order is shuffled. The shuffle matters more than it looks. The libstdc++
// std::lower_bound is branchy and g++ leaves the branches in, so a query stream
// the predictor can learn costs 6 ns per lookup while a shuffled one costs 30.
// libc++ compiles the same source to conditional moves and barely notices. Any
// measurement of a search over sorted data has to say which one it did.
//
// The count is fixed rather than proportional to N, so that the queries occupy
// the same 32 KB whatever the array size and the numbers compare like with like.
// ---------------------------------------------------------------------------
template <std::size_t N>
static std::vector<int32_t> make_queries(const std::array<int32_t, N>& hay) {
    constexpr std::size_t kQueries = 4096;
    std::vector<int32_t> q;
    q.reserve(2 * kQueries);
    for (std::size_t i = 0; i < kQueries; ++i) {
        q.push_back(hay[i % N]);        // present
        q.push_back(hay[i % N] + 1);    // absent: the array holds even values only
    }
    std::mt19937 rng(12345);
    std::shuffle(q.begin(), q.end(), rng);
    return q;
}

// Accumulating the results defeats dead code elimination; the total is printed
// at the end so that nothing in the timed loop can be optimised away.
static std::ptrdiff_t sink = 0;

template <typename F>
static double time_ns(const std::vector<int32_t>& q, F f) {
    const int trials = 5, reps = 200;
    for (auto x : q) sink += f(x);                 // warm the caches
    double best = 1e30;
    for (int t = 0; t < trials; ++t) {
        auto t0 = std::chrono::steady_clock::now();
        for (int r = 0; r < reps; ++r)
            for (auto x : q) sink += f(x);
        auto t1 = std::chrono::steady_clock::now();
        double ns = std::chrono::duration<double, std::nano>(t1 - t0).count()
                  / static_cast<double>(reps) / static_cast<double>(q.size());
        if (ns < best) best = ns;
    }
    return best;
}

template <std::size_t N>
static std::ptrdiff_t lower_bound_search(const std::array<int32_t, N>& a, int32_t t) {
    auto p = std::lower_bound(a.begin(), a.end(), t);
    return (p != a.end() && *p == t) ? p - a.begin() : -1;
}

// ---------------------------------------------------------------------------
// The five formulations of the unrolled descent that the article walks through,
// plus the loop they are meant to improve on and the library function they are
// measured against. All are marked noinline, so that each is reached through
// the same call shape.
// ---------------------------------------------------------------------------

// As first generated: the segment length is a runtime argument and the target
// arrives by const reference, so a scalar has to be materialised in memory.
template <typename T, std::size_t... Steps>
__attribute__((noinline)) static std::ptrdiff_t d_generated(
        const T* arr, const T& target, std::size_t len,
        std::integer_sequence<std::size_t, Steps...>) {
    std::size_t idx = 0;
    ((idx += (idx + Steps < len && arr[idx + Steps] <= target ? Steps : 0)), ...);
    return arr[idx] == target ? static_cast<std::ptrdiff_t>(idx) : -1;
}

#define DESCENT(name, STMT)                                                          \
    template <std::size_t Len, typename T, std::size_t... Steps>                     \
    __attribute__((noinline)) static std::ptrdiff_t name(                            \
            const T* arr, T target, std::integer_sequence<std::size_t, Steps...>) {  \
        std::size_t idx = 0;                                                         \
        STMT                                                                         \
        return arr[idx] == target ? static_cast<std::ptrdiff_t>(idx) : -1;           \
    }

DESCENT(d_templated, ((idx += (idx + Steps < Len && arr[idx + Steps] <= target ? Steps : 0)), ...);)
DESCENT(d_unguarded, ((idx += (arr[idx + Steps] <= target ? Steps : 0)), ...);)
DESCENT(d_select,    ((idx = (arr[idx + Steps] <= target ? idx + Steps : idx)), ...);)
DESCENT(d_lambda,    auto probe = [&](auto step) {
                         if (arr[idx + step.value] <= target) idx += step.value;
                     };
                     (probe(std::integral_constant<std::size_t, Steps>{}), ...);)
#undef DESCENT

template <std::size_t Len, typename T>
__attribute__((noinline)) static std::ptrdiff_t d_loop(const T* a, T t) {
    std::size_t idx = 0;
    for (std::size_t step = Len / 2; step > 0; step >>= 1)
        if (a[idx + step] <= t) idx += step;
    return a[idx] == t ? static_cast<std::ptrdiff_t>(idx) : -1;
}

template <std::size_t N>
__attribute__((noinline)) static std::ptrdiff_t d_lower(const std::array<int32_t, N>& a, int32_t t) {
    return lower_bound_search(a, t);
}

template <std::size_t N>
static bool table() {
    std::array<int32_t, N> hay;
    for (std::size_t i = 0; i < N; ++i) hay[i] = static_cast<int32_t>(i * 2);
    auto q = make_queries(hay);
    const int32_t* a = hay.data();
    constexpr auto steps = make_probe_steps<N>();

    for (auto x : q) {
        std::ptrdiff_t ref = lower_bound_search(hay, x);
        if (d_generated(a, x, N, steps)  != ref || d_templated<N>(a, x, steps) != ref ||
            d_unguarded<N>(a, x, steps)  != ref || d_select<N>(a, x, steps)    != ref ||
            d_lambda<N>(a, x, steps)     != ref || d_loop<N>(a, x)             != ref ||
            d_lower(hay, x)              != ref) {
            std::printf("  MISMATCH at N=%zu, target=%d\n", N, x);
            return false;
        }
    }

    std::printf("descent formulations, N = %zu, out of line\n\n", N);
    struct { const char* name; double ns; } rows[] = {
        { "variadic, as generated",               time_ns(q, [&](int32_t x){ return d_generated(a,x,N,steps); }) },
        { "...length templated, target by value", time_ns(q, [&](int32_t x){ return d_templated<N>(a,x,steps); }) },
        { "...redundant test also dropped",       time_ns(q, [&](int32_t x){ return d_unguarded<N>(a,x,steps); }) },
        { "...conditional update as a select",    time_ns(q, [&](int32_t x){ return d_select<N>(a,x,steps); }) },
        { "...conditional update as an if",       time_ns(q, [&](int32_t x){ return d_lambda<N>(a,x,steps); }) },
        { "plain loop",                           time_ns(q, [&](int32_t x){ return d_loop<N>(a,x); }) },
        { "std::lower_bound",                     time_ns(q, [&](int32_t x){ return d_lower(hay,x); }) },
    };
    for (auto& r : rows) std::printf("  %-38s %5.2f\n", r.name, r.ns);
    return true;
}

// ---------------------------------------------------------------------------
// The shipped entry point as a caller would actually use it: inlined, at sizes
// that are and are not powers of two.
// ---------------------------------------------------------------------------
template <std::size_t N>
static bool shipped() {
    std::array<int32_t, N> hay;
    for (std::size_t i = 0; i < N; ++i) hay[i] = static_cast<int32_t>(i * 2);
    auto q = make_queries(hay);

    for (auto x : q)
        if (unrolled_binary_search(hay, x) != lower_bound_search(hay, x)) {
            std::printf("  MISMATCH: unrolled_binary_search, N=%zu, target=%d\n", N, x);
            return false;
        }

    double lb = time_ns(q, [&](int32_t x) { return lower_bound_search(hay, x); });
    double un = time_ns(q, [&](int32_t x) { return unrolled_binary_search(hay, x); });
    std::printf("  N = %-5zu  unrolled %5.2f ns   std::lower_bound %6.2f ns   (%.1fx)\n",
                N, un, lb, lb / un);
    return true;
}

int main() {
    std::printf("%s\n", BENCH_CC);
    std::printf("nanoseconds per lookup, best of five trials after warm-up\n"
                "sorted arrays of even int32_t; half the queries absent; order shuffled\n\n");

    bool ok = table<1024>();

    std::printf("\nshipped entry point, inlined at the call site\n\n");
    ok = shipped<1000>() && ok;
    ok = shipped<1024>() && ok;
    ok = shipped<4096>() && ok;

    std::printf("\n%s (checksum %td)\n",
                ok ? "all variants agree with std::lower_bound" : "FAILURES -- see above", sink);
    return ok ? 0 : 1;
}
