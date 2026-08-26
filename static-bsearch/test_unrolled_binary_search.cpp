/*
 *
 * test_unrolled_binary_search.cpp
 *
 * Tests for unrolled_binary_search.hpp.
 *
 * Build with: g++ -std=c++17 -O2 test_unrolled_binary_search.cpp
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

#include <array>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <limits>
#include <unordered_set>
#include "unrolled_binary_search.hpp"
#include <cassert>

template<size_t N>
void init_random(std::array<int32_t, N>& arr, uint32_t seed = 12345)
{
    // Floyd's algorithm, from Bentley's More Programming Pearls, Column 13
    // ("A Sample of Brilliance", with Bob Floyd as guest author). Chooses N
    // distinct values from a universe of U in O(N) time, independent of U --
    // which is what Knuth's Algorithm S cannot do, since it must visit every
    // one of the U candidates.
    //
    // Here U is the whole 2^32 space of int32_t values, so the sample is
    // sparse and the array has gaps. Gaps are the point: without them there
    // is no way to search for a value that is absent but in range.
    constexpr uint64_t U = 1ull << 32;
    static_assert(N < U, "sample larger than universe");

    std::mt19937_64 rng(seed);
    std::unordered_set<uint32_t> chosen;
    chosen.reserve(N * 2);

    for (uint64_t j = U - N; j < U; ++j) {
        std::uniform_int_distribution<uint64_t> pick(0, j);
        uint32_t t = static_cast<uint32_t>(pick(rng));
        if (!chosen.insert(t).second)
            chosen.insert(static_cast<uint32_t>(j));
    }

    size_t i = 0;
    for (uint32_t v : chosen) arr[i++] = static_cast<int32_t>(v);
    std::sort(arr.begin(), arr.end());
}

template<size_t N>
void test_int_array( )
{
    constexpr std::array<int, 16> arr = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    for (int i = 1; i <= 16; ++i) {
        int idx = unrolled_binary_search(arr, i);
        if (idx != i - 1) {
            std::cout << "Test failed for " << i << ": got " << idx << ", expected " << (i-1) << '\n';
        }
    }
    // Test not found
    int idx = unrolled_binary_search(arr, 100);
    if (idx != -1) {
        std::cout << "Test failed for not found: got " << idx << ", expected -1\n";
    }
    std::array<int32_t, N> arrN;
    init_random<N>( arrN );
    for (int i = 0; i < N; ++i) {
        int idx = unrolled_binary_search( arrN, arrN[i] );
        if (idx != i) {
            std::cout << "Test failed for index " << i << "(" << arr[i] << "): got " << idx << ", expected " << (i-1) << '\n';
        }
    }
}

void test_string_array() {
    std::array<std::string, 4> arr = {"apple", "banana", "cherry", "date"};
    int idx = unrolled_binary_search(arr, std::string("cherry"));
    if (idx != 2) {
        std::cout << "Test failed for string: got " << idx << ", expected 2\n";
    }
    idx = unrolled_binary_search(arr, std::string("fig"));
    if (idx != -1) {
        std::cout << "Test failed for string not found: got " << idx << ", expected -1\n";
    }
}

int binary_search_1024(const std::array<int32_t, 1024>& arr, int32_t target) {
    constexpr size_t k = 1024; // Power of 2
    constexpr size_t probe = 0; // No offset needed
    const int32_t* base = &arr[0];
    size_t idx = 0;
    if (idx + 512 < k && base[idx + 512] <= target) idx += 512;
    if (idx + 256 < k && base[idx + 256] <= target) idx += 256;
    if (idx + 128 < k && base[idx + 128] <= target) idx += 128;
    if (idx + 64 < k && base[idx + 64] <= target) idx += 64;
    if (idx + 32 < k && base[idx + 32] <= target) idx += 32;
    if (idx + 16 < k && base[idx + 16] <= target) idx += 16;
    if (idx + 8 < k && base[idx + 8] <= target) idx += 8;
    if (idx + 4 < k && base[idx + 4] <= target) idx += 4;
    if (idx + 2 < k && base[idx + 2] <= target) idx += 2;
    if (idx + 1 < k && base[idx + 1] <= target) idx += 1;
    if (base[idx] == target) {
        return static_cast<int>(idx);
    }
    return -1;
}

// Test for N=1024 specialized search
void test_binary_search_1024() {
    std::array<int32_t, 1024> arrN;
    init_random<1024>(arrN);
    for (int i = 0; i < 1024; ++i) {
        int idx = binary_search_1024(arrN, arrN[i]);
        if (idx != i) {
            std::cout << "[1024-specialized] Test failed for index " << i << "(" << arrN[i] << "): got " << idx << ", expected " << i << '\n';
        }
    }
}

int binary_search_1000(const std::array<int32_t, 1000>& arr, int32_t target) {
    constexpr size_t k = 512; // Largest power of 2 <= 1000
    constexpr size_t probe = 1000 - k; // 488
    const int32_t* base = nullptr;
    if (target <= arr[probe]) {
        base = &arr[0];
    } else {
        base = &arr[probe];
    }
    size_t idx = 0;
    if (idx + 256 < k && base[idx + 256] <= target) idx += 256;
    if (idx + 128 < k && base[idx + 128] <= target) idx += 128;
    if (idx + 64 < k && base[idx + 64] <= target) idx += 64;
    if (idx + 32 < k && base[idx + 32] <= target) idx += 32;
    if (idx + 16 < k && base[idx + 16] <= target) idx += 16;
    if (idx + 8 < k && base[idx + 8] <= target) idx += 8;
    if (idx + 4 < k && base[idx + 4] <= target) idx += 4;
    if (idx + 2 < k && base[idx + 2] <= target) idx += 2;
    if (idx + 1 < k && base[idx + 1] <= target) idx += 1;
    if (base[idx] == target) {
        // Return the index in the original array
        return static_cast<int>((base - &arr[0]) + idx);
    }
    return -1;
}

// Test for N=1000 specialized search
void test_binary_search_1000() {
    std::array<int32_t, 1000> arrN;
    init_random<1000>(arrN);
    for (int i = 0; i < 1000; ++i) {
        int idx = binary_search_1000(arrN, arrN[i]);
        if (idx != i) {
            std::cout << "[1000-specialized] Test failed for index " << i << "(" << arrN[i] << "): got " << idx << ", expected " << i << '\n';
        }
    }
}

// Force a standalone instantiation so this can be disassembled and compared
// against the hand-written binary_search_1024 above. Without it the template
// is never instantiated and emits no code at all.
template int unrolled_binary_search_1024<int32_t>(const std::array<int32_t, 1024>&,
                                                 const int32_t&);

void test_unrolled_1024() {
    std::array<int32_t, 1024> arrN;
    init_random<1024>(arrN);
    for (int i = 0; i < 1024; ++i) {
        int idx = unrolled_binary_search_1024(arrN, arrN[i]);
        if (idx != i) {
            std::cout << "[unrolled-1024] failed at index " << i
                      << ": got " << idx << '\n';
        }
    }
}

// Search for values that are absent but in range. Impossible to construct
// while init_random emitted consecutive integers.
template<size_t N>
void test_absent_values() {
    std::array<int32_t, N> arrN;
    init_random<N>(arrN);
    long long checked = 0;
    for (size_t i = 0; i + 1 < N; ++i) {
        if (arrN[i+1] - arrN[i] < 2) continue;
        int32_t between = arrN[i] + 1;
        int idx = unrolled_binary_search(arrN, between);
        ++checked;
        if (idx != -1) {
            std::cout << "[absent] " << between << " reported at " << idx << '\n';
        }
    }
    std::cout << "  interior misses checked for N=" << N << ": " << checked << '\n';
}

int main() {
    constexpr size_t N = 1048576;
    constexpr auto steps = make_probe_steps<N>( );
    test_int_array<N>();
    test_int_array<1000>();
    test_binary_search_1000();
    test_binary_search_1024();
    test_string_array();
    test_unrolled_1024();
    test_absent_values<1024>();
    test_absent_values<1000>();
    std::cout << "All tests completed.\n";
    return 0;
}
