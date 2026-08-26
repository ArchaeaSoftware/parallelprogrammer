/*
 *
 * unrolled_binary_search.hpp
 *
 * Unrolled binary search for fixed-size sorted arrays.
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
#include <cstddef>
#include <type_traits>
#include <utility>

// Scalars are cheaper passed by value: taking a reference obliges the caller to
// materialise the argument in memory so that its address can be passed. Class
// types still want const&.
template <typename T>
using search_arg_t = std::conditional_t<std::is_scalar_v<T>, T, const T&>;

// Unrolled search for power-of-two segment
template <std::size_t Len, typename T, std::size_t... Steps>
std::ptrdiff_t unrolled_search_segment(const T* arr, search_arg_t<T> target, std::integer_sequence<std::size_t, Steps...>) {
    std::size_t idx = 0;
    ((idx += (idx + Steps < Len && arr[idx + Steps] <= target ? Steps : 0)), ...);
    if (arr[idx] == target) return static_cast<std::ptrdiff_t>(idx);
    return -1;
}

// Helper: largest power of 2 <= N
constexpr std::size_t floor_power_of_two(std::size_t n) {
    std::size_t p = 1;
    while (p * 2 <= n) p *= 2;
    return p;
}

template <std::size_t N, std::size_t... Steps>
constexpr auto make_probe_steps_impl(std::integer_sequence<std::size_t, Steps...>) {
    constexpr std::size_t next = (N / 2) >> sizeof...(Steps);
    if constexpr (next > 0) {
        return make_probe_steps_impl<N>(std::integer_sequence<std::size_t, Steps..., next>{});
    } else {
        return std::integer_sequence<std::size_t, Steps...>{};
    }
}

template <std::size_t N>
constexpr auto make_probe_steps() {
    return make_probe_steps_impl<N>(std::integer_sequence<std::size_t>{});
}



template <typename T, std::size_t N>
std::ptrdiff_t unrolled_binary_search(const std::array<T, N>& arr, const T& target) {
    constexpr std::size_t k = floor_power_of_two(N);
    if constexpr (N == k) {
        // Power of 2: unrolled search
        return unrolled_search_segment<N>(&arr[0], target, make_probe_steps<N>());
    } else {
        constexpr std::size_t probe = N - k;
        if (target <= arr[probe]) {
            // Unrolled search in [0, k-1]
            int res = unrolled_search_segment<k>(&arr[0], target, make_probe_steps<k>());
            return (res == -1) ? -1 : res;
        } else {
            // Unrolled search in [N-k, N-1]
            int res = unrolled_search_segment<k>(&arr[N - k], target, make_probe_steps<k>());
            return (res == -1) ? -1 : static_cast<std::ptrdiff_t>(N - k) + res;
        }
    }
}


// Example usage:
// constexpr std::array<int, 16> v = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
// std::ptrdiff_t idx = unrolled_binary_search(v, 7);   // idx == 6


// Unrolled binary search for 1024 elements, for codegen comparison
template <typename T>
std::ptrdiff_t unrolled_binary_search_1024(const std::array<T, 1024>& arr, const T& target) {
    constexpr std::size_t N = 1024;
    constexpr auto steps = make_probe_steps<N>();
    return unrolled_search_segment<N>(&arr[0], target, steps);
}