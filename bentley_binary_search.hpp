// Modern C++17+ approach: constexpr function to generate probe steps and lambda for unrolling


#include <utility>
#include <type_traits>

// Helper to generate descending powers of 2 for probe steps
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

template <typename T, std::size_t N, std::size_t... Steps>
int bentley_binary_search_impl(const std::array<T, N>& arr, const T& target, std::integer_sequence<std::size_t, Steps...>) {
    std::size_t idx = 0;
    ((idx += (idx + Steps < N && arr[idx + Steps] <= target ? Steps : 0)), ...);
    if (arr[idx] == target) return static_cast<int>(idx);
    // For non-power-of-two sizes, check the last element
    if constexpr (N > 0) {
        if (idx != N - 1 && arr[N - 1] == target) return static_cast<int>(N - 1);
    }
    return -1;
}

template <typename T, std::size_t N>
int bentley_binary_search(const std::array<T, N>& arr, const T& target) {
    constexpr auto steps = make_probe_steps<N>();
    return bentley_binary_search_impl(arr, target, steps);
}
// Jon Bentley's optimized binary search for constant-size arrays
// Loop unrolled, probes with descending powers of 2
// Templatized for any comparable type, array size, and probe steps
#include <array>
#include <cstddef>

// Returns index of target in sorted array, or -1 if not found, with explicit probe steps
template <typename T, std::size_t N, std::size_t... Steps>
int bentley_binary_search_unrolled(const std::array<T, N>& arr, const T& target) {
    std::size_t idx = 0;
    // Unroll the loop using the probe steps
    ((idx += (idx + Steps < N && arr[idx + Steps] <= target ? Steps : 0)), ...);
    if (arr[idx] == target) return static_cast<int>(idx);
    return -1;
}

// Example usage:
// constexpr std::array<int, 16> v = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
// int idx = bentley_binary_search<int, 16, 8, 4, 2, 1>(v, 7); // idx == 6
