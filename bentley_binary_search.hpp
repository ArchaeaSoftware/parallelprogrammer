// Unrolled search for power-of-two segment
template <typename T, std::size_t... Steps>
int bentley_binary_search_unrolled(const T* arr, std::size_t len, const T& target, std::integer_sequence<std::size_t, Steps...>) {
    std::size_t idx = 0;
    ((idx += (idx + Steps < len && arr[idx + Steps] <= target ? Steps : 0)), ...);
    if (arr[idx] == target) return static_cast<int>(idx);
    return -1;
}
// Modern C++17+ approach: constexpr function to generate probe steps and lambda for unrolling


#include <utility>
#include <type_traits>


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
int bentley_binary_search(const std::array<T, N>& arr, const T& target) {
    constexpr std::size_t k = floor_power_of_two(N);
    if constexpr (N == k) {
        // Power of 2: unrolled search
        return bentley_binary_search_unrolled(&arr[0], N, target, make_probe_steps<N>());
    } else {
        constexpr std::size_t probe = N - k;
        if (target <= arr[probe]) {
            // Unrolled search in [0, k-1]
            int res = bentley_binary_search_unrolled(&arr[0], k, target, make_probe_steps<k>());
            return (res == -1) ? -1 : res;
        } else {
            // Unrolled search in [N-k, N-1]
            int res = bentley_binary_search_unrolled(&arr[N - k], k, target, make_probe_steps<k>());
            return (res == -1) ? -1 : static_cast<int>(N - k + res);
        }
    }
}


// Example usage:
// constexpr std::array<int, 16> v = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
// int idx = bentley_binary_search<int, 16, 8, 4, 2, 1>(v, 7); // idx == 6


// Unrolled binary search for 1024 elements, for codegen comparison
template <typename T>
int bentley_binary_search_1024(const std::array<T, 1024>& arr, const T& target) {
    constexpr std::size_t N = 1024;
    constexpr auto steps = make_probe_steps<N>();
    return bentley_binary_search_unrolled(&arr[0], N, target, steps);
}