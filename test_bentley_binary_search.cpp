#include <array>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <limits>
#include "bentley_binary_search.hpp"
#include <cassert>

template<size_t N>
void init_random(std::array<int32_t, N>& arr, uint32_t seed = 12345)
{
    // Distinct values in ascending order, drawn from the whole int32_t range.
    // The range matters: the values must be sparse in it, so that the array
    // has gaps. Without gaps there is no way to search for a value that is
    // absent but lies between two elements, which is the case most likely to
    // expose an off-by-one in the probe sequence.
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int32_t> dist(
        std::numeric_limits<int32_t>::min(),
        std::numeric_limits<int32_t>::max());

    for (auto& v : arr) v = dist(rng);
    std::sort(arr.begin(), arr.end());
    for (;;) {
        auto last = std::unique(arr.begin(), arr.end());
        if (last == arr.end()) break;
        for (auto it = last; it != arr.end(); ++it) *it = dist(rng);
        std::sort(arr.begin(), arr.end());
    }
}

template<size_t N>
void test_int_array( )
{
    constexpr std::array<int, 16> arr = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    for (int i = 1; i <= 16; ++i) {
        int idx = bentley_binary_search(arr, i);
        if (idx != i - 1) {
            std::cout << "Test failed for " << i << ": got " << idx << ", expected " << (i-1) << '\n';
        }
    }
    // Test not found
    int idx = bentley_binary_search(arr, 100);
    if (idx != -1) {
        std::cout << "Test failed for not found: got " << idx << ", expected -1\n";
    }
    std::array<int32_t, N> arrN;
    init_random<N>( arrN );
    for (int i = 0; i < N; ++i) {
        int idx = bentley_binary_search( arrN, arrN[i] );
        if (idx != i) {
            std::cout << "Test failed for index " << i << "(" << arr[i] << "): got " << idx << ", expected " << (i-1) << '\n';
        }
    }
}

void test_string_array() {
    std::array<std::string, 4> arr = {"apple", "banana", "cherry", "date"};
    int idx = bentley_binary_search(arr, std::string("cherry"));
    if (idx != 2) {
        std::cout << "Test failed for string: got " << idx << ", expected 2\n";
    }
    idx = bentley_binary_search(arr, std::string("fig"));
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
template int bentley_binary_search_1024<int32_t>(const std::array<int32_t, 1024>&,
                                                 const int32_t&);

void test_bentley_1024() {
    std::array<int32_t, 1024> arrN;
    init_random<1024>(arrN);
    for (int i = 0; i < 1024; ++i) {
        int idx = bentley_binary_search_1024(arrN, arrN[i]);
        if (idx != i) {
            std::cout << "[bentley-1024] failed at index " << i
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
        int idx = bentley_binary_search(arrN, between);
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
    test_bentley_1024();
    test_absent_values<1024>();
    test_absent_values<1000>();
    std::cout << "All tests completed.\n";
    return 0;
}
