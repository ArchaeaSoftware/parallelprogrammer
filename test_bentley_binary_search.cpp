#include <array>
#include <iostream>
#include <cstdlib>
#include "bentley_binary_search.hpp"
#include <cassert>

template<size_t N>
void init_random(std::array<int32_t, N>& arr)
{
    // Knuth's algorithm: select N unique numbers from 1..N in sorted order
    // For this use case, M = N
    int im = 0;
    for (int in = 0; in < static_cast<int>(N) && im < static_cast<int>(N); ++in) {
        int rn = static_cast<int>(N) - in;
        int rm = static_cast<int>(N) - im;
        if (rand() % rn < rm) {
            arr[im++] = in;
        }
    }
    // Ensure all slots filled
    assert(im == static_cast<int>(N));
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

int main() {
    constexpr size_t N = 1048576;
    constexpr auto steps = make_probe_steps<N>( );
    test_int_array<N>();
    test_int_array<1000>();
    test_binary_search_1000();
    test_binary_search_1024();
    test_string_array();
    std::cout << "All tests completed.\n";
    return 0;
}
