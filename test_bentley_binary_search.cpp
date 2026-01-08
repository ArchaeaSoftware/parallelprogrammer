
#include <array>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include "bentley_binary_search.hpp"

template<size_t N>
void init_random( std::array<int32_t, N>& arr )
{
    for (std::size_t i = 0; i < arr.size(); ++i) {
        arr[i] = static_cast<int>( rand() ); // Simple deterministic initialization
    }
    std::sort(arr.begin(), arr.end());
}

void test_int_array() {
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

int main() {
    test_int_array();
    test_string_array();
    std::cout << "All tests completed.\n";
    return 0;
}
