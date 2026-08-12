cmake_minimum_required(VERSION 3.13)
project(column_block_fp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(cbfp
  src/limbs.cpp
  src/limb_column.cpp
  src/kernels_scalar.cpp
  src/kernels_dispatch.cpp
  src/column_accumulator.cpp
)

# The AVX-512 kernels get their own translation unit with their own flags, so
# the rest of the library stays baseline and the binary still runs anywhere.
# Nothing in that file may run before the CPU check in kernels_dispatch.cpp.
include(CheckCXXCompilerFlag)
set(CBFP_AVX512_FLAGS -mavx512f -mavx512dq -mavx512vl -mavx512bw)
check_cxx_compiler_flag("${CBFP_AVX512_FLAGS}" CBFP_COMPILER_HAS_AVX512)
if(CBFP_COMPILER_HAS_AVX512)
  target_sources(cbfp PRIVATE src/kernels_avx512.cpp)
  set_source_files_properties(src/kernels_avx512.cpp PROPERTIES
    COMPILE_OPTIONS "${CBFP_AVX512_FLAGS}")
  target_compile_definitions(cbfp PRIVATE CBFP_HAVE_AVX512)
endif()
target_include_directories(cbfp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_options(cbfp PRIVATE -Wall -Wextra)

enable_testing()

add_executable(test_cbfp tests/test_cbfp.cpp)
target_link_libraries(test_cbfp PRIVATE cbfp)
target_compile_options(test_cbfp PRIVATE -Wall -Wextra)
add_test(NAME cbfp COMMAND test_cbfp)

add_executable(demo examples/demo.cpp)
target_link_libraries(demo PRIVATE cbfp)
