cmake_minimum_required(VERSION 3.13)
project(truesum LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(TRUESUM_ENABLE_AVX512 "Compile the AVX-512 kernels when the compiler supports them (selected at run time)" ON)
option(TRUESUM_ENABLE_CUDA "Build the CUDA accumulation matrix when a CUDA toolkit is found" ON)
option(TRUESUM_BUILD_BENCH "Build the benchmark drivers that replicate the paper's tables" ON)

add_library(truesum
  src/limbs.cpp
  src/limb_column.cpp
  src/kernels_scalar.cpp
  src/kernels_dispatch.cpp
  src/accumulation_matrix.cpp
  src/survey.cpp
)

# The AVX-512 kernels get their own translation unit with their own flags, so
# the rest of the library stays baseline and the binary still runs anywhere.
# Nothing in that file may run before the CPU check in kernels_dispatch.cpp.
include(CheckCXXCompilerFlag)
set(TRUESUM_AVX512_FLAGS -mavx512f -mavx512dq -mavx512vl -mavx512bw
                      -mavx512cd -mavx512vpopcntdq)
check_cxx_compiler_flag("${TRUESUM_AVX512_FLAGS}" TRUESUM_COMPILER_HAS_AVX512)
if(TRUESUM_ENABLE_AVX512 AND TRUESUM_COMPILER_HAS_AVX512)
  target_sources(truesum PRIVATE src/kernels_avx512.cpp)
  set_source_files_properties(src/kernels_avx512.cpp PROPERTIES
    COMPILE_OPTIONS "${TRUESUM_AVX512_FLAGS}")
  target_compile_definitions(truesum PRIVATE TRUESUM_HAVE_AVX512)
endif()
target_include_directories(truesum PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
# The accumulator can spread its columns over worker threads. Standard C++17,
# but the linker still wants pthreads; PUBLIC so consumers inherit it.
find_package(Threads REQUIRED)
target_link_libraries(truesum PUBLIC Threads::Threads)
target_compile_options(truesum PRIVATE -Wall -Wextra)

enable_testing()

# CUDA is optional and lives in its own target: the device path is a separate
# container with device-resident storage, not a third kernel behind the CPU
# dispatch, so nothing in libtruesum depends on it.
include(CheckLanguage)
if(TRUESUM_ENABLE_CUDA)
  check_language(CUDA)
endif()
if(TRUESUM_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  add_library(truesum_cuda src/cuda_accumulation_matrix.cu)
  set_target_properties(truesum_cuda PROPERTIES CUDA_ARCHITECTURES "70;80;86")
  target_link_libraries(truesum_cuda PUBLIC truesum)

  add_executable(test_cuda tests/test_cuda.cu)
  set_target_properties(test_cuda PROPERTIES CUDA_ARCHITECTURES "70;80;86")
  target_link_libraries(test_cuda PRIVATE truesum_cuda)
  add_test(NAME truesum_cuda COMMAND test_cuda)

  if(TRUESUM_BUILD_BENCH)
    foreach(b bench_table_gpu bench_adjust_gpu)
      add_executable(${b} bench/${b}.cu)
      set_target_properties(${b} PROPERTIES CUDA_ARCHITECTURES "70;80;86")
      target_include_directories(${b} PRIVATE bench)
      target_link_libraries(${b} PRIVATE truesum_cuda)
    endforeach()
  endif()
endif()

add_executable(test_truesum tests/test_truesum.cpp)
target_link_libraries(test_truesum PRIVATE truesum)
target_compile_options(test_truesum PRIVATE -Wall -Wextra)
add_test(NAME truesum COMMAND test_truesum)

add_executable(demo examples/demo.cpp)
target_link_libraries(demo PRIVATE truesum)

# Benchmark drivers. Each replicates one table of the paper; bench/replicate.sh
# runs them all and writes their output beside bench/expected/.
if(TRUESUM_BUILD_BENCH)
  foreach(b bench_table_cpu bench_adjust_cpu bench_kernels)
    add_executable(${b} bench/${b}.cpp)
    target_include_directories(${b} PRIVATE bench)
    target_link_libraries(${b} PRIVATE truesum)
    target_compile_options(${b} PRIVATE -Wall -Wextra)
  endforeach()
endif()
