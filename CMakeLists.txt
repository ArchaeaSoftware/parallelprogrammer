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
  src/column_accumulator.cpp
)
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
