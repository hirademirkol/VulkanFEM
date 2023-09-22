cmake_minimum_required(VERSION 3.22)
project(Master)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(external/kompute ${CMAKE_CURRENT_BINARY_DIR}/kompute_build)
list(APPEND CMAKE_PREFIX_PATH "${kompute_SOURCE_DIR}/cmake")
include_directories(external/kompute/src/include)

include_directories(external/eigen)

include_directories(${CMAKE_SOURCE_DIR}/include)
add_subdirectory(src)
add_subdirectory(shaders)