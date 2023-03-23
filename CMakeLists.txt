cmake_minimum_required(VERSION 3.22)
project(Master)

set(CMAKE_CXX_STANDARD 14)

message("Applying patches for kompute")

file(GLOB PATCHES "${CMAKE_SOURCE_DIR}/patches/*.patch")
foreach(patch_file ${PATCHES})
    execute_process(COMMAND ${GIT_EXECUTABLE} apply --3way --whitespace=nowarn ${patch_file}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external/kompute)
endforeach(patch_file ${PATCHES})

add_subdirectory(external/kompute ${CMAKE_CURRENT_BINARY_DIR}/kompute_build)
include_directories(external/kompute/src/include)

list(APPEND CMAKE_PREFIX_PATH "${kompute_SOURCE_DIR}/cmake")

add_subdirectory(src)
add_subdirectory(shaders)
