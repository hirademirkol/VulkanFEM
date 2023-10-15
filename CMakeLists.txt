cmake_minimum_required(VERSION 3.22)
project(Master)

set(CMAKE_CXX_STANDARD 20)

message("Applying patches for kompute")
file(GLOB PATCHES "${CMAKE_SOURCE_DIR}/external/patches/kompute/*.patch")
foreach(patch_file ${PATCHES})
    execute_process(COMMAND ${GIT_EXECUTABLE} apply --3way --whitespace=nowarn ${patch_file}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external/kompute)
endforeach(patch_file ${PATCHES})

set(KOMPUTE_OPT_LOG_LEVEL "Warn")
add_subdirectory(external/kompute ${CMAKE_CURRENT_BINARY_DIR}/kompute_build)
list(APPEND CMAKE_PREFIX_PATH "${kompute_SOURCE_DIR}/cmake")
include_directories(external/kompute/src/include)

include_directories(external/eigen)

include_directories(${CMAKE_SOURCE_DIR}/include)
add_subdirectory(src)
add_subdirectory(shaders)