cmake_minimum_required(VERSION 3.22)
project(Master)

set(CMAKE_CXX_STANDARD 20)

message("Applying patches for kompute")

file(GLOB PATCHES "${CMAKE_SOURCE_DIR}/patches/*.patch")
foreach(patch_file ${PATCHES})
    execute_process(COMMAND ${GIT_EXECUTABLE} apply --3way --whitespace=nowarn ${patch_file}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external/kompute)
endforeach(patch_file ${PATCHES})

add_subdirectory(external/kompute ${CMAKE_CURRENT_BINARY_DIR}/kompute_build)
list(APPEND CMAKE_PREFIX_PATH "${kompute_SOURCE_DIR}/cmake")
include_directories(external/kompute/src/include ${CMAKE_SOURCE_DIR}/include)

add_subdirectory(external/SuiteSparse/CSparse)
include_directories(external/SuiteSparse/CSparse/Include)
# set( CXSparse_ROOT "${CMAKE_SOURCE_DIR}/external/SuiteSparse/CXSparse")

# set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
#     ${CMAKE_SOURCE_DIR}/external/SuiteSparse/CXSparse/cmake_modules
#     ${CMAKE_SOURCE_DIR}/external/SuiteSparse/SuiteSparse_config/cmake_modules)

# find_package ( CXSparse 4.0.0 REQUIRED )
# include_directories ( ${CXSPARSE_INCLUDE_DIR} )

# list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/SuiteSparse/SuiteSparse_config/cmake_modules")
# add_subdirectory(external/SuiteSparse/SuiteSparse_config)
# set(SUITESPARSE_CONFIG_LIBRARY "${CMAKE_SOURCE_DIR}/external/SuiteSparse/SuiteSparse_config/build")
# add_subdirectory(external/SuiteSparse/CXSparse)

add_subdirectory(src)
add_subdirectory(shaders)
