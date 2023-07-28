cmake_minimum_required(VERSION 3.25)

add_executable(${CMAKE_PROJECT_NAME} main.cpp)

target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${EIGEN_STANDARD_LIBRARIES_TO_LINK_TO})
# target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC csparse_static)

# target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE shaders kompute::kompute)

add_library(utils SHARED utils.cpp)
add_library(inout SHARED inout.cpp)
add_library(ke SHARED Ke.cpp)
add_library(fem SHARED FEM.cpp)
target_link_libraries(fem PUBLIC utils)

target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC utils inout ke fem)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

add_custom_command(TARGET ${CMAKE_PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${PROJECT_SOURCE_DIR}/data"
        "$<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>/data/")