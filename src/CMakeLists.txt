cmake_minimum_required(VERSION 3.20)

add_executable(${CMAKE_PROJECT_NAME} main.cpp)

target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE shaders kompute::kompute)
