cmake_minimum_required(VERSION 3.25)

# To add more shaders simply copy the vulkan_compile_shader command and replace it with your new shader
vulkan_compile_shader(INFILE matxVec.comp
    OUTFILE MatxVec.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE gaussSeidel.comp
    OUTFILE GaussSeidel.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE conjugateGradient.comp
    OUTFILE CG.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

# Then add it to the library, so you can access it later in your code
add_library(shaders INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/MatxVec.hpp" "${CMAKE_CURRENT_BINARY_DIR}/GaussSeidel.hpp" "${CMAKE_CURRENT_BINARY_DIR}/CG.hpp")

target_include_directories(shaders INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
