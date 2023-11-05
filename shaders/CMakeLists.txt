cmake_minimum_required(VERSION 3.24)

# To add more shaders simply copy the vulkan_compile_shader command and replace it with your new shader
vulkan_compile_shader(INFILE MatxVec.comp
    OUTFILE MatxVec.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE FixedNodes.comp
    OUTFILE FixedNodes.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE VecDotVec.comp
    OUTFILE VecDotVec.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE ConjugateGradient_1.comp
    OUTFILE ConjugateGradient_1.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE ConjugateGradient_2.comp
    OUTFILE ConjugateGradient_2.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")
    
vulkan_compile_shader(INFILE CWiseMult.comp
    OUTFILE CWiseMult.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE PreSmoothAndAXPY.comp
    OUTFILE PreSmoothAndAXPY.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE Smooth.comp
    OUTFILE Smooth.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE Restrict.comp
    OUTFILE Restrict.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

vulkan_compile_shader(INFILE Interpolate.comp
    OUTFILE Interpolate.hpp
    NAMESPACE "shaders"
    RELATIVE_PATH "${kompute_SOURCE_DIR}/cmake")

# Then add it to the library, so you can access it later in your code
add_library(shaders INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/MatxVec.hpp" 
                              "${CMAKE_CURRENT_BINARY_DIR}/FixedNodes.hpp" 
                              "${CMAKE_CURRENT_BINARY_DIR}/VecDotVec.hpp" 
                              "${CMAKE_CURRENT_BINARY_DIR}/ConjugateGradient_1.hpp" 
                              "${CMAKE_CURRENT_BINARY_DIR}/ConjugateGradient_2.hpp"
                              "${CMAKE_CURRENT_BINARY_DIR}/CWiseMult.hpp" 
                              "${CMAKE_CURRENT_BINARY_DIR}/PreSmoothAndAXPY.hpp" 
                              "${CMAKE_CURRENT_BINARY_DIR}/Smooth.hpp" 
                              "${CMAKE_CURRENT_BINARY_DIR}/Restrict.hpp" 
                              "${CMAKE_CURRENT_BINARY_DIR}/Interpolate.hpp")

target_include_directories(shaders INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
