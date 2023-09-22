cmake_minimum_required(VERSION 3.24)

add_executable(${CMAKE_PROJECT_NAME} main.cpp FEM.cpp Ke.cpp inout.cpp utils.cpp)

target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${EIGEN_STANDARD_LIBRARIES_TO_LINK_TO})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE shaders kompute::kompute)

add_custom_command(TARGET ${CMAKE_PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${PROJECT_SOURCE_DIR}/data"
        "$<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>/data/")