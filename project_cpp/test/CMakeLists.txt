#Minimum version of CMake
cmake_minimum_required(VERSION 2.8)

include_directories(${PROJECT_SOURCE_DIR}/src)

#EXECUTABLES

add_executable(test_compilation test_compilation.cpp)
add_executable(test_filters test_filters.cpp)

# LINK LIBRARIES

target_link_libraries(test_compilation SteerablePyramid ${OpenLibraries})
target_link_libraries(test_filters SteerablePyramid ${OpenLibraries})

# Transform into tests

add_test(TestCompilation test_compilation)
