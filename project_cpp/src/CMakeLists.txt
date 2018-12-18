#Minimum version of CMake
cmake_minimum_required(VERSION 2.8)

include_directories(library)

add_library(SteerablePyramid SHARED SteerablePyramid.cpp)
add_executable(main main.cpp)
target_link_libraries(main SteerablePyramid ${OpenCV_LIBS})
