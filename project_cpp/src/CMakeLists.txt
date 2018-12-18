#Minimum version of CMake
cmake_minimum_required(VERSION 2.8)

include_directories(library)

add_library(SteerablePyramid SHARED SteerablePyramid.cpp)

add_executable(main main.cpp)
# lopencv_imgproc, pour la version 1 de OpenCV
# à changer en lopencv_imgcodecs si utilisation de la v2
target_link_libraries(main SteerablePyramid -lopencv_imgproc ${OpenCV_LIBS})
