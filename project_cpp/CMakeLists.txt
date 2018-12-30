#Minimum version of CMake
cmake_minimum_required(VERSION 2.8)

project(Wave_project)

include(CTest)

FIND_PACKAGE(OpenMP REQUIRED)
message(STATUS $ENV{OPENMP_FOUND})
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package( OpenCV REQUIRED )
message(STATUS $ENV{OPENCV_FOUND})

if(NOT "${OPENCV_VERSION}" STRGREATER "2")
    set(OpenLibraries "-lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs")
else()
    set(OpenLibraries "-lopencv_imgproc")
endif()

#SOURCES
add_subdirectory(src)

#TESTS
add_subdirectory(test)
