cmake_minimum_required(VERSION 3.9)
project(PatternRefinement)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)

include_directories(
        ${OPENCV_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIRS}
        ${Ceres_INCLUDE_DIRS}
)

add_executable(PatternRefinement
        main.cpp
        intensity_functor.cpp
        intensity_functor.h)

target_link_libraries(PatternRefinement
        ceres
        ${OpenCV_LIBS})