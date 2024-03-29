cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})
project(msu_smooth CXX C)
message("CXX is ${CXX}")


set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)

message("eigen3_include_dir is ${EIGEN3_INCLUDE_DIR}")

include_directories(SYSTEM include ${GITHOME_MSU}/commonutils/software/include ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})

file(GLOB smooth_src src/*.cc)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(smooth STATIC ${smooth_src})
