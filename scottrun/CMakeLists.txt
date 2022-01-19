project(smooth)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)
add_custom_target(extern_b3d2 COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/hp/alice/software/b3d2)
add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/smooth/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3 -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-result -std=c++11 -mfpmath=sse")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/hp/alice/software/b3d2/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include ${GITHOME_MSU}/smooth/software/include)        

link_directories(${GITHOME_MSU}/smooth/software/lib ${GITHOME_MSU}/hp/alice/software/b3d2/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/coral/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib commonutils)
find_library(${GITHOME_MSU}/coral/software/lib coral)
find_library(${GITHOME_MSU}/hp/alice/software/b3d2/lib b3d_alice)
find_library(${GITHOME_MSU}/smooth/software/lib smooth)

add_dependencies(extern_smooth extern_b3d2 extern_coral extern_commonutils)

#################################################
set(SOURCE_FILES ./smoothmain.cc)
add_executable(smoothy ${CMAKE_CURRENT_SOURCE_DIR}/smoothmain.cc)
add_dependencies(smoothy extern_smooth extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(smoothy ${GITHOME_MSU}/smooth/software/lib/libsmooth.a ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
