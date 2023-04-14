cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(smooth)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)
add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/smooth/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -std=c++17")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include ${GITHOME_MSU}/smooth/software/include)        

link_directories(${GITHOME_MSU}/smooth/software/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/coral/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/coral/software/lib msu_coral)
find_library(${GITHOME_MSU}/smooth/software/lib smooth)

#################################################
set(SOURCE_FILES ./smoothmain.cc)
add_executable(smoothy ${CMAKE_CURRENT_SOURCE_DIR}/smoothmain.cc)
add_dependencies(smoothy extern_smooth extern_coral extern_commonutils)
target_link_libraries(smoothy ${GITHOME_MSU}/smooth/software/lib/libsmooth.a ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################
set(SOURCE_FILES ./keepscore.cc)
add_executable(keepscore ${CMAKE_CURRENT_SOURCE_DIR}/keepscore.cc)
add_dependencies(keepscore extern_smooth extern_coral extern_commonutils)
target_link_libraries(keepscore ${GITHOME_MSU}/smooth/software/lib/libsmooth.a ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################
set(SOURCE_FILES ./coultest.cc)
add_executable(coultest ${CMAKE_CURRENT_SOURCE_DIR}/coultest.cc)
add_dependencies(coultest extern_smooth extern_coral extern_commonutils)
target_link_libraries(coultest ${GITHOME_MSU}/smooth/software/lib/libsmooth.a ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
