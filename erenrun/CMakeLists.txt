cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_VERBOSE_MAKEFILE on)
project(smooth)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/smooth/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -std=c++17")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/smooth/software/include)      

link_directories(${GITHOME_MSU}/smooth/software/lib ${GITHOME_MSU}/commonutils/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/smooth/software/lib smooth)

#################################################
set(SOURCE_FILES ./simplextest.cc)
add_executable(simplextest ${CMAKE_CURRENT_SOURCE_DIR}/simplextest.cc)
add_dependencies(simplextest extern_smooth extern_commonutils)
target_link_libraries(simplextest ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################
set(SOURCE_FILES ./fakemodel.cc)
add_executable(fakemodel ${CMAKE_CURRENT_SOURCE_DIR}/fakemodel.cc)
add_dependencies(fakemodel extern_smooth extern_commonutils)
target_link_libraries(fakemodel ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################
set(SOURCE_FILES ./smoothy.cc)
add_executable(smoothy ${CMAKE_CURRENT_SOURCE_DIR}/smoothy.cc)
add_dependencies(smoothy extern_smooth extern_commonutils)
target_link_libraries(smoothy ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################
set(SOURCE_FILES ./smoothy_writecoefficients.cc)
add_executable(smoothy_writecoefficients ${CMAKE_CURRENT_SOURCE_DIR}/smoothy_writecoefficients.cc)
add_dependencies(smoothy_writecoefficients extern_smooth extern_commonutils)
target_link_libraries(smoothy_writecoefficients ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################
set(SOURCE_FILES ./smoothy_readcoefficients.cc)
add_executable(smoothy_readcoefficients ${CMAKE_CURRENT_SOURCE_DIR}/smoothy_readcoefficients.cc)
add_dependencies(smoothy_readcoefficients extern_smooth extern_commonutils)
target_link_libraries(smoothy_readcoefficients ${GITHOME_MSU}/smooth/software/lib/libsmooth.a  ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
