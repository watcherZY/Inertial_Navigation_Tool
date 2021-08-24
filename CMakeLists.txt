cmake_minimum_required( VERSION 2.8 )
project(SINS)

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-O3" )

# 添加Eigen头文件
include_directories( "/usr/include/eigen3" )

aux_source_directory(. DIR_SRCS)
add_executable( SINS ${DIR_SRCS})
