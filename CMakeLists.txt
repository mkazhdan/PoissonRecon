## -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)
project(AdaptiveMultigridSolvers)


set(CMAKE_VERBOSE_MAKEFILE ON)

## -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++14 -pthread -Wno-invalid-offsetof -fPIC -Wno-deprecated")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -Ofast -DRELEASE -funroll-loops -ffast-math")


## -----------------------------------------------------------------------------
find_package(OpenMP)


## -----------------------------------------------------------------------------
# all include paths
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/Src)
include_directories(${PROJECT_SOURCE_DIR}/PNG)
include_directories(${PROJECT_SOURCE_DIR}/JPEG)
include_directories(${PROJECT_SOURCE_DIR}/ZLIB)

## -----------------------------------------------------------------------------
## build dependencies as libraries
## -----------------------------------------------------------------------------
set(TRG_PNG  pr_png)
set(TRG_JPEG pr_jpeg)
set(TRG_ZLIB pr_zlib)

file(GLOB_RECURSE SRCS_PNG "PNG/*.c")
file(GLOB_RECURSE SRCS_JPEG "JPEG/*.cpp")
file(GLOB_RECURSE SRCS_ZLIB "ZLIB/*.c")

add_library(${TRG_PNG}  STATIC ${SRCS_PNG})
add_library(${TRG_JPEG} SHARED ${SRCS_JPEG})
add_library(${TRG_ZLIB} SHARED ${SRCS_ZLIB})


## -----------------------------------------------------------------------------
## build poisson reconstruction
## -----------------------------------------------------------------------------
add_executable(poisson_recon
               ${PROJECT_SOURCE_DIR}/Src/PoissonRecon.cpp)

target_link_libraries(poisson_recon
                      ${TRG_PNG} ${TRG_JPEG} ${TRG_ZLIB})

## -----------------------------------------------------------------------------
