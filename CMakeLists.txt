## -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)
project(AdaptiveMultigridSolvers)

set(CMAKE_VERBOSE_MAKEFILE ON)

## -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set custom compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CFLAGS "-std=c++14 -pthread -Wno-deprecated -Wno-invalid-offsetof")
    set(LFLAGS "-lstdc++ -lpthread")
else()
    set(CFLAGS, "-std=c++14 -pthread -Wno-deprecated -Wno-invalid-offsetof -Wno-dangling-else")
    set(LFLAGS, "-lstdc++")
endif()


set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} ${CFLAGS} -g3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CFLAGS} -Ofast -DRELEASE -funroll-loops -ffast-math")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG   "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${LFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${LFLAGS} -Ofast")

## -----------------------------------------------------------------------------
find_package(OpenMP REQUIRED)

## -----------------------------------------------------------------------------
find_package(PNG)        # PNG depends on zlib and adds zlib to path already
if(PNG_FOUND)
    #message("-- png lib found: include_dir=(${PNG_INCLUDE_DIRS}), libs=(${PNG_LIBRARIES})")
else()

    set(PNG_PATH "${PROJECT_SOURCE_DIR}/PNG")
    message("-- png lib not found: will build from (${PNG_PATH})")

    set(PNG_INCLUDE_DIRS ${PNG_PATH})
    set(PNG_LIBRARIES    "png")  # will result in libpng

    file(GLOB_RECURSE SRCS_PNG "${PNG_PATH}/*.c")
    add_library(${PNG_LIBRARIES} SHARED ${SRCS_PNG})
endif()

## -----------------------------------------------------------------------------
find_package(ZLIB)      # dont really need to check for zlib if png was success
if(ZLIB_FOUND)
    #message("-- zlib found: include_dir=(${ZLIB_INCLUDE_DIRS}), libs=(${ZLIB_LIBRARIES})")
else()
    set(ZLIB_PATH "${PROJECT_SOURCE_DIR}/ZLIB")
    message("-- zlib not found: will build from (${ZLIB_PATH})")

    set(ZLIB_INCLUDE_DIRS ${ZLIB_PATH})
    set(ZLIB_LIBRARIES    "z")  # will result in libz

    file(GLOB_RECURSE SRCS_ZLIB "${ZLIB_PATH}/*.c")
    add_library(${ZLIB_LIBRARIES} SHARED ${SRCS_ZLIB})
endif()

## -----------------------------------------------------------------------------
find_package(JPEG)
if(JPEG_FOUND)
    #message("-- jpeg lib found: include_dir=(${JPEG_INCLUDE_DIRS}), libs=(${JPEG_LIBRARIES})")
else()
    set(JPEG_PATH "${PROJECT_SOURCE_DIR}/JPEG")
    message("-- jpeg lib not found: will build from (${JPEG_PATH})")

    set(JPEG_INCLUDE_DIRS ${JPEG})
    set(JPEG_LIBRARIES    "jpeg")  # will result in libjpeg

    file(GLOB_RECURSE SRCS_JPEG "${JPEG_PATH}/*.c")
    add_library(${JPEG_LIBRARIES} SHARED ${SRCS_JPEG})
endif()


## -----------------------------------------------------------------------------
# all include paths
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/Src)
include_directories(${PNG_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS} ${JPEG_INCLUDE_DIRS})


## -----------------------------------------------------------------------------
## build poisson reconstruction library and executable
## -----------------------------------------------------------------------------
add_library(poisson_recon_lib SHARED
            ${PROJECT_SOURCE_DIR}/Src/PoissonReconLib.h
            ${PROJECT_SOURCE_DIR}/Src/PoissonReconLib.cpp)

target_link_libraries(poisson_recon_lib
                      ${PNG_LIBRARIES} ${ZLIB_LIBRARIES} ${JPEG_LIBRARIES})


add_executable(poisson_recon
               ${PROJECT_SOURCE_DIR}/Src/PoissonRecon.cpp)

target_link_libraries(poisson_recon
                      poisson_recon_lib
                      ${PNG_LIBRARIES} ${ZLIB_LIBRARIES} ${JPEG_LIBRARIES})

## -----------------------------------------------------------------------------
