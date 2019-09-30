cmake_minimum_required(VERSION 3.15)

project(PoissonRecon VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-bigobj)
add_definitions(-DNOMINMAX)

find_package(OpenMP)
if (OPENMP_FOUND)
	add_definitions(-D_OPENMP)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(include ${PROJECT_SOURCE_DIR})

file(GLOB PNG_SRC PNG/*.c)
add_library(libpng ${PNG_SRC})

file(GLOB JPEG_SRC JPEG/*.cpp)
add_library(libjpeg ${JPEG_SRC})

file(GLOB ZLIB_SRC ZLIB/*.c)
add_library(zlib ${ZLIB_SRC})

add_executable(AdaptiveTreeVisualization Src/AdaptiveTreeVisualization.cpp)
target_link_libraries(AdaptiveTreeVisualization libpng libjpeg zlib)

add_executable(ChunkPLY Src/ChunkPLY.cpp)
target_link_libraries(ChunkPLY libpng libjpeg zlib)

add_executable(EDTInHeat Src/EDTInHeat.cpp)
target_link_libraries(EDTInHeat libpng libjpeg zlib)

add_executable(ImageStitching Src/ImageStitching.cpp)
target_link_libraries(ImageStitching libpng libjpeg zlib)

add_executable(PointInterpolant Src/PointInterpolant.cpp)
target_link_libraries(PointInterpolant libpng libjpeg zlib)

add_executable(PoissonRecon Src/PoissonRecon.cpp)
target_link_libraries(PoissonRecon libpng libjpeg zlib)

add_executable(SSDRecon Src/SSDRecon.cpp)
target_link_libraries(SSDRecon libpng libjpeg zlib)

add_executable(SurfaceTrimmer Src/SurfaceTrimmer.cpp)
target_link_libraries(SurfaceTrimmer libpng libjpeg zlib)
