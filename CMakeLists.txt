cmake_minimum_required(VERSION 3.20)
project(PoissonRecon)

option(BUILD_SHARED_LIB "Builds PoissonRecon as a shared library artifact" OFF)
option(BUILD_EXECUTABLE "Builds PoissonRecon as a command-line application" ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Threads REQUIRED)

if (MSVC)
    string(JOIN " " CMAKE_CXX_FLAGS_RELEASE
        ${CMAKE_CXX_FLAGS_RELEASE}
        -DRELEASE
        -DNOMINMAX
        -DWIN32_LEAN_AND_MEAN
        /bigobj
    )
else()
    string(JOIN " " CMAKE_CXX_FLAGS_RELEASE
        ${CMAKE_CXX_FLAGS_RELEASE}
        -DRELEASE
        -Ofast
        -funroll-loops
        -ffast-math
        -Wno-deprecated
        -Wno-invalid-offsetof
        -Wno-dangling-else
    )
endif()

string(JOIN " " CMAKE_EXE_LINKER_FLAGS_RELEASE
    ${CMAKE_EXE_LINKER_FLAGS_RELEASE}
    ${CMAKE_CXX_FLAGS_RELEASE}
)

file(GLOB_RECURSE SRCS_ZLIB "ZLIB/*.c")
add_library(ZLIB STATIC ${SRCS_ZLIB})

file(GLOB_RECURSE SRCS_PNG "PNG/*.c")
add_library(PNG STATIC ${SRCS_PNG})

file(GLOB_RECURSE SRCS_JPEG "JPEG/*.cpp")
add_library(JPEG STATIC ${SRCS_JPEG})

include_directories(${CMAKE_CURRENT_LIST_DIR} JPEG ZLIB PNG)

if (BUILD_SHARED_LIB)
    message(STATUS "PoissonRecon: Building shared library")
    add_library(PoissonRecon SHARED Src/PoissonRecon.cpp)
    target_link_libraries(PoissonRecon PRIVATE PNG ZLIB JPEG Threads::Threads)
    target_compile_definitions(PoissonRecon PRIVATE -Dmain=PoissonRecon)
    install(
        TARGETS PoissonRecon
        EXPORT PoissonRecon
        INCLUDES DESTINATION include
        FRAMEWORK DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
    install(
        EXPORT PoissonRecon
        NAMESPACE PoissonRecon::
        DESTINATION lib/cmake/PoissonRecon
    )
endif()

if (BUILD_EXECUTABLE)
    message(STATUS "PoissonRecon: Building executable")
    add_executable(PoissonReconApp Src/PoissonRecon.cpp)
    target_link_libraries(PoissonReconApp PRIVATE PNG ZLIB JPEG Threads::Threads)
    install(TARGETS PoissonReconApp)
endif()
