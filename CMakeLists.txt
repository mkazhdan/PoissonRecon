cmake_minimum_required(VERSION 2.8)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

project( QPOISSON_RECON_LIB )

#to remove very annoying "NO ARGUMENTS THAT DEPEND ON A TEMPLATE PARAMETER" errors with G++!
if( UNIX OR MINGW )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -w")
elseif( MSVC )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

set( SAMPLE_PROJECT_VERSION_MAJOR 12 )
set( SAMPLE_PROJECT_VERSION_MINOR 0 )

file( GLOB header_list Src_CC_wrap/*.h )
file( GLOB source_list Src_CC_wrap/*.cpp )

add_library( ${PROJECT_NAME} STATIC ${header_list} ${source_list} )

# Add preprocessor definitions
set_property( TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS _CRT_SECURE_NO_DEPRECATE _CRT_SECURE_NO_WARNINGS NOMINMAX )

option( POISSON_RECON_WITH_OPEN_MP "Check to compile PoissonRecon plugin with OpenMP support" OFF )
if ( POISSON_RECON_WITH_OPEN_MP )
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set_property( TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS WITH_OPENMP )
	endif()
endif ()
