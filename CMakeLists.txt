cmake_minimum_required( VERSION 3.10 )

project( QPOISSON_RECON_LIB )

add_library( ${PROJECT_NAME} STATIC )

# To remove very annoying "NO ARGUMENTS THAT DEPEND ON A TEMPLATE PARAMETER" errors with G++!
if( UNIX OR MINGW )
	target_compile_options( ${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-fpermissive -w> )
elseif( MSVC )
	target_compile_options( ${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj> )
endif()

target_sources( ${PROJECT_NAME}
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/Src_CC_wrap/PoissonReconLib.h
		${CMAKE_CURRENT_SOURCE_DIR}/Src_CC_wrap/PointData.h
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/Src_CC_wrap/PoissonReconLib.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Src_CC_wrap/PointData.cpp
)

target_include_directories( ${PROJECT_NAME}
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/Src_CC_wrap
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/Src
)

target_compile_definitions( ${PROJECT_NAME}
	PRIVATE
		_CRT_SECURE_NO_DEPRECATE
		_CRT_SECURE_NO_WARNINGS
		NOMINMAX
)

option( POISSON_RECON_WITH_OPEN_MP "Compile PoissonRecon plugin with OpenMP support" OFF )
if( POISSON_RECON_WITH_OPEN_MP )
	find_package( OpenMP )
	
	if( OpenMP_CXX_FOUND ) 
		target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX )
		target_compile_definitions( ${PROJECT_NAME} PRIVATE WITH_OPENMP )
	endif()
endif ()
