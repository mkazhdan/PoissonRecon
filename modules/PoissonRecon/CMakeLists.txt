cmake_minimum_required(VERSION 3.7)
project(PoissonRecon)

if(MSVC)
  set(DLLIMPORT "__declspec(dllimport)")
  set(DLLEXPORT "__declspec(dllexport)")
  set(DISABLE "__pragma(warning(push,0))")
  set(ENABLE "__pragma(warning(pop))")
endif()

set(PROJECT_SRCS
${PROJECT_SOURCE_DIR}/src/CmdLineParser.cpp
${PROJECT_SOURCE_DIR}/src/Factor.cpp
${PROJECT_SOURCE_DIR}/src/Geometry.cpp
${PROJECT_SOURCE_DIR}/src/MarchingCubes.cpp
${PROJECT_SOURCE_DIR}/src/PlyFile.cpp
)

set(PROJECT_INCS
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Allocator.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Array.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BinaryNode.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BSplineData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/CmdLineParser.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Factor.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FunctionData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Geometry.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MarchingCubes.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MAT.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MemoryUsage.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MultiGridOctreeData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MyTime.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Octree.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Ply.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PointStream.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Polynomial.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PPolynomial.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrix.h
)

source_group("include" FILES ${PROJECT_INCS})
source_group("src"     FILES ${PROJECT_SRCS})

add_library(${PROJECT_NAME} SHARED ${PROJECT_SRCS} ${PROJECT_INCS})

target_include_directories(${PROJECT_NAME} PUBLIC
                          $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                          $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)



set_property(TARGET ${PROJECT_NAME} PROPERTY
            FOLDER ${CMAKE_PROJECT_NAME}/modules)

install(TARGETS ${PROJECT_NAME}
        EXPORT ${CMAKE_PROJECT_NAME}Targets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

install(FILES 
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Allocator.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Array.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BinaryNode.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BSplineData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/CmdLineParser.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Factor.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FunctionData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Geometry.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MarchingCubes.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MAT.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MemoryUsage.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MultiGridOctreeData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MyTime.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Octree.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Ply.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PointStream.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Polynomial.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PPolynomial.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrix.h
        DESTINATION include/Mesh/PoissonRecon)

export(TARGETS ${PROJECT_NAME}
      APPEND FILE ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}Targets.cmake)

file(APPEND
    ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake
    "set(POISSONRECON_LIBRARIES \${POISSONRECON_LIBRARIES} ${PROJECT_NAME})\n")
