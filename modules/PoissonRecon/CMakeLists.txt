cmake_minimum_required(VERSION 3.7)
project(PoissonRecon)

set(PROJECT_SRCS
${PROJECT_SOURCE_DIR}/src/PlyFile.cpp
)

set(PROJECT_INCS

${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Allocator.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Array.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Array.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BinaryNode.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BMPStream.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BSplineData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BSplineData.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/CmdLineParser.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/CmdLineParser.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Factor.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.Evaluation.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.Initialize.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.IsoSurface.specialized.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.SortedTreeNodes.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.System.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.WeightedSamples.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FunctionData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FunctionData.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Geometry.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Geometry.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Image.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/JPEG.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/JPEG.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/LinearSolvers.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MarchingCubes.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MAT.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MAT.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MyMiscellany.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Ply.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PlyFile.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PNG.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PNG.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PointStream.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PointStream.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PointStreamData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Polynomial.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Polynomial.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PPolynomial.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PPolynomial.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/RegularTree.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/RegularTree.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrix.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrix.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrixInterface.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrixInterface.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Window.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Window.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PoissonRecon.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SurfaceTrimmer.h

)

source_group("include" FILES ${PROJECT_INCS})
source_group("src"     FILES ${PROJECT_SRCS})


if((NOT APPLE) OR (NOT CMAKE_COMPILER_IS_GNUCXX) OR (GCC_VERSION VERSION_GREATER 4.2.1) OR (CMAKE_COMPILER_IS_CLANG))
  find_package(OpenMP)
endif()
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  message (STATUS "Found OpenMP")
  if(MSVC)
    if(MSVC_VERSION EQUAL 1500)
      set(OPENMP_DLL VCOMP90)
    elseif(MSVC_VERSION EQUAL 1600)
      set(OPENMP_DLL VCOMP100)
    elseif(MSVC_VERSION EQUAL 1700)
      set(OPENMP_DLL VCOMP110)
    elseif(MSVC_VERSION EQUAL 1800)
      set(OPENMP_DLL VCOMP120)
    elseif(MSVC_VERSION EQUAL 1900)
      set(OPENMP_DLL VCOMP140)
    elseif(MSVC_VERSION EQUAL 1910)
      set(OPENMP_DLL VCOMP140)
    endif()
    if(OPENMP_DLL)
      set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /DELAYLOAD:${OPENMP_DLL}D.dll")
      set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DELAYLOAD:${OPENMP_DLL}.dll")
    else(OPENMP_DLL)
      message(WARNING "Delay loading flag for OpenMP DLL is invalid.")
    endif(OPENMP_DLL)
  endif(MSVC)
else(OPENMP_FOUND)
  message (STATUS "Not found OpenMP")
endif()

add_library(${PROJECT_NAME} ${PROJECT_SRCS} ${PROJECT_INCS})

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

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
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Array.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BinaryNode.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BMPStream.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BSplineData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/BSplineData.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/CmdLineParser.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/CmdLineParser.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Factor.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.Evaluation.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.Initialize.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.IsoSurface.specialized.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.SortedTreeNodes.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.System.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FEMTree.WeightedSamples.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FunctionData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/FunctionData.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Geometry.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Geometry.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Image.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/JPEG.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/JPEG.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/LinearSolvers.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MarchingCubes.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MAT.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MAT.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/MyMiscellany.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Ply.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PlyFile.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PNG.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PNG.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PointStream.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PointStream.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PointStreamData.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Polynomial.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Polynomial.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PPolynomial.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PPolynomial.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/RegularTree.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/RegularTree.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrix.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrix.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrixInterface.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SparseMatrixInterface.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Window.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/Window.inl
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/PoissonRecon.h
${PROJECT_SOURCE_DIR}/include/Mesh/PoissonRecon/SurfaceTrimmer.h
        DESTINATION include/Mesh/PoissonRecon)

export(TARGETS ${PROJECT_NAME}
      APPEND FILE ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}Targets.cmake)

file(APPEND
    ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake
    "set(POISSONRECON_LIBRARIES \${POISSONRECON_LIBRARIES} ${PROJECT_NAME})\n")
