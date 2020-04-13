cmake_minimum_required(VERSION 3.7)

project(PointInterpolant)

set(PROJECT_SRCS ${PROJECT_SOURCE_DIR}/src/PointInterpolant.cpp)

source_group("src" FILES ${PROJECT_SRCS})

find_package(OpenMP REQUIRED)

add_executable(${PROJECT_NAME} ${PROJECT_SRCS})

target_link_libraries(${PROJECT_NAME} PoissonRecon OpenMP::OpenMP_CXX)

set_property(TARGET ${PROJECT_NAME} PROPERTY FOLDER ${CMAKE_PROJECT_NAME}/apps)

install(TARGETS ${PROJECT_NAME} EXPORT ${CMAKE_PROJECT_NAME}Targets DESTINATION
                                                                    bin)
