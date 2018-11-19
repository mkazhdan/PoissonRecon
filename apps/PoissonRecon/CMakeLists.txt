cmake_minimum_required(VERSION 3.7)
project(Poisson)

set(PROJECT_SRCS
${PROJECT_SOURCE_DIR}/src/PoissonReconApp.cpp
)

source_group("src" FILES ${PROJECT_SRCS})

find_package( tinyxml2 REQUIRED )

add_executable(${PROJECT_NAME} ${PROJECT_SRCS})

target_link_libraries(${PROJECT_NAME} tinyxml2::tinyxml2 PoissonRecon)

set_property(TARGET ${PROJECT_NAME} PROPERTY FOLDER ${CMAKE_PROJECT_NAME}/apps)

install(TARGETS ${PROJECT_NAME}
        EXPORT ${CMAKE_PROJECT_NAME}Targets
        DESTINATION bin)
