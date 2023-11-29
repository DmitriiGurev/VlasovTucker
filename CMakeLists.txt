cmake_minimum_required(VERSION 3.10)

project(VlasovTucker)

include_directories(
    ./include
    ./src
    ./include/eigen-3.4.0
    ./libs/MshIO-main/include
    )

find_package(OpenMP REQUIRED)

add_library(
    vlasov_tucker
    ./src/mesh.cpp
    ./src/poisson.cpp
    ./src/velocity_grid.cpp
    ./src/plasma_parameters.cpp
    ./src/solver.cpp
    ./src/vtk.cpp
    ./src/smoother.cpp
    # TODO: Add tucker.cpp
    )

add_library(libmshio STATIC IMPORTED)
set_property(
    TARGET libmshio
    PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/MshIO-main/build/libmshio.a
    )

add_executable(run ./src/run.cpp)
target_link_libraries(run vlasov_tucker libmshio)

# Tests
add_subdirectory(./test)