cmake_minimum_required(VERSION 3.10)
add_compile_options(-mincoming-stack-boundary=2)

project(VlasovTucker)

include_directories(
    ./include
    ./src
    ./include/eigen-3.4.0
    ./libs/MshIO-main/include
    )

add_library(
    vlasov_tucker
    ./src/mesh.cpp
    ./src/poisson.cpp
    ./src/velocity_grid.cpp
    ./src/plasma_parameters.cpp
    ./src/solver.cpp
    ./src/vtk.cpp
    ./src/tucker.cpp
    )

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(vlasov_tucker PUBLIC OpenMP::OpenMP_CXX)
endif()

add_library(libmshio STATIC IMPORTED)
set_property(
    TARGET libmshio
    PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/MshIO-main/build/libmshio.a
    )

add_executable(run ./src/run.cpp)
target_link_libraries(run vlasov_tucker libmshio)

# Tests
add_subdirectory(./test)