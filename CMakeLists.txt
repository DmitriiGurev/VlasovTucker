cmake_minimum_required(VERSION 3.10)

project(VlasovTucker)

# add_compile_options(-O3 -DNDEBUG -fno-math-errno -march=native)
add_compile_options(-O3 -fno-math-errno -march=native)
add_compile_options(-mincoming-stack-boundary=2)

include_directories(
    ./include
    ./src
    ./include/eigen-3.4.0
    ./libs/MshIO-main/include
    )

add_library(
    vlasov_tucker
    ./src/primitives.cpp
    ./src/mesh.cpp
    ./src/poisson.cpp
    ./src/velocity_grid.cpp
    ./src/particle_data.cpp
    ./src/solver.cpp
    ./src/multicomponent_solver.cpp
    ./src/vtk.cpp
    ./src/full.cpp
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

# Tests
add_subdirectory(./test)

# Examples
add_subdirectory(./examples)