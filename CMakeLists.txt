cmake_minimum_required( VERSION 3.10 )

project( vlasov-t )

include_directories( ./include ./src ./test ./include/eigen-3.4.0 ./libs/MshIO-main/include )

add_executable( tucker_test ./src/tucker.cpp ./test/tucker_test.cpp )
add_executable( mesh_test ./src/mesh.cpp ./test/mesh_test.cpp )
add_executable( poisson_test ./src/mesh.cpp ./src/poisson.cpp ./test/poisson_test.cpp ./src/vtk.cpp )

add_library( libmshio STATIC IMPORTED )
set_property( TARGET libmshio PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/MshIO-main/build/libmshio.a )
target_link_libraries( mesh_test libmshio )
target_link_libraries( poisson_test libmshio )