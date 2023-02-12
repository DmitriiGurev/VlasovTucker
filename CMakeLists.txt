cmake_minimum_required( VERSION 3.10 )

project( vlasov-t )

include_directories( ./include ./src ./include/eigen-3.4.0 ./libs/MshIO-main/include )

add_executable( tucker_test ./src/tucker.cpp ./test/tucker_test.cpp )
add_executable( mshio_test ./test/mshio_test.cpp )

add_library( libmshio STATIC IMPORTED )
set_property( TARGET libmshio PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/MshIO-main/build/libmshio.a )
target_link_libraries( mshio_test libmshio )