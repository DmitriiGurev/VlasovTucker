cmake_minimum_required( VERSION 3.10 )
# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

project( vlasov-t )

include_directories( ./include ./src ./include/eigen-3.4.0 ./libs/MshIO-main/include )

add_executable( tucker_test ./src/tucker.cpp ./test/tucker_test.cpp )

# add_subdirectory( ./libs/MshIO-main )
# target_link_libraries( vlasov_t mshio )

# add_library(libmshio STATIC IMPORTED)
# set_property(TARGET libmshio PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/MshIO-main/build/libmshio.a)
# target_link_libraries( vlasov_t libmshio )

# file(GLOB_RECURSE SRC_FILES src/*.cpp)
# add_executable(vlasov_t ${SRC_FILES})

# target_include_directories(vlasov_t PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# add_library(libmshio.a STATIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/MshIO-main/build)
# set_target_properties(libmshio.a PROPERTIES LINKER_LANGUAGE CXX)