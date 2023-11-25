# add_executable(mesh_test ./mesh_test.cpp)
# target_link_libraries(mesh_test vlasov_tucker libmshio)

add_executable(poisson_test ./poisson_test.cpp)
target_link_libraries(poisson_test vlasov_tucker libmshio)

# add_executable(tucker_test ../src/tucker.cpp ./tucker_test.cpp)