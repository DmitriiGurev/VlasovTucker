add_executable(oscillations ./oscillations.cpp)
target_link_libraries(oscillations vlasov_tucker libmshio)

add_executable(sheath ./sheath.cpp)
target_link_libraries(sheath vlasov_tucker libmshio)