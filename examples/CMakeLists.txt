add_executable(oscillations ./oscillations.cpp)
target_link_libraries(oscillations vlasov_tucker libmshio)

add_executable(landau_damping ./landau_damping.cpp)
target_link_libraries(landau_damping vlasov_tucker libmshio)

add_executable(instability ./instability.cpp)
target_link_libraries(instability vlasov_tucker libmshio)

add_executable(bump_on_tail ./bump_on_tail.cpp)
target_link_libraries(bump_on_tail vlasov_tucker libmshio)

add_executable(sheath ./sheath.cpp)
target_link_libraries(sheath vlasov_tucker libmshio)