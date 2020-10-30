# Project to apply error estimation for hybrid H1 approximations
add_executable(HybridH1vsMixed main_HybridH1vsMixed.cpp)
target_link_libraries(HybridH1vsMixed Methods Tools)

# Project to apply error estimation for hybrid H1 approximations
add_executable(HybridH1 main_HybridH1.cpp)
target_link_libraries(HybridH1 Methods Tools)

# Project to apply error estimation for hybrid H1 approximations
add_executable(Mixed main_Mixed.cpp)
target_link_libraries(Mixed Methods Tools)