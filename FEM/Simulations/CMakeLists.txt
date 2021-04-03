
add_pz_target(
  NAME HybridH1vsMixed
SOURCES main_HybridH1vsMixed.cpp
)

# Project to apply error estimation for hybrid H1 approximations
#add_executable(HybridH1vsMixed main_HybridH1vsMixed.cpp)
target_link_libraries(HybridH1vsMixed PRIVATE Methods Tools)

# Project to apply error estimation for hybrid H1 approximations
add_pz_target(
  NAME HybridH1
SOURCES main_HybridH1.cpp)
#add_executable(HybridH1 main_HybridH1.cpp)
target_link_libraries(HybridH1 PRIVATE Methods Tools)

add_pz_target(
  NAME ContributeTest
SOURCES main_ContributeTest.cpp)
#add_executable(ContributeTest main_ContributeTest.cpp)
target_link_libraries(ContributeTest PRIVATE Methods Tools)

# Project to apply error estimation for hybrid H1 approximations
add_pz_target(
  NAME Mixed
SOURCES main_Mixed.cpp)

#add_executable(Mixed main_Mixed.cpp)
target_link_libraries(Mixed PRIVATE Methods Tools)
