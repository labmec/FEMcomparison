# Adds materials to ErrorEstimationLib
target_sources(Tools PUBLIC
    TPZMatLaplacianHybrid.cpp
    TPZMatLaplacianHybrid.h
    LCC_MatLaplacianHybrid.cpp
    LCC_MatLaplacianHybrid.h
    LCC_MixedPoisson.h
    LCC_MixedPoisson.cpp
)

target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_LIST_DIR})
