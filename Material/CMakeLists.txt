# Adds materials to ErrorEstimationLib
target_sources(Tools PUBLIC
    LCC_MatLaplacianHybrid.cpp
    LCC_MatLaplacianHybrid.h
    LCC_TPZMixedDarcyFlow.h
    LCC_TPZMixedDarcyFlow.cpp
    LCC_LagrangeMultiplier.h
    LCC_LagrangeMultiplier.cpp
)

target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_LIST_DIR})
