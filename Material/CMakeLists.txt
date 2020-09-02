# Adds materials to ErrorEstimationLib
target_sources(Tools PUBLIC
    TPZMatLaplacianHybrid.cpp
    TPZMatLaplacianHybrid.h
)

target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_LIST_DIR})
