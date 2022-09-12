target_sources(Tools PUBLIC
        TPZCreateHybridizedMixedSpaces.cpp
        TPZCreateHybridizedMixedSpaces.h
        TPZCreateMultiphysicsSpace.cpp
        TPZCreateMultiphysicsSpace.h
        )

target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_LIST_DIR})