#This file creates a library with tools
add_library(Tools
    DataStructure.h
    TPZCreateMultiphysicsSpace.cpp
    TPZCreateMultiphysicsSpace.h
    Tools.h
    Tools.cpp
)

target_link_libraries(Tools pz)
target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
