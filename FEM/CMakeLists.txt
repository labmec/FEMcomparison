# Add subdirectories
add_library(Methods
        InputTreatment.h
        InputTreatment.cpp
        MeshInit.h
        MeshInit.cpp
        Output.cpp
        Output.h
        Solver.cpp
        Solver.h
        )

target_link_libraries(Methods pz Tools)
target_include_directories(Methods PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(Simulations)

