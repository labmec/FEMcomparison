cmake_minimum_required (VERSION 3.13)

project (FEMcomparison)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#this allow us to define custom FindXXX.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(cmake/StandardProjectSettings.cmake)

## Sanitizer options if supported by compiler
#include(cmake/Sanitizers.cmake)
#add_library(project_options INTERFACE)
#enable_sanitizers(project_options)

# Allow for static analysis options
include(cmake/StaticAnalyzers.cmake)

# Enables NeoPZ library
include(cmake/EnableNeoPZ.cmake)
enable_pz()

# Enables log4cxx library
include(cmake/EnableLog4CXX.cmake)
enable_log4cxx()

find_package(Threads)

find_package(TBB)

#enabling MKL library
option(USING_MKL "Whether the MKL library will be linked in" OFF)
if(USING_MKL)
    include(cmake/EnableMKL.cmake)
    enable_mkl(Tools)
endif(USING_MKL)

# This option enables a lot of warnings and treat them as errors, to ensure
# good programming practices are used. Since its behaviour is extreme, it
# should be turned off by default.
option(ENABLE_ALL_WARNINGS "Enable a lot of compiler warnings and treat them as errors" OFF)
if (ENABLE_ALL_WARNINGS)
    add_library(project_warnings INTERFACE FEM/MeshInit.cpp FEM/MeshInit.h FEM/Solver.cpp FEM/Solver.h FEM/Output.cpp FEM/Output.h FEM/InputTreatment.cpp FEM/InputTreatment.h)
    include(cmake/CompilerWarnings.cmake)
    set_project_warnings(project_warnings)
    set(CMAKE_CXX_FLAGS "")
endif ()

# Add subdirectories
add_subdirectory(Tools)
add_subdirectory(Material)
add_subdirectory(FEM)
