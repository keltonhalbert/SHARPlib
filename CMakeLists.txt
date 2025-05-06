cmake_minimum_required(VERSION 3.15)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(SHARPlib VERSION 1.0.0 LANGUAGES CXX )
include(GNUInstallDirs)
include(CTest)

option(BUILD_PYBIND "Build python bindings with pybind11" OFF)
option(NO_QC "Disable QC checks for missing/NaN values in SHARPlib" OFF)
set(NO_QC ${NO_QC} CACHE BOOL "Disable QC checks for missing/NaN values in SHARPlib tests" FORCE)

## Bring the headers into the project
include_directories(include)
include_directories(external/fmt/include)

## add all CPP files as sources
file(GLOB SOURCES "src/SHARPlib/*.cpp" "src/SHARPlib/params/*.cpp")

## Add the library build command
## The library target is SHARP, we are building a static lib,
## and the sources are proviced
add_library(SHARPlib STATIC ${SOURCES})

## We are telling the build process that the SHARP library
## should use the C++ linker, and setting compile options
set_target_properties(SHARPlib PROPERTIES LINKER_LANGUAGE CXX)

## Compiler-specific options
if (MSVC) # For MSVC (cl) compiler
    # Add compile options for MSVC
    target_compile_options(SHARPlib PRIVATE
        /W3           # Warning level 3
        /WX           # Treat warnings as errors
        /utf-8        # Set source and execution character set to UTF-8
        /wd4244       # Disable warning C4244 (possible loss of data)
        /wd4305       # Disable warning C4305 (truncation from double to float)
    )

    # Add debug and release specific options
    target_compile_options(SHARPlib PRIVATE
        $<$<CONFIG:Debug>:/Zi /Od /MDd>      # Debug: Enable debug info, disable optimization, debug runtime
        $<$<CONFIG:Release>:/O2 /MD /DNDEBUG> # Release: Optimize for speed, release runtime, disable assertions
    )
else() # For GCC/Clang
    # Add compile options for GCC/Clang
    target_compile_options(SHARPlib PRIVATE
        -Wall           # Enable common warnings
        -Wextra         # Enable extra warnings
        -Wpedantic      # Enforce strict compliance with standard
        -Werror         # Treat warnings as errors
    )

    # Add debug and release specific options
    target_compile_options(SHARPlib PRIVATE
        $<$<CONFIG:Debug>:-g -O3 -ftree-vectorize -funroll-loops> # Debug: Enable debug info, optimizations, vectorization
        $<$<CONFIG:Release>:-O3 -ftree-vectorize -funroll-loops -fPIC> # Release: Optimize, vectorize, generate PIC
    )
endif()

## Check if the build type is set and default to
## release, and then set the compiler flags
## for release and debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(NO_QC)
    target_compile_definitions(SHARPlib PRIVATE NO_QC)
endif()

## If an install prefix hasn't been set, use the 
## current source directory
if (NOT CMAKE_INSTALL_PREFIX) 
    set(LOCAL_PREFIX "${CMAKE_SOURCE_DIR}")
    set(CMAKE_INSTALL_PREFIX "${LOCAL_PREFIX}")
endif()

install(TARGETS SHARPlib DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

## This ensures that tests and python bindings arent 
## unintentionally added or built by other CMake projects
## that link to this library
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    add_subdirectory(tests/unit)
    add_subdirectory(tests/benchmark)
endif()

if (BUILD_PYBIND)
    add_subdirectory(src/nanobind)
endif()

enable_testing()
