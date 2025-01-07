cmake_minimum_required(VERSION 3.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(SHARPlib CXX)
include(GNUInstallDirs)
include(CTest)

## Useful macro for setting the CXX standard
macro(use_cxx17)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  else ()
    set (CMAKE_CXX_STANDARD 17)
  endif ()
endmacro(use_cxx17)

## This is where we set the CXX standard
use_cxx17()

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
    # Warning level 4 and treat warnings as errors
    # target_compile_options(SHARPlib PRIVATE /W4 /WX /utf-8)
    target_compile_options(SHARPlib PRIVATE /utf-8)

    # Debug and Release flags
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /MDd")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /MD /DNDEBUG")
else() # For GCC/Clang
    target_compile_options(SHARPlib PRIVATE -Wall -Wextra -Wpedantic -Werror)

    # Debug and Release flags
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O3 -ftree-vectorize -funroll-loops")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize -funroll-loops")
endif()

## Check if the build type is set and default to
## release, and then set the compiler flags
## for release and debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

## used for passing custom build options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(LOCAL_PREFIX "${CMAKE_SOURCE_DIR}")
set(CMAKE_INSTALL_PREFIX "${LOCAL_PREFIX}")

## install it to the PROJECT_ROOT/lib directory
install(TARGETS SHARPlib DESTINATION lib)


## This ensures that tests arent unintentionally added or built
## by other CMake projects that link to this library
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    add_subdirectory(tests/unit)
    add_subdirectory(tests/benchmark)
endif ()

enable_testing()
