cmake_minimum_required(VERSION 2.8.12.2)
project(SHARP)
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
include_directories(external/fmt/include EXCLUDE_FROM_ALL)
#add_subdirectory(external/fmt EXCLUDE_FROM_ALL)

## add all CPP files as sources
file(GLOB SOURCES "src/SHARPlib/*.cpp")

## Add the library build command
## The library target is SHARP, we are building a static lib,
## and the sources are proviced
add_library(SHARP STATIC ${SOURCES})

##target_link_libraries(SHARP PRIVATE fmt)

## We are telling the build process that the SHARP library
## should use the C++ linker, and setting compile options
set_target_properties(SHARP PROPERTIES LINKER_LANGUAGE CXX)
add_compile_options(-Wall -Wextra -Wpedantic)

## Check if the build type is set and default to
## release, and then set the compiler flags
## for release and debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -O3 -Wall -ftree-vectorize -funroll-loops -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -ffast-math -ftree-vectorize -funroll-loops -march=native")

## used for passing custom build options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(LOCAL_PREFIX "${CMAKE_SOURCE_DIR}")
set(CMAKE_INSTALL_PREFIX "${LOCAL_PREFIX}")

## install it to the PROJECT_ROOT/lib directory
install(TARGETS SHARP DESTINATION lib)

## add unit test directory
add_subdirectory(tests)
add_test(NAME sharp_thermo COMMAND test_sharp_thermo)
add_test(NAME sharp_interp COMMAND test_sharp_interp)
add_test(NAME sharp_parcel COMMAND test_sharp_parcel)
add_test(NAME sharp_layer COMMAND test_sharp_layer)
add_test(NAME sharp_winds COMMAND test_sharp_winds)
add_test(NAME sharp_algs COMMAND test_sharp_algs)

