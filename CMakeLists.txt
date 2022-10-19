cmake_minimum_required(VERSION 2.8.12.2)
project(SHARP)

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

## add all CPP files as sources
file(GLOB SOURCES "src/sharptab/*.cpp")

## Add the library build command
add_library(SHARP STATIC ${SOURCES})
set_target_properties(SHARP PROPERTIES LINKER_LANGUAGE CXX)


## Make the output from the makefiles verbose
#set(CMAKE_VERBOSE_MAKEFILE ON)

## set install prefix to ../lib
set(LOCAL_PREFIX "../")
set(CMAKE_INSTALL_PREFIX "${LOCAL_PREFIX}")

## install it to the local lib directory
install(TARGETS SHARP DESTINATION lib)

## add unit test directory
add_subdirectory(tests)
