cmake_minimum_required(VERSION 2.8.9)
project(SHARPlib)

## Bring the headers into the project
include_directories(include)

## add all CPP files as sources
file(GLOB SOURCES "src/sharptab/*.cpp")

## Add the library build command
add_library(SHARP STATIC ${SOURCES})

## set install prefix to ../lib
set(LOCAL_PREFIX "../")
set(CMAKE_INSTALL_PREFIX "${LOCAL_PREFIX}")
## install it to the local lib directory
install(TARGETS SHARP DESTINATION lib)
