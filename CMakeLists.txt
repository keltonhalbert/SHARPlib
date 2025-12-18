cmake_minimum_required(VERSION 3.15)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GitVersion)

get_project_version()

configure_file(
    include/SHARPlib/version.h.in
    include/SHARPlib/version.h
)

project(SHARPlib VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH} LANGUAGES CXX)
include(GNUInstallDirs)
include(CTest)

# Options
option(BUILD_PYBIND "Build python bindings with nanobind" OFF)
option(NO_QC "Disable QC checks for missing/NaN values in SHARPlib" OFF)
set(NO_QC ${NO_QC} CACHE BOOL "Disable QC checks for missing/NaN values in SHARPlib tests" FORCE)

# Sources
file(GLOB SOURCES "src/SHARPlib/*.cpp" "src/SHARPlib/params/*.cpp")

add_library(SHARPlib STATIC ${SOURCES})

target_include_directories(SHARPlib PUBLIC include external/fmt/include)

set_target_properties(SHARPlib PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(SHARPlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(MSVC)
    target_compile_options(SHARPlib PRIVATE
        /W3 /WX
        /utf-8
        /wd4244 /wd4305
        /fp:precise 
    )
    
    target_compile_options(SHARPlib PRIVATE
        $<$<CONFIG:Debug>:/Zi /Od /MDd>
        $<$<CONFIG:Release>:/O2 /MD /DNDEBUG>
    )

elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    target_compile_options(SHARPlib PRIVATE
        -Wall -Wextra -Wpedantic -Werror
        -fp-model=precise 
        -ftree-vectorize
        -funroll-loops
    )

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(SHARPlib PRIVATE
        -Wall -Wextra -Wpedantic -Werror
        -ffp-model=precise 
        -ftree-vectorize
        -funroll-loops
    )

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(SHARPlib PRIVATE
        -Wall -Wextra -Wpedantic -Werror
        -ftree-vectorize
        -funroll-loops
    )

endif()

if(NOT MSVC)
    target_compile_options(SHARPlib PRIVATE
        $<$<CONFIG:Debug>:-g -O0> 
        $<$<CONFIG:Release>:-O3>
    )
endif()


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(NO_QC)
    target_compile_definitions(SHARPlib PRIVATE NO_QC)
endif()

if (NOT CMAKE_INSTALL_PREFIX) 
    set(LOCAL_PREFIX "${CMAKE_SOURCE_DIR}")
    set(CMAKE_INSTALL_PREFIX "${LOCAL_PREFIX}")
endif()

install(TARGETS SHARPlib DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    add_subdirectory(tests/unit)
    add_subdirectory(tests/benchmark)
endif()

if (BUILD_PYBIND)
    add_subdirectory(src/nanobind)
endif()

enable_testing()
