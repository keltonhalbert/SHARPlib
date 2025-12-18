cmake_minimum_required(VERSION 3.15)
project(SHARPLIB LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GitVersion)
get_project_version()

configure_file(
    include/SHARPlib/version.h.in
    include/SHARPlib/version.h
)

if(NOT PROJECT_VERSION)
    set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
endif()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(CTest)

option(BUILD_PYBIND "Build python bindings with nanobind" OFF)
option(NO_QC "Disable QC checks for missing/NaN values in SHARPlib" OFF)
option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON) 
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS 
    "src/SHARPlib/*.cpp" 
)

add_library(SHARPlib ${SOURCES})

target_include_directories(SHARPlib PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/fmt/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(SHARPlib PUBLIC cxx_std_17)
set_target_properties(SHARPlib PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
)

if(MSVC)
    target_compile_options(SHARPlib PRIVATE
        /W3 /utf-8 /fp:precise
        /wd4244 /wd4305
        $<$<BOOL:${WARNINGS_AS_ERRORS}>:/WX>
    )
else()
    target_compile_options(SHARPlib PRIVATE
        -Wall -Wextra -Wpedantic
        $<$<BOOL:${WARNINGS_AS_ERRORS}>:-Werror>
        $<$<CONFIG:Release>:-ftree-vectorize -funroll-loops>
    )
    
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")
        target_compile_options(SHARPlib PRIVATE -ffp-model=precise)
    endif()
endif()

if(NO_QC)
    target_compile_definitions(SHARPlib PRIVATE NO_QC)
endif()

install(TARGETS SHARPlib
    EXPORT SHARPlibTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/SHARPlib 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT SHARPlibTargets
    FILE SHARPlibTargets.cmake
    NAMESPACE SHARPlib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SHARPlib
)

set(PKG_CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/SHARPlibConfig.cmake")
set(PKG_VER_FILE    "${CMAKE_CURRENT_BINARY_DIR}/SHARPlibConfigVersion.cmake")

write_basic_package_version_file(
    "${PKG_VER_FILE}"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SHARPlibConfig.cmake.in"
    "${PKG_CONFIG_FILE}"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SHARPlib
)

install(FILES 
    "${PKG_CONFIG_FILE}" 
    "${PKG_VER_FILE}"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SHARPlib
)

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
    enable_testing()
    add_subdirectory(tests/unit)
    add_subdirectory(tests/benchmark)
endif()

if (BUILD_PYBIND)
    add_subdirectory(src/nanobind)
endif()
