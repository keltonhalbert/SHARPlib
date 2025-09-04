add_custom_target(pybind)

if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

## Bring in the Python headers for the pybind11 bindings
find_package(Python 3.10 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_subdirectory(nanobind)

include_directories(${Python3_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/external/fmt/include)
# for version.h
include_directories("${CMAKE_BINARY_DIR}/include")


nanobind_add_module(calc LTO calc.cpp)
target_link_libraries(calc PRIVATE SHARPlib)
if (MSVC) # For MSVC (cl) compiler
    # Add compile options for MSVC
    target_compile_options(calc PRIVATE
        /utf-8        # Set source and execution character set to UTF-8
    )
endif()
set_target_properties(calc PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp" # add this too
)

# Also override per-config on MSVC
foreach(cfg IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
    set_target_properties(calc PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY_${cfg} "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp"
        RUNTIME_OUTPUT_DIRECTORY_${cfg} "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp"
        ARCHIVE_OUTPUT_DIRECTORY_${cfg} "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp"
    )
endforeach()

nanobind_add_stub(
    calc_stub
    RECURSIVE
    MODULE nwsspc.sharp.calc
    OUTPUT 
        __init__.pyi 
        interp.pyi 
        layer.pyi
        params.pyi 
        parcel.pyi
        thermo.pyi
        winds.pyi
    OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp"
    PYTHON_PATH "${CMAKE_CURRENT_BINARY_DIR}"
    MARKER_FILE "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp/py.typed"
    DEPENDS calc
)

install(TARGETS calc DESTINATION nwsspc/sharp)

install(
    DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/nwsspc/sharp/"
    DESTINATION "nwsspc/sharp"
    FILES_MATCHING
        PATTERN "*.pyi"
        PATTERN "py.typed"
)
