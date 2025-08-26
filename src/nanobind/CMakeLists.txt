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


nanobind_add_module(calc LTO calc.cpp)
nanobind_add_stub(
    calc_stub
    MODULE calc
    RECURSIVE
    OUTPUT "calc" 
    PYTHON_PATH "$<TARGET_FILE_DIR:calc>"
    MARKER_FILE "py.typed"
    DEPENDS calc
)
if (MSVC) # For MSVC (cl) compiler
    # Add compile options for MSVC
    target_compile_options(calc PRIVATE
        /utf-8        # Set source and execution character set to UTF-8
    )
endif()
target_link_libraries(calc PRIVATE SHARPlib)
install(TARGETS calc DESTINATION nwsspc/sharp)
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/py.typed"
    DESTINATION "nwsspc/sharp"
)
install(
    DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/calc/"
    DESTINATION "nwsspc/sharp"
)
