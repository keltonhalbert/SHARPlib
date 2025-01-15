add_custom_target(pybind)

if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

## Bring in the Python headers for the pybind11 bindings
find_package(Python3 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_subdirectory(nanobind)

include_directories(${Python3_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/include)

nanobind_add_module(thermo LTO thermo.cpp)
nanobind_add_module(layer LTO layer.cpp)
nanobind_add_module(interp LTO interp.cpp)

target_link_libraries(thermo PRIVATE SHARPlib)
target_link_libraries(layer PRIVATE SHARPlib)
target_link_libraries(interp PRIVATE SHARPlib)

install(TARGETS thermo layer interp DESTINATION nwsspc/sharp/calc)

