add_custom_target(pybind)

## Bring in the Python headers for the pybind11 bindings
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
## Bring in pybind11 from the symbolic link to our git submodule
add_subdirectory(pybind11)

include_directories(${Python3_INCLUDE_DIRS})
include_directories(pybind11/include)
include_directories(${CMAKE_SOURCE_DIR}/include)

pybind11_add_module(_thermo MODULE thermo.cpp)

target_compile_options(_thermo PRIVATE)
target_link_libraries(_thermo PRIVATE SHARPlib)

install(TARGETS _thermo DESTINATION ${CMAKE_INSTALL_PREFIX}/nwsspc/sharp/calc/)

