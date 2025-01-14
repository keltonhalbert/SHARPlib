add_custom_target(pybind)

## Bring in the Python headers for the pybind11 bindings
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
# find_package(pybind11 CONFIG REQUIRED)
## Bring in pybind11 from the symbolic link to our git submodule
add_subdirectory(pybind11)

include_directories(${Python3_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/include)

pybind11_add_module(thermo MODULE thermo.cpp)
pybind11_add_module(layer MODULE layer.cpp)

target_link_libraries(thermo PRIVATE SHARPlib pybind11::headers)
target_link_libraries(layer PRIVATE SHARPlib pybind11::headers)

install(TARGETS thermo layer DESTINATION nwsspc/sharp/calc)

