# SHARPlib
[![C++ CI (Linux, MacOS, Windows)](https://github.com/keltonhalbert/SHARPlib/actions/workflows/cmake.yml/badge.svg)](https://github.com/keltonhalbert/SHARPlib/actions/workflows/cmake.yml)
[![Python CI (Linux, MacOS, Windows)](https://github.com/keltonhalbert/SHARPlib/actions/workflows/python.yml/badge.svg?branch=develop)](https://github.com/keltonhalbert/SHARPlib/actions/workflows/python.yml)
[![Build Docs](https://github.com/keltonhalbert/SHARPlib/actions/workflows/doxygen-gh-pages.yml/badge.svg)](https://github.com/keltonhalbert/SHARPlib/actions/workflows/doxygen-gh-pages.yml)


**Sounding and Hodograph Analysis and Research Program (SHARP)** C++ library for conducting analysis of atmospheric sounding profiles. Based on the NSHARP routines written by John Hart and Rich Thompson at the NWS Storm Prediction Center in Norman, Oklahoma. 

### C++ Version
Code is designed to target the **C++17** standard `-std=c++17`, and generally should not use **C++2X** features. This library is written to work on **RHEL8** systems, and with careful configuration, backwards to **RHEL7**, and should not attempt to use bleeding edge features at risk of breaking compatibility. The C++ 17 standard is not supported by the default RHEL7 compilers, but a version of GCC that supports the C++17 standard can be acquired through the **__devtoolset__** RHEL channel. 

### Style Guide
Though there are likely to be instances where it will need to be deviated from, this code generally attempts to abide by the [ISO C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines). While line width requirements are generally archaic, where possible we attempt to __keep line lengths to a maximum of 80 characters__ in order to preserve split screen code editing.  

Another note to the Style Guide is that, where possible/appropriate, full or verbose variable names are preferred to abbreviated ones when working with function parameters. For example, `temperature` or `pressure` is preferable to `temp` or `pres` when defining function arguments, so that it is abundantly clear to the code reader what is being passed through. This is especially the case with temperature, as `temp` is commonly used to refer to temporary variables, leading to confusion. 

### BEFORE YOU BUILD
SHARPlib has some light-weight dependencies for testing, benchmarking, documentation building, string formatting, and creating python bindings. These can easily be downloaded for building by runing the following command to download the dependencies from GitHub over SSH:
```bash
git submodule update --init --recursive 
```

### Building SHARPlib (C++) 
To build SHARPlib, execute the following commands in the project root directory:
```bash
cmake -B build .
cmake --build build -j N_BUILD_PROCESSES
```

This will build the static library in the ```{$POJECT_ROOT}/build``` directory in parallel with N_BUILD_PROCESSES. This isn't terribly useful by itself, so to install the static library, you can execute: 
```bash
cmake -B build . --install-prefix=/path/to/where/you/want/SHARPlib
cmake --build build -j N_BUILD_PROCESSES
cmake --install build
```

If you wish to create a debug build, simply pass the following arguments to CMake:
```bash
cmake --build build -j N_BUILD_PROCESSES --config Debug
```

### Building SHARPlib (Python)
SHARPlib C++ code can be called from Python using [nanobind](https://github.com/wjakob/nanobind) to handle the wrapping.
Building SHARPlib with its Python bindings is quite easy-- you can simply install it via pip from the current directory: 
```bash
pip install .
```

If you desire to manually build the SHARPlib library + python bindings, you may execute:
```bash
cmake -B build . -DBUILD_PYBIND=ON
cmake --build build -j N_BUILD_PROCESSES
```


### Unit Testing Framework
For unit tests, we make use of the [doctest singe header source library](https://github.com/doctest/doctest) found in the `tests` directory. In order to build and run the tests, execute the following commands from the project root directory:

```bash
cmake -B build . 
cmake --build build -j N_BUILD_PROCESSES --target SHARPlib_tests
ctest --test-dir build
```


### Building the Documentation
To build the HTML documentation pages, the following can be executed in the terminal from the project root: 
```bash
doxygen Doxyfile
```

This will generate the HTML pages using the docstring in the header files. Obviously, it requires that Doxygen be installed, [which can be found here.](https://doxygen.nl/) 

Documentation is automatically built on push, [and can be found here](https://keltonhalbert.github.io/SHARPlib/)
