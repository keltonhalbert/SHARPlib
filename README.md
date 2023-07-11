# SHARP-calc
[![C++ CI](https://github.com/keltonhalbert/SHARP-calc/actions/workflows/cmake.yml/badge.svg)](https://github.com/keltonhalbert/SHARP-calc/actions/workflows/cmake.yml)


**Sounding and Hodograph Analysis and Research Program (SHARP)** C++ library for conducting analysis of atmospheric sounding profiles. Based on the NSHARP routines written by John Hart and Rich Thompson at the NWS Storm Prediction Center in Norman, Oklahoma. 

### C++ Version
Code is designed to target the **C++17** standard `-std=c++17`, and generally should not use **C++2X** features. This library is written to work on **RHEL8** systems, and with careful configuration, backwards to **RHEL7**, and should not attempt to use bleeding edge features at risk of breaking compatibility. The C++ 17 standard is not supported by the default RHEL7 compilers, but a version of GCC that supports the C++17 standard can be acquired through the **__devtoolset__** RHEL channel. 

### Style Guide
Though there are likely to be instances where it will need to be deviated from, this code generally attempts to abide by the [ISO C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines). While line width requirements are generally archaic, where possible we attempt to __keep line lengths to a maximum of 80 characters__ in order to preserve split screen code editing.  

Another note to the Style Guide is that, where possible/appropriate, full or verbose variable names are preferred to abbreviated ones when working with function parameters. For example, `temperature` or `pressure` is preferable to `temp` or `pres` when defining function arguments, so that it is abundantly clear to the code reader what is being passed through. This is especially the case with temperature, as `temp` is commonly used to refer to temporary variables, leading to confusion. 


### Unit Testing Framework
For unit tests, we make use of the [doctest singe header source library](https://github.com/doctest/doctest) found in the `tests` directory. In order to build and run the tests, execute the following commands from the project root directory:
```
mkdir build; cd build
cmake ..
make SHARPlb_tests
make test
```
NOTE: Right now, on Apple CLang, (and potentially CLang as a whole) some of the kinematics unit tests don't pass. Don't be alarmed if that's the case, it's just a difference in how the compilers treat certain things and will be remedied in the future. 

### Building the Static Library
To build the static library, simply run the following commands from the project root directory:
```
mkdir build; cd build
cmake .. 
make
make install
```
It will install the static library to PROJECT_ROOT/lib.

The library will be built in release mode by default, which sets the `-O3` optimization flag. To build in debug mode, clean out the build directory (CMake will cache certain things) and run the following command and proceed to build normally:
```
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

If you are planning on using the library with gridded data, you can turn off quality control checks traditionally only needed for observed sounding profile data. This can be disabled by running CMake with the following build flag:
```
cmake -DCMAKE_CXX_FLAGS="-DNO_QC" ..
```
This will also disable missing value checks in the tests directory.

If you would like a verbose compile process, run the make command in the following manner:
```
make VERBOSE=1
```

### Building the Python bindings
To build the python bindings, simply nagivate to the swig directory and install via pip. The python bindings require having SWIG installed, as well as a C++17 compatible compiler. The bindings will be installed in the following namespace tree: ```nwsspc.sharp.calc```.
```
cd swig
CC=/path/to/compilers/g++ pip install .
```


### Building the Docs
To build the HTML documentation pages, simply navigate your terminal to the `docs` directory and run:

```
cd docs
git submodule update --init --recursive ## download the CSS for the documentation
doxygen Doxyfile
```

This will generate the HTML pages using the docstring in the header files. Obviously, it requires that Doxygen be installed, [which can be found here.](https://doxygen.nl/) 

