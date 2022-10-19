# SHARP-calc
**Sounding and Hodograph Analysis and Research Program (SHARP)** C++ library for conducting analysis of atmospheric sounding profiles. Based on the NSHARP routines written by John Hart and Rich Thompson at the NWS Storm Prediction Center in Norman, Oklahoma. 

### C++ Version
Code is designed to target the **C++17** standard `-std=c++17`, and generally should not use **C++2X** features. This library is written to work on **RHEL8** systems, and with careful configuration, backwards to **RHEL7**, and should not attempt to use bleeding edge features at risk of breaking compatibility. 

### Style Guide
Though there are likely to be instances where it will need to be deviated from, this code generally attempts to abide by the [Google C++ style guide](https://google.github.io/styleguide/cppguide.html). While line width requirements are generally archaic, where possible we attempt to __keep line lengths to a maximum of 80 characters__ in order to preserve split screen code editing.  

Another note to the Style Guide is that, where possible/appropriate, full or verbose variable names are preferred to abbreviated ones when working with function parameters. For example, `temperature` or `pressure` is preferable to `temp` or `pres` when defining function arguments, so that it is abundantly clear to the code reader what is being passed through. This is especially the case with temperature, as `temp` is commonly used to refer to temporary variables, leading to confusion. 

### Building the Static Library
To build the static library, simply run the following commands from the project root directory:
```
mkdir build; cd build
cmake .. 
make
```
It will install the static library to PROJECT_ROOT/lib

### Building the Docs
To build the HTML documentation pages, simply navigate your terminal to the `docs` directory and run:

```doxygen Doxyfile```

This will generate the HTML pages using the docstring in the header files. Obviously, it requires that Doxygen be installed, [which can be found here.](https://doxygen.nl/) 

