Build and Install SHARPlib (C++)
================================

.. _cpp_dependencies:

Dependencies
------------
The only external dependencies are a compatible C++ compiler, and CMake. 

.. _cpp_version_req:

C++ Version
-----------
Code is designed to target the **C++17** standard ``-std=c++17``, and generally should not use **C++2X** features. This library is written to work on **RHEL8** systems, and with careful configuration, backwards to **RHEL7**, and should not attempt to use bleeding edge features at risk of breaking compatibility. The C++ 17 standard is not supported by the default RHEL7 compilers, but a version of GCC that supports the C++17 standard can be acquired through the **devtoolset** RHEL channel.

.. _cloning:

Clone the Repository
--------------------
The best way to get the source code is through the git CLI. Downloading the zip file from the repository will not work due to some light-weight dependencies for testing, benchmarking, documentation, string formatting, and python binding generation that are managed through git submodules. To clone the repository:

.. code-block:: bash 

    ## if using https
    git clone https://github.com/keltonhalbert/SHARPlib.git 

    ## if using SSH
    git@github.com:keltonhalbert/SHARPlib.git
    
    ## clone the submodule dependencies
    git submodule update --init --recursive

This will clone the repository into a directory called SHARPlib that is a child of whichever parent directory you executed the clone in.

.. _building:

Building and Installing SHARPlib
--------------------------------
As long as you have a C++17 compatible compiler (generally, GCC >= 8), building is rather straightforward. To build SHARPlib, execute the following commands in the project root (SHARPlib) directory:

.. code-block:: bash

   cmake -B build .
   cmake --build build -j N_BUILD_PROCESSES

This will build the static library in the ``{$POJECT_ROOT}/build`` directory in parallel with N_BUILD_PROCESSES. This isn't terribly useful by itself, so to install the static library, you can execute:

.. code-block:: bash

   cmake -B build . --install-prefix=/path/to/where/you/want/SHARPlib
   cmake --build build -j N_BUILD_PROCESSES
   cmake --install build

If you wish to create a debug build, simply pass the following arguments to CMake:

.. code-block:: bash

   cmake --build build -j N_BUILD_PROCESSES --config Debug

.. _testing:

Testing SHARPlib
----------------------
For unit tests, we make use of the `doctest single header source library <https://github.com/doctest/doctest>`_ found in the ``tests`` directory. In order to build and run the tests, execute the following commands from the project root directory:

.. code-block:: bash

   cmake -B build .
   cmake --build build -j N_BUILD_PROCESSES --target SHARPlib_tests
   ctest --test-dir build

