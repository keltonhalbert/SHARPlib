Build and Install SHARPlib (Python)
===================================

.. _dependencies:

Dependencies 
------------
SHARPlib depends on ``numpy >= 1.23.0`` and supports ``python>=3.10``

.. _installing_python:

Installing SHARPlib 
----------------------------
SHARPlib is available via pip/PyPI, and can be installed for Linux, MacOS, and Windows 64-bit systems using:

.. code-block:: bash

   pip install sharplib


SHARPlib is also available via conda-forge, and can be installed with ``mamba`` or ``conda``:

.. code-block:: bash

   conda install sharplib


If you do not have conda-forge set as your default channel:

.. code-block:: bash

   conda install -c conda-forge sharplib

Building SHARPlib 
--------------------------
SHARPlib C++ code can be called from Python using `nanobind <https://github.com/wjakob/nanobind>`_ to handle the wrapping. Building SHARPlib with its Python bindings is quite easy-- you can simply clone this repository and install it via pip from the current directory:

.. code-block:: bash

   # clone from https 
   https://github.com/keltonhalbert/SHARPlib.git

   # or, from ssh
   git@github.com:keltonhalbert/SHARPlib.git

   # Ensure the submodule dependencies are cloned
   git submodule update --init --recursive
   pip install .

If you desire to manually build the SHARPlib library + python bindings, you may execute:

.. code-block:: bash

   cmake -B build . -DBUILD_PYBIND=ON
   cmake --build build -j N_BUILD_PROCESSES

.. NOTE::
   You will need a C++17 compatible compiler installed (i.e. cl, gcc, clang).

Testing SHARPlib 
-------------------------
Python tests are facilitated using ``pytest``. The ``pytest`` package must be installed in your environment to execute tests...

.. code-block:: bash

   # via pip 
   pip install pytest 

   # via conda 
   conda install -c conda-forge pytest

To execute the tests:

.. code-block:: bash

   pytest

in the root of this git repository to execute tests.
