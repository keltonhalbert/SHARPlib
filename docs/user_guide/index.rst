User Guide
==========

This User Guide can be used to show you how to build or install precompiled binaries of the Python library on Linux, MacOS, and Windows, and how to build the source C++ library for usage in other compiled projects.   

SHARPlib is purely a mathematical computations library, and does not provide any means of visualizing data. If working with the Python library, it is recommended you use `matplotlib <https://matplotlib.org/>`_ or `MetPy <https://unidata.github.io/MetPy/latest/index.html>`_.

SHARPlib is intended to have minimal dependencies, and be buildable on RHEL8 or newer systems. Builds are continually tested on Linux, MacOS, and Windows to ensure consistent ease of building regardless of platform. The library also sticks to the C++17 standard to strike a balance between modern features and the widest possible system compatibility. The Python library supports the same versions of Python that NumPy supports, since NumPy is a core dependency for arrays. 

.. toctree::
   :maxdepth: 2
   :hidden:

   build_install_cpp.rst
   build_install_python.rst
