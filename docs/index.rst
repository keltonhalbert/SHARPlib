.. SHARPlib documentation master file

========
SHARPlib
========

.. |cpp_ci_badge| image:: https://github.com/keltonhalbert/SHARPlib/actions/workflows/cmake.yml/badge.svg
   :target: https://github.com/keltonhalbert/SHARPlib/actions/workflows/cmake.yml
   :alt: C++ CI (Linux, MacOS, Windows)

.. |python_ci_badge| image:: https://github.com/keltonhalbert/SHARPlib/actions/workflows/python.yml/badge.svg
   :target: https://github.com/keltonhalbert/SHARPlib/actions/workflows/python.yml
   :alt: Python CI (Linux, MacOS, Windows)

.. |wheels_badge| image:: https://github.com/keltonhalbert/SHARPlib/actions/workflows/wheels.yml/badge.svg
   :target: https://github.com/keltonhalbert/SHARPlib/actions/workflows/wheels.yml
   :alt: Build Wheels

.. |docs_badge| image:: https://github.com/keltonhalbert/SHARPlib/actions/workflows/docs.yml/badge.svg
   :target: https://github.com/keltonhalbert/SHARPlib/actions/workflows/docs.yml
   :alt: Build Docs

.. |conda_version_badge| image:: https://anaconda.org/conda-forge/sharplib/badges/version.svg
   :target: https://anaconda.org/conda-forge/sharplib
   :alt: Conda Forge Version

.. |pypi_version_badge| image:: https://img.shields.io/pypi/v/SHARPlib
   :target: https://pypi.org/project/SHARPlib/
   :alt: PyPI Version

.. |python_version_badge| image:: https://img.shields.io/pypi/pyversions/SHARPlib
   :alt: PyPI - Python Version

.. |conda_platforms_badge| image:: https://anaconda.org/conda-forge/sharplib/badges/platforms.svg
   :target: https://anaconda.org/conda-forge/sharplib
   :alt: Conda Forge Platforms

|cpp_ci_badge| |python_ci_badge| |wheels_badge| |docs_badge| |conda_version_badge| |pypi_version_badge| |python_version_badge| |conda_platforms_badge|

**Sounding and Hodograph Analysis and Research Program (SHARP)** C++ library for conducting analysis of atmospheric sounding profiles. Based on the NSHARP routines written by John Hart and Rich Thompson at the NWS Storm Prediction Center in Norman, Oklahoma.

.. _about_sharplib:

About SHARPlib
==============
Since the 1990s, the National Weather Service (NWS) Storm Prediction Center (SPC) has actively researched, developed, and maintained various software packages and libraries in-house for the post-processing and visualization of atmospheric sounding data. Initially developed by John Hart, SHARP was developed to display and compute derived atmospheric indices from observed weather balloons, and vertical profiles from model forecast data. Eventually it was included in N-AWIPS/GEMPAK/AWIPS2 as NSHARP, used to process archive research data as SHARPTab, used to generate web graphics as SHARPGIF, used to post-process model and mesoanalysis gridded data, ported to Python as SHARPpy, and is used currently in SPC operations as BigSHARP.

SHARPpy sought to open source these computations and visualizations in order to facilitate reproducible open science, provide a cross-platform means of visualizing sounding data, and internally simplify the number of SHARP derivatives with unique code bases. While many of the goals were met by SHARPpy, it struggled in a few key areas:

* Pure python is not a performant solution when post-processing gridded data, particularly when dealing with parcel lifting routines.
* Maintaining an interactive data visualization application, complete with live access data feeds, and with a core computational library, was challenging.
* The computational library component wasn't very generalizable, in large part because I wrote the code while taking an OOP class, and everything became a nail to the object hammer.

SHARPlib seeks to continue with and improve upon the successes of SHARPpy, while having the benefit of more experience and hindsight. It is most analogous to the ``sharppy.sharptab`` import, but written in C++ for performance and wrapped for Python/Numpy using `nanobind <https://github.com/wjakob/nanobind>`_. It is separate from any visualization software and dependencies, generalized to be more composable where appropriate, and optimised for performance.

License
=======
This software was developed by the US government:

National Oceanic and Atmospheric Administration

National Weather Service 

National Centers for Environmental Prediction 

Storm Prediction Center

https://www.spc.noaa.gov/

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>

.. toctree::
   :maxdepth: 2
   :hidden:

   user_guide/index
   reference/index
   faq

