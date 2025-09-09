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

.. |docs_badge| image:: https://github.com/keltonhalbert/SHARPlib/actions/workflows/doxygen-gh-pages.yml/badge.svg
   :target: https://github.com/keltonhalbert/SHARPlib/actions/workflows/doxygen-gh-pages.yml
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
* Maintaining an interactive data visualization application complete with live access data feeds along with a core computational library was challenging.
* The computational library component wasn't very generalizable, in large part because I wrote the code while taking an OOP class, and everything became a nail to the object hammer.

SHARPlib seeks to take the successes of SHARPpy, while having the benefit of more experience and hindsight. It is most analogous to the ``sharppy.sharptab`` import, but written in C++ for performance and wrapped for Python/Numpy using `nanobind <https://github.com/wjakob/nanobind>`_. It is separate from any visualization software and dependencies, generalized to be more composable where appropriate, and optimised for performance.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   build_install_cpp.rst
   build_install_python.rst

.. toctree::
   :maxdepth: 2
   :caption: Reference:

   python_api/index.rst
   cpp_api/index.rst
