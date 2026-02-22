winds
=====

Vectors and Components
----------------------

.. doxygenstruct:: sharp::WindVector
   :members:

.. doxygenstruct:: sharp::WindComponents
   :members:

.. doxygenfunction:: sharp::components_to_vector(WindComponents)
.. doxygenfunction:: sharp::components_to_vector(float, float)
.. doxygenfunction:: sharp::vector_to_components(WindVector)
.. doxygenfunction:: sharp::vector_to_components(float, float)
.. doxygenfunction:: sharp::vector_magnitude
.. doxygenfunction:: sharp::vector_magnitude_precise
.. doxygenfunction:: sharp::vector_angle
.. doxygenfunction:: sharp::u_component
.. doxygenfunction:: sharp::v_component

Kinematic Variables
-------------------

.. doxygenfunction:: sharp::helicity
.. doxygenfunction:: sharp::wind_shear
.. doxygenfunction:: sharp::mean_wind
.. doxygenfunction:: sharp::max_wind(PressureLayer lyr, const float pressure[], const float uwin[], const float vwin[], const std::ptrdiff_t N, std::size_t* lvl_max)
.. doxygenfunction:: sharp::max_wind(HeightLayer lyr, const float pressure[], const float uwin[], const float vwin[], const std::ptrdiff_t N, std::size_t* lvl_max)
