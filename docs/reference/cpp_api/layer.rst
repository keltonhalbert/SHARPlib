layer
=====

Layer Definitions
-----------------

.. doxygenstruct:: sharp::PressureLayer
   :members:

.. doxygenstruct:: sharp::HeightLayer
   :members:

.. doxygenstruct:: sharp::LayerIndex 
   :members:

Layer Conversions
-----------------

.. doxygenfunction:: sharp::pressure_layer_to_height
.. doxygenfunction:: sharp::height_layer_to_pressure
.. doxygenfunction:: sharp::get_layer_index(PressureLayer&, const float[], const std::ptrdiff_t)
.. doxygenfunction:: sharp::get_layer_index(HeightLayer&, const float[], const std::ptrdiff_t)

Layer Calculations
------------------

.. doxygenfunction:: sharp::layer_min
.. doxygenfunction:: sharp::layer_max
.. doxygenfunction:: sharp::layer_mean(PressureLayer, const float[], const float[], const std::ptrdiff_t)
.. doxygenfunction:: sharp::layer_mean(HeightLayer, const float[], const float[], const float[], const std::ptrdiff_t, const bool)
.. doxygenfunction:: sharp::integrate_layer_trapz
