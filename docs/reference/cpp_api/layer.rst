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

QC builds, the default, treat a level whose data is ``sharp::MISSING`` or NaN
as missing. ``layer_min`` and ``layer_max`` skip missing levels, and
interpolate the layer bottom and top across them as ``interp_height`` and
``interp_pressure`` do. A layer with no valid data returns ``MISSING``. So
does a layer that lies wholly outside the profile, and the level it reports is
the layer's endpoint nearest the profile. Builds configured with
``-DNO_QC=ON`` skip these checks.

.. doxygenfunction:: sharp::layer_min
.. doxygenfunction:: sharp::layer_max
.. doxygenfunction:: sharp::layer_mean(PressureLayer, const float[], const float[], const std::ptrdiff_t)
.. doxygenfunction:: sharp::layer_mean(HeightLayer, const float[], const float[], const float[], const std::ptrdiff_t, const bool)
.. doxygenfunction:: sharp::integrate_layer_trapz
