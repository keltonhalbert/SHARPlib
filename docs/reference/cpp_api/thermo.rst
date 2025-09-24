thermo
======

Dry Variables
-------------

.. doxygenfunction:: sharp::theta
.. doxygenfunction:: sharp::theta_level 
.. doxygenfunction:: sharp::lapse_rate(HeightLayer, const float[], const float[], const std::ptrdiff_t)
.. doxygenfunction:: sharp::lapse_rate(PressureLayer, const float[], const float[], const float[], const std::ptrdiff_t)
.. doxygenfunction:: sharp::lapse_rate_max(HeightLayer, const float, const float[], const float[], const std::ptrdiff_t, HeightLayer*)
.. doxygenfunction:: sharp::lapse_rate_max(PressureLayer, const float, const float[], const float[], const float[], const std::ptrdiff_t, PressureLayer*)

Dry Ascent
----------

.. doxygenfunction:: sharp::drylift 
.. doxygenfunction:: sharp::lcl_temperature 

Moist Variables 
---------------
.. doxygenfunction:: sharp::mixratio(float)
.. doxygenfunction:: sharp::mixratio(float, float)
.. doxygenfunction:: sharp::mixratio_ice
.. doxygenfunction:: sharp::vapor_pressure 
.. doxygenfunction:: sharp::vapor_pressure_ice
.. doxygenfunction:: sharp::specific_humidity 
.. doxygenfunction:: sharp::relative_humidity
.. doxygenfunction:: sharp::temperature_at_mixratio
.. doxygenfunction:: sharp::virtual_temperature
.. doxygenfunction:: sharp::density
.. doxygenfunction:: sharp::buoyancy(const float, const float)
.. doxygenfunction:: sharp::buoyancy(const float[], const float[], float[], std::ptrdiff_t)
.. doxygenfunction:: sharp::thetae
.. doxygenfunction:: sharp::wetbulb(Lft, float, float, float)
.. doxygenfunction:: sharp::wetbulb(lifter_wobus, float, float, float)
.. doxygenfunction:: sharp::wetbulb(lifter_cm1, float, float, float)
.. doxygenfunction:: sharp::theta_wetbulb(Lft, float, float, float)
.. doxygenfunction:: sharp::theta_wetbulb(lifter_wobus, float, float, float)
.. doxygenfunction:: sharp::theta_wetbulb(lifter_cm1, float, float, float)

Moist Ascent
------------

.. doxygenenum:: sharp::adiabat
.. doxygenfunction:: sharp::wobf
.. doxygenfunction:: sharp::wetlift 
.. doxygenfunction:: sharp::saturated_lift 
.. doxygenfunction:: sharp::moist_adiabat_cm1


