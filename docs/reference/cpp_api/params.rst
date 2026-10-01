params
======

Convective
----------

.. doxygenfunction:: sharp::effective_inflow_layer(Lifter&, const float[], const float[], const float[], const float[], const float[], float[], float[], const std::ptrdiff_t, const float, const float, Parcel*)
.. doxygenfunction:: sharp::storm_motion_bunkers(const float[], const float[], const float[], const float[], const std::ptrdiff_t, HeightLayer, HeightLayer, const bool, const bool)
.. doxygenfunction:: sharp::storm_motion_bunkers(const float[], const float[], const float[], const float[], const std::ptrdiff_t, PressureLayer, const Parcel&, const bool)
.. doxygenfunction:: sharp::mcs_motion_corfidi
.. doxygenfunction:: sharp::effective_bulk_wind_difference
.. doxygenfunction:: sharp::energy_helicity_index
.. doxygenfunction:: sharp::supercell_composite_parameter
.. doxygenfunction:: sharp::significant_tornado_parameter
.. doxygenfunction:: sharp::significant_hail_parameter
.. doxygenfunction:: sharp::derecho_composite_parameter
.. doxygenfunction:: sharp::large_hail_parameter
.. doxygenfunction:: sharp::hail_growth_layer
.. doxygenfunction:: sharp::convective_temperature(Lifter&, const float[], const float[], const float[], const float[], const float[], float[], float[], const std::ptrdiff_t, float)
.. doxygenfunction:: sharp::precipitable_water

Fire
----

.. doxygenfunction:: sharp::equilibrium_moisture_content
.. doxygenfunction:: sharp::fosberg_fire_index
.. doxygenfunction:: sharp::pft_plume_potential_temperature
.. doxygenfunction:: sharp::pft_plume_mixratio
.. doxygenfunction:: sharp::pyrocumulonimbus_firepower_threshold

Winter
------

.. doxygenfunction:: sharp::dendritic_layer
.. doxygenfunction:: sharp::snow_squall_parameter

Precipitation Type
~~~~~~~~~~~~~~~~~~

Precipitation-type probabilities from the modified Bourgouin method
(Birk et al. 2021, https://doi.org/10.1175/WAF-D-20-0118.1).

.. doxygenvariable:: sharp::BOURGOUIN_PRESSURE_MIN

.. doxygenstruct:: sharp::BourgouinEnergy
   :members:

.. doxygenstruct:: sharp::PrecipTypeProbabilities
   :members:

.. Wet-bulb melting and refreezing energies from a sounding

.. doxygenfunction:: sharp::bourgouin_energy

.. Precipitation generation layer from a sounding

.. doxygenfunction:: sharp::precipitation_generation_layer

.. Probability of ice, and precipitation-type probabilities from energies

.. doxygenfunction:: sharp::probability_of_ice
.. doxygenfunction:: sharp::modified_bourgouin(const BourgouinEnergy&, const float, const float)

.. Precipitation type from a full sounding

.. doxygenfunction:: sharp::modified_bourgouin(const float[], const float[], const float[], const float[], const float[], const std::ptrdiff_t, const float, const float, const float)

