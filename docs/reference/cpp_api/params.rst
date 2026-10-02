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

The defaults follow the paper. For 1 Hz soundings and other noisy,
high-resolution profiles, start with ``min_depth = 100`` m and
``min_energy = 2`` J/kg. 2 J/kg is the melting-layer minimum of the original
Bourgouin (2000) method. Both options are deviations from the paper.
``sharp::bourgouin_energy`` and ``sharp::precipitation_generation_layer``
describe when they help, the measurements behind these values, and what the
options cost. Neither option changes the total melting energy.

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


Spectral bin classifier
^^^^^^^^^^^^^^^^^^^^^^^

Precipitation type from the spectral bin classifier (SBC) of Reeves et al.
(2016, https://doi.org/10.1175/JAMC-D-16-0044.1). The classifier follows a
spectrum of drop sizes from the cloud top to the surface, computes the liquid
fraction of each through melting and refreezing layers, and returns one of
seven categories at the surface. It gives no probabilities.

SHARPlib ports the 2023 version of the algorithm. It follows the Python
reference by D. Tripp, which the authors consider authoritative, and the C++
MRMS code by A. Rosenow and D. Tripp. The authors gave permission for the
port, and its documentation notes each place where it departs from the
paper. Drop-size distribution diameters are in mm.

.. Result types and the drop-size distribution

.. doxygenenum:: sharp::PrecipType

.. doxygenstruct:: sharp::SpectralBinResult
   :members:

.. doxygenvariable:: sharp::SBC_MAX_BINS
.. doxygenvariable:: sharp::SBC_ICE_NUCLEATION_TEMPERATURE

.. doxygenstruct:: sharp::SpectralBinDSD
   :members:

.. doxygenfunction:: sharp::spectral_bin_dsd
.. doxygenfunction:: sharp::spectral_bin_dsd_default

.. Cloud top from a sounding

.. Precipitation type from a given cloud top: pre-classifier

.. Microphysics: frozen cloud tops and melting

.. Microphysics: refreezing

.. Microphysics: liquid cloud tops

.. Precipitation type from a full sounding
