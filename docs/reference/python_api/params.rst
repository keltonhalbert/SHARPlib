params
======

.. automodule:: nwsspc.sharp.calc.params 

   Convective
   ----------
   Convective-based parameters, largely geared towards severe thunderstorms and tornadoes.

   .. autofunction:: nwsspc.sharp.calc.params.effective_inflow_layer
   .. autofunction:: nwsspc.sharp.calc.params.storm_motion_bunkers
   .. autofunction:: nwsspc.sharp.calc.params.mcs_motion_corfidi
   .. autofunction:: nwsspc.sharp.calc.params.effective_bulk_wind_difference
   .. autofunction:: nwsspc.sharp.calc.params.energy_helicity_index
   .. autofunction:: nwsspc.sharp.calc.params.significant_tornado_parameter
   .. autofunction:: nwsspc.sharp.calc.params.supercell_composite_parameter
   .. autofunction:: nwsspc.sharp.calc.params.significant_hail_parameter
   .. autofunction:: nwsspc.sharp.calc.params.derecho_composite_parameter
   .. autofunction:: nwsspc.sharp.calc.params.large_hail_parameter
   .. autofunction:: nwsspc.sharp.calc.params.hail_growth_layer
   .. autofunction:: nwsspc.sharp.calc.params.convective_temperature
   .. autofunction:: nwsspc.sharp.calc.params.precipitable_water

   Fire
   ---- 
   Fire-weather parameters

   .. autofunction:: nwsspc.sharp.calc.params.equilibrium_moisture_content
   .. autofunction:: nwsspc.sharp.calc.params.fosberg_fire_index
   .. autofunction:: nwsspc.sharp.calc.params.pyrocumulonimbus_firepower_threshold

   Winter
   ------
   Winter-weather parameters.

   .. autofunction:: nwsspc.sharp.calc.params.dendritic_layer
   .. autofunction:: nwsspc.sharp.calc.params.snow_squall_parameter

   Precipitation Type
   ~~~~~~~~~~~~~~~~~~
   Precipitation-type probabilities from the modified Bourgouin method (Birk et al. 2021, https://doi.org/10.1175/WAF-D-20-0118.1).

   The defaults follow the paper. For 1 Hz soundings and other noisy, high-resolution profiles, start with ``min_depth=100.0`` (m) and ``min_energy=2.0`` (J/kg). 2 J/kg is the melting-layer minimum of the original Bourgouin (2000) method. Both options are deviations from the paper. :func:`~nwsspc.sharp.calc.params.bourgouin_energy` and :func:`~nwsspc.sharp.calc.params.precipitation_generation_layer` describe when they help, the measurements behind these values, and what the options cost. Neither option changes the total melting energy.

   The melting and refreezing energies ignore levels above ``pressure_min``. Its default, 25000 Pa (250 hPa), is available as ``params.BOURGOUIN_PRESSURE_MIN``.

   .. autoclass:: nwsspc.sharp.calc.params.BourgouinEnergy
      :members: melting_energy_total, melting_energy_aloft, refreezing_energy

   .. autoclass:: nwsspc.sharp.calc.params.PrecipTypeProbabilities
      :members: rain, snow, freezing_rain, ice_pellets

   .. Wet-bulb melting and refreezing energies from a sounding

   .. autofunction:: nwsspc.sharp.calc.params.bourgouin_energy

   .. Precipitation generation layer from a sounding

   .. autofunction:: nwsspc.sharp.calc.params.precipitation_generation_layer

   .. Probability of ice, and precipitation-type probabilities from energies

   .. autofunction:: nwsspc.sharp.calc.params.probability_of_ice
   .. autofunction:: nwsspc.sharp.calc.params.modified_bourgouin

   .. Precipitation type from a full sounding

   :func:`~nwsspc.sharp.calc.params.modified_bourgouin` also computes the probabilities from the profiles of a full sounding. Both overloads are documented together above.

   Spectral bin classifier
   ^^^^^^^^^^^^^^^^^^^^^^^
   Precipitation type from the spectral bin classifier (SBC) of Reeves et al. (2016, https://doi.org/10.1175/JAMC-D-16-0044.1). The classifier follows a spectrum of drop sizes from the cloud top to the surface, computes the liquid fraction of each through melting and refreezing layers, and returns one of seven categories at the surface. It gives no probabilities.

   SHARPlib ports the 2023 version of the algorithm. It follows the Python reference by D. Tripp, which the authors consider authoritative, and the C++ MRMS code by A. Rosenow and D. Tripp. The authors gave permission for the port, and its documentation notes each place where it departs from the paper. Drop-size distribution diameters are in mm.

   .. Result types and the drop-size distribution

   .. autoclass:: nwsspc.sharp.calc.params.PrecipType
      :members:

   .. autoclass:: nwsspc.sharp.calc.params.SpectralBinResult
      :members: precip_type, liquid_fraction, supercooled_liquid_height

   A drop-size distribution holds at most ``params.SBC_MAX_BINS`` (64) bins. The default ice nucleation temperature, 267.15 K (-6 C), is available as ``params.SBC_ICE_NUCLEATION_TEMPERATURE``.

   .. autoclass:: nwsspc.sharp.calc.params.SpectralBinDSD
      :members: nbins, rime_factor, diameter, concentration

   .. autofunction:: nwsspc.sharp.calc.params.spectral_bin_dsd
   .. autofunction:: nwsspc.sharp.calc.params.spectral_bin_dsd_default

   .. Cloud top from a sounding

   .. Precipitation type from a given cloud top: pre-classifier

   .. Microphysics: frozen cloud tops and melting

   .. Microphysics: refreezing

   .. Microphysics: liquid cloud tops

   .. Precipitation type from a full sounding
