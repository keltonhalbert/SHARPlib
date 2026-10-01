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
