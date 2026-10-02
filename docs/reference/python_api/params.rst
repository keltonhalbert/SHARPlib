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

   It runs in three steps. :func:`~nwsspc.sharp.calc.params.spectral_bin_cloud_top` finds the cloud top from the dewpoint depression and the relative humidity. A pre-classifier decides columns at or below 0 C throughout, columns above 0 C throughout, and columns with a cloud top warmer than Tice over a surface below 0 C. The bin microphysics integrates the other columns, level by level from the cloud top down. :func:`~nwsspc.sharp.calc.params.spectral_bin_classifier` without a cloud top runs all three steps. Given a cloud top, it runs the last two.

   SHARPlib ports the 2023 version of the algorithm. It follows the Python reference by D. Tripp (sbc_alg_2023Aug31.py and run_sbc.py, 2023-08-31), which the authors consider authoritative, and the C++ MRMS code by A. Rosenow and D. Tripp (versions 2.0.0 to 2.0.3), whose comments credit H. Reeves with the pre-classifier. The authors gave permission for the port. The work is NOAA-funded, and the port is in the public domain. Its documentation notes each place where it departs from the paper or from the reference.

   .. rubric:: Departures from Reeves et al. (2016)

   The port follows the 2023 code, which departs from the paper in these steps. Line references like classify:404 point into sbc_alg_2023Aug31.py.

   .. list-table::
      :header-rows: 1
      :widths: 16 38 46

      * - Step
        - Paper
        - 2023 code
      * - Cloud top
        - If column-max RH > 80 %, the level of highest RH; otherwise RA or SN from the surface Tw
        - First level from the top with T - Td <= 6 K and RH > 60 %. If below it max(T - Td) > 10 K or min RH < 40 %, search again below the driest level. Fallback: first level below the top with RH >= 80 %
      * - Column entirely at or below 0 C
        - SN if Tw at cloud top <= Tice, else FZRA
        - SN if Tw at cloud top < Tice **and** min Tw from the level nearest 3 km AGL down to the surface < Tice, else FZRA ("non-classical freezing rain")
      * - Cloud top warmer than Tice, surface below 0 C
        - integrate; refreezing can give PL
        - always FZRA (pre-classifier)
      * - Column entirely above 0 C
        - integrate
        - RA
      * - Surface decision
        - Fig. 2 flowchart: if Pi = 0 and Pw > 0, RA over a warm surface and FZRA over a cold one, at any crossing count. Otherwise, over a cold surface: PL if Pw/Pi < 0.15, FZRA if Pi/Pw < 0.15, else an FZRA-PL mix; over a warm surface: PL if Nc > 1, else RASN. (The text in §2d(3) instead says a warm surface gives RA if Pi/Pw < 0.15, otherwise PL.)
        - liquid-fraction thresholds of 0.85, 0.60, and 0.15 by crossing count and surface Tw. Adds RAPL. SN is possible over a warm surface
      * - Tice
        - -6 C
        - -6 C. After a subfreezing level whose previous level has every bin fully melted, Tice becomes -10 C for the rest of the column
      * - Riming
        - f_rim = 1
        - 1. After the first level handled by the subfreezing branch G (classify:404), and not after the frozen cloud-top branches A and B, melting layers switch to f_rim = 5 (graupel) and ice-pellet fall speeds (classify:264-266, 300-301)
      * - DSD
        - DSD25 with 0.1 mm bins (18 bins)
        - 4 bins

   .. rubric:: Quirks of the reference

   To match the reference, the port keeps its quirks: one refreeze level shared by all bins, quantities that read back as 0 where the reference does not set them, an ice density of 0.917 g cm^-3 for the surface flux and 0.918 for the class of each bin, the air density from the wet-bulb temperature, and a term of the refreezing rate that is always 0. The microphysics sections below describe each one, with the others it keeps.

   .. rubric:: Inputs

   * The profiles run from the surface up, and index 0 must be the surface (2 m) level. The reference's user guide requires the 2 m level before the cloud-top search. Height must strictly increase, AGL or MSL.
   * Pressure (Pa), height (m), temperature, dewpoint, and wet-bulb temperature (K), and relative humidity over liquid water (fraction). Relative humidity is an array, like the model relative humidity that the reference reads, rather than a value computed from temperature and dewpoint. That is why the 80 % fallback of the cloud-top rule can fire.
   * The caller supplies the wet-bulb temperature. The reference computes it with its own approximation (``wetbulb_calc`` in run_sbc.py), and its user guide says that the classifier performs best with that formulation. The golden data that SHARPlib is tested against use it too.

   .. rubric:: Outputs

   :class:`~nwsspc.sharp.calc.params.SpectralBinResult` holds the category, the liquid fraction of the precipitation mass reaching the surface, and the height of the lowest supercooled liquid water (m AGL). The liquid fraction is not a probability. No cloud, invalid inputs, or too few valid levels give a missing result. :class:`~nwsspc.sharp.calc.params.PrecipType` has stable integer values for gridded output, and 0 is left unused for a possible "no precipitation" category:

   .. list-table::
      :header-rows: 1
      :widths: 40 12 48

      * - ``params.PrecipType``
        - Value
        - Category
      * - ``missing``
        - -9999
        - missing, the value of ``constants.MISSING``
      * - ``rain``
        - 1
        - rain (RA)
      * - ``snow``
        - 2
        - snow (SN)
      * - ``rain_snow``
        - 3
        - rain and snow (RASN)
      * - ``freezing_rain``
        - 4
        - freezing rain (FZRA)
      * - ``ice_pellets``
        - 5
        - ice pellets (PL)
      * - ``freezing_rain_ice_pellets``
        - 6
        - freezing rain and ice pellets (FZRAPL)
      * - ``rain_ice_pellets``
        - 7
        - rain and ice pellets (RAPL)

   With ``return_profile=True``, the classifier also returns the liquid fraction of every bin at every level, mainly for validation.

   .. rubric:: Drop-size distribution, Tice, and riming

   :func:`~nwsspc.sharp.calc.params.spectral_bin_dsd_default` is the 4-bin distribution of the Python reference. To use another, such as the bins of the C++ MRMS code or a finer spacing, build it once with :func:`~nwsspc.sharp.calc.params.spectral_bin_dsd` from the diameters (mm) and concentrations, and pass it to every call.

   ``ice_nucleation_temperature`` is Tice, -6 C by default. It decides whether the cloud top is frozen, enters two rules of the pre-classifier, and is the temperature at or below which every bin refreezes. The -10 C that replaces it below a level where every bin has melted is fixed, as in the reference. The riming factor, from 1 (none) to 5 (graupel), belongs to the distribution: :func:`~nwsspc.sharp.calc.params.spectral_bin_dsd` takes it as ``rime_factor``, with the reference's value of 1 as the default.

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

   .. autofunction:: nwsspc.sharp.calc.params.spectral_bin_cloud_top

   .. Precipitation type from a given cloud top: pre-classifier

   .. autofunction:: nwsspc.sharp.calc.params.spectral_bin_classifier

   .. Microphysics: frozen cloud tops and melting

   .. rubric:: Microphysics: frozen cloud tops and melting

   A column that the pre-classifier does not decide runs the microphysics of the reference (``classify`` in sbc_alg_2023Aug31.py) for every bin of the drop-size distribution, level by level from the cloud top down to the surface. A layer spans two adjacent valid levels. The function first counts the 0 C crossings Nc between adjacent levels, where a level at exactly 0 C counts as subfreezing, and finds the first crossing.

   With Tw at or below Tice at the cloud top, every bin starts as snow, with the diameter, density, and fall-speed coefficients aa and bb of the drop-size distribution. Above the first crossing the snow falls unchanged. Its fall speed is the raindrop fall speed near the ground times sqrt(rho_0 / rho) / aa, with the air density rho and rho_0 = 1.292e-3 g cm^-3. At each level with Tw at or above 0 C, a bin that is not all liquid melts. Its liquid fraction grows with the heat flux from the air across the layer above, which depends on Tw and the relative humidity. A bin melts with the aspect ratio of a raindrop when its class at the level above is liquid-like, ice pellets included. Otherwise its aspect ratio is 0.8. The class of a bin is rain, snow, or both, from the liquid share of its mass flux and a threshold of 0.15.

   At the surface, the liquid and ice mass fluxes are Pw = sum(m0 fw v N) and Pi = sum(m0 (1 - fw) v N) / 0.917 over the bins, with the mass m0, liquid fraction fw, fall speed v, and concentration N of each bin. ``liquid_fraction`` is Pw / (Pw + Pi). The reference rounds it to 0.1 %, and this function does not. A surface with Tw above 0 C is warm:

   ===== ======= ================================================
   Nc    Surface Category by liquid_fraction
   ===== ======= ================================================
   1     warm    RA above 0.85, SN below 0.60, otherwise RASN
   > 1   warm    PL below 0.15, RA above 0.85, otherwise RAPL
   any   cold    PL below 0.15, FZRA above 0.85, otherwise FZRAPL
   ===== ======= ================================================

   FZRA and FZRAPL set ``supercooled_liquid_height`` to 0 m. This decision departs from Fig. 2 of Reeves et al. (2016), which compares Pw and Pi with a ratio of 0.15, has no RAPL, and gives RA, RASN, or PL, never SN, over a warm surface.

   The function keeps these quirks of the reference:

   * The air density p / (R_d Tw), with R_d = 287 J kg^-1 K^-1, uses the wet-bulb temperature.
   * aa and bb keep their cloud-top values all the way down. The reference recomputes them at melting levels but never reads the new values.
   * The surface ice flux uses an ice density of 0.917 g cm^-3, and the class of each bin 0.918.
   * A quantity that the reference does not set at a level reads back as 0 at the next level, because its arrays start at 0. The snow of a frozen cloud top has a liquid fraction of 0, and a bin that has melted completely has no snow mass and no ice or snow volume.
   * The reference weights the mass flux of each bin by the bin width and by the ratio of its fall speed at the surface to its fall speed at the level. That ratio is always 1 where it is used, and both factors cancel in every ratio of fluxes, so the function leaves them out.

   With ``return_profile=True``, the profile holds the liquid fraction of every bin at every level of the integration.

   .. Microphysics: refreezing

   .. rubric:: Microphysics: refreezing

   Every other level with Tw below 0 C refreezes. With a frozen cloud top, these are the levels below the first crossing. There a bin refreezes when it holds ice at the level above, with fw below 1, and every bin refreezes when Tw is at or below Tice. A refreezing bin loses liquid mass as it gives heat to the air across the layer above, by conduction and by vapor exchange. The rate depends on Tw and the relative humidity and follows Kumjian et al. (2012), as in Reeves et al. (2016). The fall speed of the bin is v_PL + (v_r - v_PL) fw / fw_r. Here fw is its liquid fraction at the level above, v_r and fw_r are its fall speed and liquid fraction at the refreeze level, and v_PL is the ice pellet fall speed near the ground times sqrt(rho_0 / rho). Its class comes from the liquid share of its mass flux and a threshold of 0.15. It is freezing rain FZRA, ice pellets PL, or the mix FZRAPL, with FZDZ and FZDZPL for drizzle, the bins below 0.6 mm. A later melting level melts all of these classes, PL included, with the aspect ratio of a raindrop.

   A bin that does not refreeze falls unchanged. Rain becomes FZRA or FZDZ, and rain mixed with snow or ice pellets becomes FZRAPL or FZDZPL. ``supercooled_liquid_height`` takes the height of each level where a bin refreezes into a class with liquid, or where one of these conversions happens, so it ends at the lowest one. A bin that falls unchanged with any other class, FZRA included, leaves it alone.

   Two rules of the reference depart from Reeves et al. (2016):

   * The Tice switch. At a level of this kind where every bin is all liquid at the level above, Tice becomes 263.15 K, or -10 C, for the rest of the column. The new value decides the refreezing test and the refreezing rate at this level, and the branch at each level below.
   * The graupel switch. Below the first level of this kind, a melting level uses a riming factor of 5 instead of the one of the drop-size distribution, and the ice pellet fall speed near the ground times sqrt(rho_0 / rho).

   The function keeps these quirks of the reference:

   * One refreeze level serves all bins. Whenever a bin starts refreezing, at its first refreezing level since the cloud top or the last melting level, the refreeze level moves to the level above. The move applies to every bin, including bins that started refreezing earlier. At the level of the move, it reaches the bins after the mover in the drop-size distribution but not those before it.
   * The denominator of the refreezing rate has a term xsi D / D_w that is always 0, because the reference reads the diameter D of the level before it sets it.
   * A refreeze level of 0, the cloud top, means that nothing has refrozen. A melting level classes ice as ice pellets once something has refrozen, and as snow before. Refreezing that starts just below the cloud top sets the refreeze level to 0, so its ice still counts as snow.
   * The absolute humidity at ice saturation uses the dry-air density p / (R_d T), from the temperature rather than Tw.
   * A refreezing bin has no liquid, ice, or snow volume at the next level, and a bin that falls unchanged has no snow-core diameter there.

   .. Microphysics: liquid cloud tops

   .. rubric:: Microphysics: liquid cloud tops

   With Tw above Tice at the cloud top, every bin starts as a liquid drop of its melted diameter D. It falls at the raindrop fall speed near the ground times sqrt(rho_0 / rho). A drop below 0.6 mm is freezing drizzle and a larger one freezing rain. The cloud top becomes the supercooled-liquid height, even when it is warmer than 0 C, as in the reference.

   Above the first crossing, the drops fall unchanged at the fall speed of Foote and du Toit (1969), -0.193 + 4.96 D - 0.904 D^2 + 0.0566 D^3, times exp(z / 20 km):

   * Below a cloud top at or below 0 C, at each level warmer than Tice, unless some level above the first crossing warms through Tice. The drops stay supercooled, and each of these levels lowers the supercooled-liquid height.
   * Below a cloud top warmer than 0 C, at each level warmer than Tice. The drops are rain.

   Other levels above the first crossing follow the rules for the levels below it. Both conditions compare Tw with the Tice in use, which refreezing can switch to -10 C. Whether a level warms through Tice uses the Tice passed in, as in the reference.

   z is the height as passed, AGL or MSL. The Python reference passes height AGL and the C++ MRMS code height MSL. The fall speed changes the result only when some drops refreeze below. Moving z by 1500 m changed no result in tests of the reference.

   Reeves et al. (2016) only state that every drop of a cloud top warmer than Tice starts as liquid.

   .. Precipitation type from a full sounding

   :func:`~nwsspc.sharp.calc.params.spectral_bin_classifier` also runs from a full sounding, without a cloud top. Both overloads are documented together above.
