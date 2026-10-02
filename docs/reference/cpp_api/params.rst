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

It runs in three steps. ``sharp::spectral_bin_cloud_top`` finds the cloud top
from the dewpoint depression and the relative humidity. A pre-classifier
decides columns at or below 0 C throughout, columns above 0 C throughout, and
columns with a cloud top warmer than Tice over a surface below 0 C. The bin
microphysics integrates the other columns, level by level from the cloud top
down. ``sharp::spectral_bin_classifier`` without a cloud top runs all three
steps. Given a cloud top, it runs the last two.

SHARPlib ports the 2023 version of the algorithm. It follows the Python
reference by D. Tripp (sbc_alg_2023Aug31.py and run_sbc.py, 2023-08-31),
which the authors consider authoritative, and the C++ MRMS code by A. Rosenow
and D. Tripp (versions 2.0.0 to 2.0.3), whose comments credit H. Reeves with
the pre-classifier. The authors gave permission for the port. The work is
NOAA-funded, and the port is in the public domain. Its documentation notes
each place where it departs from the paper or from the reference.

.. rubric:: Departures from Reeves et al. (2016)

The port follows the 2023 code, which departs from the paper in these steps.
Line references like classify:404 point into sbc_alg_2023Aug31.py.

.. list-table::
   :header-rows: 1
   :widths: 16 38 46

   * - Step
     - Paper
     - 2023 code
   * - Cloud top
     - If column-max RH > 80 %, the level of highest RH; otherwise RA or SN
       from the surface Tw
     - First level from the top with T - Td <= 6 K and RH > 60 %. If below
       it max(T - Td) > 10 K or min RH < 40 %, search again below the driest
       level. Fallback: first level below the top with RH >= 80 %
   * - Column entirely at or below 0 C
     - SN if Tw at cloud top <= Tice, else FZRA
     - SN if Tw at cloud top < Tice **and** min Tw from the level nearest
       3 km AGL down to the surface < Tice, else FZRA ("non-classical
       freezing rain")
   * - Cloud top warmer than Tice, surface below 0 C
     - integrate; refreezing can give PL
     - always FZRA (pre-classifier)
   * - Column entirely above 0 C
     - integrate
     - RA
   * - Surface decision
     - Fig. 2 flowchart: if Pi = 0 and Pw > 0, RA over a warm surface and
       FZRA over a cold one, at any crossing count. Otherwise, over a cold
       surface: PL if Pw/Pi < 0.15, FZRA if Pi/Pw < 0.15, else an FZRA-PL
       mix; over a warm surface: PL if Nc > 1, else RASN. (The text in
       §2d(3) instead says a warm surface gives RA if Pi/Pw < 0.15,
       otherwise PL.)
     - liquid-fraction thresholds of 0.85, 0.60, and 0.15 by crossing count
       and surface Tw. Adds RAPL. SN is possible over a warm surface
   * - Tice
     - -6 C
     - -6 C. After a subfreezing level whose previous level has every bin
       fully melted, Tice becomes -10 C for the rest of the column
   * - Riming
     - f_rim = 1
     - 1. After the first level handled by the subfreezing branch G
       (classify:404), and not after the frozen cloud-top branches A and B,
       melting layers switch to f_rim = 5 (graupel) and ice-pellet fall
       speeds (classify:264-266, 300-301)
   * - DSD
     - DSD25 with 0.1 mm bins (18 bins)
     - 4 bins

.. rubric:: Quirks of the reference

To match the reference, the port keeps its quirks: one refreeze level shared
by all bins, quantities that read back as 0 where the reference does not set
them, an ice density of 0.917 g cm^-3 for the surface flux and 0.918 for the
class of each bin, the air density from the wet-bulb temperature, and a term
of the refreezing rate that is always 0. The microphysics sections of
``sharp::spectral_bin_classifier`` describe each one, with the others it
keeps.

.. rubric:: Inputs

* The profiles run from the surface up, and index 0 must be the surface
  (2 m) level. The reference's user guide requires the 2 m level before the
  cloud-top search. Height must strictly increase, AGL or MSL.
* Pressure (Pa), height (m), temperature, dewpoint, and wet-bulb temperature
  (K), and relative humidity over liquid water (fraction). Relative humidity
  is an array, like the model relative humidity that the reference reads,
  rather than a value computed from temperature and dewpoint. That is why
  the 80 % fallback of the cloud-top rule can fire.
* The caller supplies the wet-bulb temperature. The reference computes it
  with its own approximation (wetbulb_calc in run_sbc.py), and its user
  guide says that the classifier performs best with that formulation. The
  golden data that SHARPlib is tested against use it too.

.. rubric:: Outputs

``sharp::SpectralBinResult`` holds the category, the liquid fraction of the
precipitation mass reaching the surface, and the height of the lowest
supercooled liquid water (m AGL). The liquid fraction is not a probability.
No cloud, invalid inputs, or too few valid levels give a missing result.
``sharp::PrecipType`` has stable integer values for gridded output, and 0 is
left unused for a possible "no precipitation" category:

.. list-table::
   :header-rows: 1
   :widths: 40 12 48

   * - ``sharp::PrecipType``
     - Value
     - Category
   * - ``missing``
     - -9999
     - missing, the value of ``sharp::MISSING``
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

The optional ``liquid_fraction_profile`` receives the liquid fraction of
every bin at every level, mainly for validation.

.. rubric:: Drop-size distribution, Tice, and riming

``sharp::spectral_bin_dsd_default`` is the 4-bin distribution of the Python
reference. To use another, such as the bins of the C++ MRMS code or a finer
spacing, build it once with ``sharp::spectral_bin_dsd`` from the diameters
(mm) and concentrations of up to ``sharp::SBC_MAX_BINS`` bins, and pass it to
every call.

``ice_nucleation_temperature`` is Tice, -6 C by default
(``sharp::SBC_ICE_NUCLEATION_TEMPERATURE``). It decides whether the cloud top
is frozen, enters two rules of the pre-classifier, and is the temperature at
or below which every bin refreezes. The -10 C that replaces it below a level
where every bin has melted is fixed, as in the reference. The riming factor,
from 1 (none) to 5 (graupel), belongs to the distribution:
``sharp::spectral_bin_dsd`` takes it, with the reference's value of 1 as the
default.

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

.. doxygenfunction:: sharp::spectral_bin_cloud_top

.. Precipitation type from a given cloud top: pre-classifier

.. doxygenfunction:: sharp::spectral_bin_classifier(const float[], const float[], const float[], const float[], const float[], const float[], const std::ptrdiff_t, const float, const SpectralBinDSD&, const float, float[])

.. Microphysics: frozen cloud tops and melting

.. Microphysics: refreezing

.. Microphysics: liquid cloud tops

.. Precipitation type from a full sounding

.. doxygenfunction:: sharp::spectral_bin_classifier(const float[], const float[], const float[], const float[], const float[], const float[], const std::ptrdiff_t, const SpectralBinDSD&, const float, float[])
