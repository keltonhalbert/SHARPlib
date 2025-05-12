#ifndef SHARPLIB_PARAMS_BINDINGS_H
#define SHARPLIB_PARAMS_BINDINGS_H

// clang-format off
#include <nanobind/nanobind.h>

// clang-format on
#include <SHARPlib/layer.h>
#include <SHARPlib/params/convective.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/winds.h>

#include "SHARPlib/params/winter.h"
#include "sharplib_types.h"

namespace nb = nanobind;

inline void make_params_bindings(nb::module_ m) {
    nb::module_ m_params =
        m.def_submodule("params",
                        "Sounding and Hodograph Analysis and Research Program "
                        "Library (SHARPlib) :: Derived Parameters");

    m_params.def(
        "effective_inflow_layer",
        [](sharp::lifter_wobus& lifter, const_prof_arr_t pressure,
           const_prof_arr_t height, const_prof_arr_t temperature,
           const_prof_arr_t dewpoint, const_prof_arr_t virtemp,
           float cape_thresh, float cinh_thresh, sharp::Parcel* mupcl) {
            float* pcl_vtmp = new float[height.size()];
            float* pcl_buoy = new float[height.size()];

            sharp::PressureLayer eil = sharp::effective_inflow_layer(
                lifter, pressure.data(), height.data(), temperature.data(),
                dewpoint.data(), virtemp.data(), pcl_vtmp, pcl_buoy,
                height.size(), cape_thresh, cinh_thresh, mupcl);

            delete[] pcl_vtmp;
            delete[] pcl_buoy;

            return eil;
        },
        nb::arg("lifter"), nb::arg("pressure"), nb::arg("height"),
        nb::arg("temperature"), nb::arg("dewpoint"), nb::arg("virtemp"),
        nb::arg("cape_thresh") = 100.0f, nb::arg("cinh_thresh") = -250.0f,
        nb::arg("mupcl") = nb::none(),
        R"pbdoc(
Computes the Effective Inflow Layer, or the layer of the atmosphere
beliefed to be the primary source of inflow for supercell thunderstorms. 
The Effective Inflow Layer, and its use in computing shear and storm 
relative helicity, is described by Thompson et al. 2007:
    https://www.spc.noaa.gov/publications/thompson/effective.pdf

Standard/default values for cape_thresh and cinh_thresh have been 
experimentally determined to be cape_thresh = 100 J/kg and 
cinh_thresh = -250.0 J/kg. If an empty parcel object is passed via the 
'mupcl' kwarg, the Most Unstable parcel found during the EIL search will 
be returned. 

Parameters:
    lifter: An instantiated lifter_cm1() or lifter_wobus()
    pressure: A 1D NumPy array of pressure values (Pa)
    height: A 1D NumPy array of height values (Pa)
    temperature: A 1D NumPy array of temperature values (K)
    dewpoint: A 1D NumPy array of dewpoint values (K)
    virtemp: A 1D NumPy array of virtual temperature values (K)


    )pbdoc");

    m_params.def(
        "effective_inflow_layer",
        [](sharp::lifter_cm1& lifter, const_prof_arr_t pressure,
           const_prof_arr_t height, const_prof_arr_t temperature,
           const_prof_arr_t dewpoint, const_prof_arr_t virtemp,
           float cape_thresh, float cinh_thresh, sharp::Parcel* mupcl) {
            float* pcl_vtmp = new float[height.size()];
            float* pcl_buoy = new float[height.size()];

            sharp::PressureLayer eil = sharp::effective_inflow_layer(
                lifter, pressure.data(), height.data(), temperature.data(),
                dewpoint.data(), virtemp.data(), pcl_vtmp, pcl_buoy,
                height.size(), cape_thresh, cinh_thresh, mupcl);

            delete[] pcl_vtmp;
            delete[] pcl_buoy;

            return eil;
        },
        nb::arg("lifter"), nb::arg("pressure"), nb::arg("height"),
        nb::arg("temperature"), nb::arg("dewpoint"), nb::arg("virtemp"),
        nb::arg("cape_thresh") = 100.0f, nb::arg("cinh_thresh") = -250.0f,
        nb::arg("mupcl") = nb::none(),
        R"pbdoc(
Computes the Effective Inflow Layer, or the layer of the atmosphere
beliefed to be the primary source of inflow for supercell thunderstorms. 
The Effective Inflow Layer, and its use in computing shear and storm 
relative helicity, is described by Thompson et al. 2007:
    https://www.spc.noaa.gov/publications/thompson/effective.pdf

Standard/default values for cape_thresh and cinh_thresh have been 
experimentally determined to be cape_thresh = 100 J/kg and 
cinh_thresh = -250.0 J/kg. If an empty parcel object is passed via the 
'mupcl' kwarg, the Most Unstable parcel found during the EIL search will 
be returned. 

Parameters:
    lifter: An instantiated lifter_cm1() or lifter_wobus()
    pressure: A 1D NumPy array of pressure values (Pa)
    height: A 1D NumPy array of height values (Pa)
    temperature: A 1D NumPy array of temperature values (K)
    dewpoint: A 1D NumPy array of dewpoint values (K)
    virtemp: A 1D NumPy array of virtual temperature values (K)


    )pbdoc");

    m_params.def(
        "storm_motion_bunkers",
        [](const_prof_arr_t pressure, const_prof_arr_t height,
           const_prof_arr_t u_wind, const_prof_arr_t v_wind,
           sharp::HeightLayer mean_wind_layer_agl,
           sharp::HeightLayer wind_shear_layer_agl, const bool leftMover,
           const bool pressureWeighted) {
            sharp::WindComponents storm_mtn = sharp::storm_motion_bunkers(
                pressure.data(), height.data(), u_wind.data(), v_wind.data(),
                height.size(), mean_wind_layer_agl, wind_shear_layer_agl,
                leftMover, pressureWeighted);

            return storm_mtn;
        },
        nb::arg("pressure"), nb::arg("height"), nb::arg("u_wind"),
        nb::arg("v_wind"), nb::arg("mean_wind_layer_agl"),
        nb::arg("wind_shear_layer_agl"), nb::arg("leftMover") = false,
        nb::arg("pressureWeighted") = false,
        R"pbdoc(
Estimates the supercell storm motion using the Bunkers et al. 2000 method 
described in the following paper:
    https://doi.org/10.1175/1520-0434(2000)015%3C0061:PSMUAN%3E2.0.CO;2
        
This does not use any of the updated methods described by Bunkers et al. 2014, 
which uses Effective Inflow Layer metrics to get better estimates of storm 
motion, especially when considering elevated convection. 

Params:
    pressure: 1D NumPy array of pressure values (Pa)
    height: 1D NumPy array of height values (meters)
    u_wind: 1D NumPy array of U wind component values (m/s)
    v_wind: 1D NumPy array of V wind compnent values (m/s)
    mean_wind_layer_agl: HeightLayer (AGL) for computing the mean wind 
    wind_shear_layer_agl: HeightLayer (AGL) for computing wind_shear
    leftMover: Whether to compute left mover supercell motion (default: False)
    pressureWeighted: Whether to use the pressure weighted mean wind (default: False)

Returns:
    U, V wind components of storm motion (m/s)
    )pbdoc");

    m_params.def(
        "storm_motion_bunkers",
        [](const_prof_arr_t pressure, const_prof_arr_t height,
           const_prof_arr_t u_wind, const_prof_arr_t v_wind,
           sharp::PressureLayer eff_infl_lyr, sharp::Parcel& mupcl,
           const bool leftMover) {
            sharp::WindComponents storm_mtn = sharp::storm_motion_bunkers(
                pressure.data(), height.data(), u_wind.data(), v_wind.data(),
                height.size(), eff_infl_lyr, mupcl, leftMover);

            return storm_mtn;
        },
        nb::arg("pressure"), nb::arg("height"), nb::arg("u_wind"),
        nb::arg("v_wind"), nb::arg("eff_infl_lyr"), nb::arg("mupcl"),
        nb::arg("leftMover") = false,
        R"pbdoc(
Estimates supercell storm motion using the Bunkers et al. 2014 
method described in the following paper:
    http://dx.doi.org/10.15191/nwajom.2014.0211
    
This method is parcel based, using a mean-wind vector defined as the 
pressure-weighted mean wind between the Effective Inflow Layer surface 
(see effective_inflow_layer routine) and 65% of the depth between that 
surface and the most unstable parcel's Equilibrium Level. This method 
produces the same storm motion estimate for surface based supercells, 
and captures the motion of elevated supercells better than the 
Bunkers 2000 method. 

The input parameters of eff_infl_lyr and mupcl (effective inflow layer 
pressure bounds and the most unstable parcel, respectively) are required
to be precomputed and passed to this routine. These are expensive 
operations that are presumed to be computed at some other point 
in the analysis pipeline. 

Parameters:
    pressure: 1D NumPy array of pressure values (Pa)
    height: 1D NumPy array of height values (meters)
    u_wind: 1D NumPy array of U wind component values (m/s)
    v_wind: 1D NumPy array of V wind component values (m/s)
    eff_infl_lyr: Effective Inflow Layer PressureLayer
    mupcl: Most Unstable Parcel 
    leftMover: Whether or not to compute left moving supercell motion (default: False)

Returns:
    U, V wind components of storm motion (m/s)

    )pbdoc");

    m_params.def("energy_helicity_index", &sharp::energy_helicity_index,
                 nb::arg("cape"), nb::arg("helicity"),
                 R"pbdoc(
Computes the Energy Helicity Index.

Parameters:
    cape: Convective Available Potential Energy (J/kg)
    helicity: Storm Relative Helicity (m^2 / s^2 a.k.a J/kg)

Returns:
    Energy Helicity Index (umitless)
    )pbdoc");

    m_params.def("significant_tornado_parameter",
                 &sharp::significant_tornado_parameter, nb::arg("pcl"),
                 nb::arg("lcl_hght_agl"), nb::arg("storm_relative_helicity"),
                 nb::arg("bulk_wind_difference"),
                 R"pbdoc(
Computes the Significant Tornado Parameter.

References:
    Thompson et al 2012: https://www.spc.noaa.gov/publications/thompson/waf-env.pdf

Parameters:
    pcl: For effective-layer STP, a mixed-layer parcel, and for fixed-layer STP, a surface-based parcel 
    lcl_hght_agl: The parcel LCL height in meters
    storm_relative_helicity: For effecitve-layer STP, effecitve SRH, and for fixed-layer, 0-1 km SRH (m^2 / s^2)
    bulk_wind_difference: For effective-layer STP, effecitve BWD, and for fixed-layer STP, 0-6 km BWD (m/s)
    )pbdoc");

    m_params.def("supercell_composite_parameter",
                 &sharp::supercell_composite_parameter, nb::arg("mu_cape"),
                 nb::arg("eff_srh"), nb::arg("eff_shear"),
                 R"pbdoc(
Computes the Supercell Composite Parameter. 

References:
    Thompson et al 2003: https://www.spc.noaa.gov/publications/thompson/ruc_waf.pdf
    Thompson et al 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf
    Thompson et al 2012: https://www.spc.noaa.gov/publications/thompson/waf-env.pdf
        

Parameters:
    mu_cape: The CAPE of the Most Unstable Parcel (J/kg)
    eff_srh: Effecive inflow layer Storm Relative Helicity (m^2/s^2) 
    eff_shear: Effective layer shear (m/s)

Returns:
    Supercell Composite Parameter (unitless)
    )pbdoc");

    m_params.def("significant_hail_parameter",
                 &sharp::significant_hail_parameter, nb::arg("mu_pcl"),
                 nb::arg("lapse_rate_700_500mb"), nb::arg("tmpk_500mb"),
                 nb::arg("freezing_level_agl"), nb::arg("shear_0_6km"),
                 R"pbdoc(
Compute the significant hail parameter, given a precomputed most-unstable parcel, 
the 700-500 mb lapse rate, the 500mb temperature, the height (AGL) of the 
freezing level, and the 0-6 km shear magnitude.

The Sig. Hail Parameter (SHIP) was developed using a large database of 
surface-modified, observed severe hail proximity soundings. It is based on 
parameters, and is meant to delineate between SIG (>=2" diameter) and NON-SIG
(<2" diameter) hail environments.

SHIP = [(MUCAPE j/kg) * (Mixing Ratio of MU PARCEL g/kg) *  
        (700-500mb LAPSE RATE c/km) * (-500mb TEMP C) *
        (0-6km Shear m/s) ] / 42,000,000

0-6 km shear is confined to a range of 7-27 m s-1, mixing ratio is confined to 
a range of 11-13.6 g kg-1, and the 500 mb temperature is set to -5.5 C for 
any warmer values.

Once the initial version of SHIP is calculated, the values are modified in 
the following scenarios:

1) If MUCAPE < 1300 J kg-1, SHIP = SHIP * (MUCAPE/1300); 2) if 700-500 mb 
lapse rate < 5.8 C km-1, SHIP = SHIP * (lr75/5.8); 3) if freezing 
level < 2400 m AGL, SHIP = SHIP * (fzl/2400)

It is important to note that SHIP is NOT a forecast hail size.

Since SHIP is based on the RAP depiction of MUCAPE - unrepresentative MUCAPE 
"bullseyes" may cause a similar increase in SHIP values. This typically occurs 
when bad surface observations get into the RAP model.

Developed in the same vein as the STP and SCP parameters, values of SHIP 
greater than 1.00 indicate a favorable environment for SIG hail. Values greater 
than 4 are considered very high. In practice, maximum contour values of 1.5-2.0 
or higher will typically be present when SIG hail is going to be reported. 

Parameters:
    mu_pcl: A precomputed Most Unstable parcel
    lapse_rate_700_500mb: The 700-500 mb lapse rate (K/km)
    tmpk_500mb: The 500mb temperature (K)
    freezing_level_agl: The height of the freezing level (AGL, meters)
    shear_0_6km: The 0-6 km shear vector magnitude (m/s)

Returns:
    The significant hail parameter
    )pbdoc");

    m_params.def(
        "precipitable_water",
        [](sharp::PressureLayer& layer, const_prof_arr_t pres,
           const_prof_arr_t mixr) {
            return sharp::precipitable_water(layer, pres.data(), mixr.data(),
                                             pres.size());
        },

        nb::arg("layer"), nb::arg("pres"), nb::arg("mixr"),
        R"pbdoc(
Given a PressureLayer to integrate over, compute the precipitable water 
from the given pressure and mixing ratio arrays.

Parameters:
    layer: a PressureLayer over which to integrate (Pa)
    pres: 1D NumPy array of presssure values (Pa)
    mixr: 1D NumPy array of water vapor mixing ratio values (unitless)

Returns:
    Precipitable water content (mm)
    )pbdoc");

    m_params.def(
        "dendritic_layer",
        [](const_prof_arr_t pres, const_prof_arr_t tmpk) {
            return sharp::dendritic_layer(pres.data(), tmpk.data(),
                                          pres.size());
        },
        nb::arg("pressure"), nb::arg("temperature"),
        R"pbdoc(
Search for and return the PressureLayer of the lowest altitude 
dendritic growth zone. If none is found, the top and bottom of the 
PressureLayer are set to MISSING. 

Parameters:
    pressure: 1D NumPy array of pressure values (Pa)
    temperature: 1D NumPy array of temperature values (K)

Return:
    The PressureLayer containing the dendritic growth zone
    )pbdoc");
}

#endif
