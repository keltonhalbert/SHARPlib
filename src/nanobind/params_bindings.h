#ifndef SHARPLIB_PARAMS_BINDINGS_H
#define SHARPLIB_PARAMS_BINDINGS_H

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>

// clang-format on
#include <SHARPlib/layer.h>
#include <SHARPlib/params/convective.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/winds.h>

#include "SHARPlib/params/fire.h"
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
            if ((pressure.shape(0) != height.shape(0)) ||
                (pressure.shape(0) != temperature.shape(0)) ||
                (pressure.shape(0) != dewpoint.shape(0)) ||
                (pressure.shape(0) != virtemp.shape(0))) {
                throw nb::buffer_error(
                    "All input arrays must have the same size!");
            }
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

References 
----------
Thompson et al. 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf

Parameters 
----------
lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
pressure : numpy.ndarray[dtype=float32]
    A 1D NumPy array of pressure values (Pa)
height : numpy.ndarray[dtype=float32] 
    A 1D NumPy array of height values (Pa)
temperature : numpy.ndarray[dtype=float32] 
    A 1D NumPy array of temperature values (K)
dewpoint : numpy.ndarray[dtype=float32] 
    A 1D NumPy array of dewpoint values (K)
virtemp : numpy.ndarray[dtype=float32] 
    A 1D NumPy array of virtual temperature values (K)
cape_thresh : float, default = 100.0
    The CAPE threshold used to compute the Effective Inflow Layer 
cinh_thresh : float, default = -250.0 
    The CINH threshold used to compute the Effective Inflow Layer
muplc : None or nwsspc.sharp.calc.parcel.Parcel, optional

Returns 
-------
nwsspc.sharp.calc.layer.PressureLayer
    The Effective Inflow Layer
    )pbdoc");

    m_params.def(
        "effective_inflow_layer",
        [](sharp::lifter_cm1& lifter, const_prof_arr_t pressure,
           const_prof_arr_t height, const_prof_arr_t temperature,
           const_prof_arr_t dewpoint, const_prof_arr_t virtemp,
           float cape_thresh, float cinh_thresh, sharp::Parcel* mupcl) {
            if ((pressure.shape(0) != height.shape(0)) ||
                (pressure.shape(0) != temperature.shape(0)) ||
                (pressure.shape(0) != dewpoint.shape(0)) ||
                (pressure.shape(0) != virtemp.shape(0))) {
                throw nb::buffer_error(
                    "All input arrays must have the same size!");
            }
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

References 
----------
Thompson et al. 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf

Parameters 
----------
lifter : nwsspc.sharp.calc.parcel.lifter_cm1 
pressure : numpy.ndarray[dtype=float32]
    A 1D NumPy array of pressure values (Pa)
height : numpy.ndarray[dtype=float32] 
    A 1D NumPy array of height values (Pa)
temperature : numpy.ndarray[dtype=float32] 
    A 1D NumPy array of temperature values (K)
dewpoint : numpy.ndarray[dtype=float32] 
    A 1D NumPy array of dewpoint values (K)
virtemp : numpy.ndarray[dtype=float32] 
    A 1D NumPy array of virtual temperature values (K)
cape_thresh : float, default = 100.0
    The CAPE threshold used to compute the Effective Inflow Layer 
cinh_thresh : float, default = -250.0 
    The CINH threshold used to compute the Effective Inflow Layer
muplc : None or nwsspc.sharp.calc.parcel.Parcel, optional

Returns 
-------
nwsspc.sharp.calc.layer.PressureLayer
    The Effective Inflow Layer
    )pbdoc");

    m_params.def(
        "storm_motion_bunkers",
        [](const_prof_arr_t pressure, const_prof_arr_t height,
           const_prof_arr_t u_wind, const_prof_arr_t v_wind,
           sharp::HeightLayer mean_wind_layer_agl,
           sharp::HeightLayer wind_shear_layer_agl, const bool leftMover,
           const bool pressureWeighted) {
            if ((pressure.shape(0) != height.shape(0)) ||
                (pressure.shape(0) != u_wind.shape(0)) ||
                (pressure.shape(0) != v_wind.shape(0))) {
                throw nb::buffer_error(
                    "All input arrays must have the same size!");
            }
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

References 
----------

Buners et al. 2000: https://doi.org/10.1175/1520-0434(2000)015%3C0061:PSMUAN%3E2.0.CO;2

Parameters 
----------
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of pressure values (Pa)
height : numpy.ndarray[dtype=float32] 
    1D NumPy array of height values (meters)
u_wind : numpy.ndarray[dtype=float32]
    1D NumPy array of U wind component values (m/s)
v_wind : numpy.ndarray[dtype=float32]
    1D NumPy array of V wind compnent values (m/s)
mean_wind_layer_agl : nwsspc.sharp.calc.layer.HeightLayer 
    HeightLayer (AGL) for computing the mean wind 
wind_shear_layer_agl : nwsspc.sharp.calc.layer.HeightLayer 
    HeightLayer (AGL) for computing wind_shear
leftMover : bool 
    Whether to compute left mover supercell motion (default: False)
pressureWeighted : bool 
    Whether to use the pressure weighted mean wind (default: False)

Returns
-------
nwsspc.sharp.calc.winds.WindComponents
    U, V wind components of storm motion (m/s)
    )pbdoc");

    m_params.def(
        "storm_motion_bunkers",
        [](const_prof_arr_t pressure, const_prof_arr_t height,
           const_prof_arr_t u_wind, const_prof_arr_t v_wind,
           sharp::PressureLayer eff_infl_lyr, sharp::Parcel& mupcl,
           const bool leftMover) {
            if ((pressure.shape(0) != height.shape(0)) ||
                (pressure.shape(0) != u_wind.shape(0)) ||
                (pressure.shape(0) != v_wind.shape(0))) {
                throw nb::buffer_error(
                    "All input arrays must have the same size!");
            }
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

References
----------
Bunkers et al. 2014: http://dx.doi.org/10.15191/nwajom.2014.0211

Parameters 
----------
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure values (Pa)
height : numpy.ndarray[dtype=float32]
    1D NumPy array of height values (meters)
u_wind : numpy.ndarray[dtype=float32] 
    1D NumPy array of U wind component values (m/s)
v_wind : numpy.ndarray[dtype=float32] 
    1D NumPy array of V wind component values (m/s)
eff_infl_lyr : nwsspc.sharp.calc.layer.PressureLayer 
    Effective Inflow Layer PressureLayer
mupcl : nwsspc.sharp.calc.parcel.Parcel 
    Most Unstable Parcel 
leftMover : bool 
    Whether or not to compute left moving supercell motion (default: False)

Returns
-------
nwsspc.sharp.calc.winds.WindComponents
    U, V wind components of storm motion (m/s)

    )pbdoc");

    m_params.def(
        "mcs_motion_corfidi",
        [](const_prof_arr_t pressure, const_prof_arr_t height,
           const_prof_arr_t u_wind, const_prof_arr_t v_wind)
            -> std::pair<sharp::WindComponents, sharp::WindComponents> {
            auto pres = pressure.view();
            auto hght = height.view();
            auto uwin = u_wind.view();
            auto vwin = v_wind.view();

            if ((pres.shape(0) != hght.shape(0)) ||
                (pres.shape(0) != uwin.shape(0)) ||
                (pres.shape(0) != vwin.shape(0))) {
                throw nb::buffer_error(
                    "pressure, height, u_wind, and v_wind must have the "
                    "same sizes!");
            }

            return sharp::mcs_motion_corfidi(pres.data(), hght.data(),
                                             uwin.data(), vwin.data(),
                                             pressure.size());
        },
        nb::arg("pressure"), nb::arg("height"), nb::arg("u_wind"),
        nb::arg("v_wind"), R"pbdoc(
Compute the Corfidi upshear and downshear MCS motion vectors.

Estimates the mesoscale convective system (MCS) motion vectors for upshear 
and downshear propagating convective systems as in Corfidi et al. 2003.
The method is based on observations that MCS motion is a function of 
1) the advection of existing cells by the mean wind and 
2) the propagation of new convection relative to existing storms.

References
----------
Corfidi et al. 2003: https://www.spc.noaa.gov/publications/corfidi/mcs2003.pdf

Parameters 
----------
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure values (Pa)
height : numpy.ndarray[dtype=float32]
    1D NumPy array of height values (meters)
u_wind : numpy.ndarray[dtype=float32]
    1D NumPy array of u-wind components (m/s)
v_wind : numpy.ndarray[dtype=float32]
    1D NumPy array of v-wind components (m/s)

Returns 
-------
tuple[nwsspc.sharp.calc.winds.WindComponents, nwsspc.sharp.calc.winds.WindComponents]
    (upshear, downshear)
    )pbdoc");

    m_params.def(
        "effective_bulk_wind_difference",
        [](const_prof_arr_t pres, const_prof_arr_t hght, const_prof_arr_t uwin,
           const_prof_arr_t vwin, sharp::PressureLayer eil,
           const float eql_pres) {
            if ((pres.shape(0) != hght.shape(0)) ||
                (pres.shape(0) != uwin.shape(0)) ||
                (pres.shape(0) != vwin.shape(0))) {
                throw nb::buffer_error(
                    "pressure, height, u_wind, and v_wind must have the "
                    "same sizes!");
            }
            return sharp::effective_bulk_wind_difference(
                pres.data(), hght.data(), uwin.data(), vwin.data(),
                pres.shape(0), eil, eql_pres);
        },
        nb::arg("pressure"), nb::arg("height"), nb::arg("u_wind"),
        nb::arg("v_wind"), nb::arg("effective_inflow_layer"),
        nb::arg("equilibrium_level_pressure"),
        R"pbdoc(
Compute the Effective Bulk Wind Difference 

The effective bulk wind difference is the wind shear between 
the bottom height of the effective inflow layer, and 50% of 
the equilibrium level depth. This is analogous to the usage 
of 0-6 km wind shear, but allows more flexibility for elevated 
convection. Returns MISSING if the effective inflow layer or 
equilibrium level pressure are MISSING.

    )pbdoc");

    m_params.def("energy_helicity_index", &sharp::energy_helicity_index,
                 nb::arg("cape"), nb::arg("helicity"),
                 R"pbdoc(
Computes the Energy Helicity Index.

EHI is a composite parameter based on the premise that 
storm rotation shoudl be maximized when CAPE is large 
and SRH is large. Typically, the layers used for helicity 
are either 0-1 km AGL or 0-3 km AGL.
 
References 
----------
https://doi.org/10.1175/1520-0434(2003)18%3C530:RSATFP%3E2.0.CO;2

Parameters 
----------
CAPE : float 
    Convective Available Potential Energy (J/kg)
helicity : float 
    Storm Relative Helicity (m^2 / s^2 a.k.a J/kg)

Returns
-------
float
    Energy Helicity Index (umitless)
    )pbdoc");

    m_params.def("significant_tornado_parameter",
                 &sharp::significant_tornado_parameter, nb::arg("pcl"),
                 nb::arg("lcl_hght_agl"), nb::arg("storm_relative_helicity"),
                 nb::arg("bulk_wind_difference"),
                 R"pbdoc(
Computes the Significant Tornado Parameter.


The Significant Tornado Parameter is used to diagnose environments
where tornadoes are favored. STP traditionally comes in two flavors:
fixed-layer, and effective-layer. Fixed-layer STP expects surface-based
CAPE, the surface-based LCL, 0-1 km storm-relative helicity, the
0-6 km bulk wind difference, and the surface-based CINH. For the
effective inflow layer based STP, use 100mb mixed-layer CAPE,
100mb mixed-layer LCL height AGL, effective-layer srh, the
effective layer bulk wind difference, and the 100mb mixed-layer
CINH. NOTE: The effective bulk wind difference is the shear between
the bottom of the effective inflow layer and 50% of the height of the
equilibrium level of the most unstable parcel.

References
----------
Thompson et al 2012: https://www.spc.noaa.gov/publications/thompson/waf-env.pdf

Parameters 
----------
pcl : nwsspc.sharp.calc.parcel.Parcel 
    For effective-layer STP, a mixed-layer parcel, and for fixed-layer STP, a surface-based parcel 
lcl_hght_agl : float 
    The parcel LCL height in meters
storm_relative_helicity : float 
    For effective-layer STP, effective SRH, and for fixed-layer, 0-1 km SRH (m^2 / s^2)
bulk_wind_difference : float 
    For effective-layer STP, effective BWD, and for fixed-layer STP, 0-6 km BWD (m/s)

Returns
-------
float 
    The Significant Tornado Parameter
    )pbdoc");

    m_params.def("supercell_composite_parameter",
                 &sharp::supercell_composite_parameter, nb::arg("mu_cape"),
                 nb::arg("eff_srh"), nb::arg("eff_shear"),
                 R"pbdoc(
Computes the Supercell Composite Parameter. 

The supercell composite parameter is used to diagnose environments
where supercells are favored. Requires computing most unstable
CAPE, effective layer storm relative helicity, and effective
bulk shear. Effective bulk shear is the vector difference between
the winds at the bottom of the effective inflow layer, and 50% of
the equilibrium level height. It is similar to the 0-6 km shear
vector, but allows for elevated supercell thunderstorms.

The left-moving supercell composite parameter can be computed by
providing effective SRH calculated using the bunkers left-moving
storm motion, and will return negative values.

References
----------
Thompson et al 2003: https://www.spc.noaa.gov/publications/thompson/ruc_waf.pdf

Thompson et al 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf

Thompson et al 2012: https://www.spc.noaa.gov/publications/thompson/waf-env.pdf

        

Parameters 
----------
mu_cape : float 
    The CAPE of the Most Unstable Parcel (J/kg)
eff_srh : float 
    Effective inflow layer Storm Relative Helicity (m^2/s^2) 
eff_shear : float 
    Effective layer shear (m/s)

Returns
-------
float
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

Parameters 
----------
mu_pcl : nwsspc.sharp.calc.parcel.Parcel 
    A precomputed Most Unstable parcel
lapse_rate_700_500mb : float
    The 700-500 mb lapse rate (K/km)
tmpk_500mb : float 
    The 500mb temperature (K)
freezing_level_agl : float 
    The height of the freezing level (AGL, meters)
shear_0_6km : float
    The 0-6 km shear vector magnitude (m/s)

Returns
-------
float
    The significant hail parameter
    )pbdoc");

    m_params.def("derecho_composite_parameter",
                 &sharp::derecho_composite_parameter, nb::arg("dcape"),
                 nb::arg("mucape"), nb::arg("shear_0_6km"),
                 nb::arg("mean_wind_0_6km"),
                 R"pbdoc(
The Derecho Composite Parameter (DCP) is based on a dataset of 113 derecho 
events compiled by the Evans and Boswell (2001) study. It is intended to 
hightlight environments favorable for cold pool/outflow driven convective 
events. The physical mechanisms behind this parameter focus on cold pool 
production (DCAPE), ability to sustain strong convection (MUCAPE), 
convective organization (0 - 6 km shear), and sufficient deep-layer flow 
within the environment.

References
----------
Evans and Doswell 2001: https://doi.org/10.1175/1520-0434(2001)016%3C0329:EODEUP%3E2.0.CO;2

Parameters
----------
dcape : float 
    Downdraft Convective Available Potential Energy (J/kg)
mucape : float 
    Most Unstable Parcel Convective Available Potential Energy (J/kg)
shear_0_6km : float 
    Shear magnitude in the 0 - 6 km layer AGL (m/s)
mean_wind_0_6km : float 
    Mean wind magnitude in the 0 - 6 km layer AGL (m/s)

Returns
-------
float 
    The Derecho Composite Parameter
)pbdoc");

    m_params.def(
        "large_hail_parameter",
        [](const sharp::Parcel& mu_pcl, const float lapse_rate_700_500mb,
           sharp::PressureLayer hail_growth_zone,
           const sharp::WindComponents storm_motion, const_prof_arr_t pres,
           const_prof_arr_t hght, const_prof_arr_t uwin,
           const_prof_arr_t vwin) {
            if ((pres.shape(0) != hght.shape(0)) ||
                (pres.shape(0) != uwin.shape(0)) ||
                (pres.shape(0) != vwin.shape(0))) {
                throw nb::buffer_error(
                    "pressure, height, u_wind, and v_wind must have the "
                    "same sizes!");
            }

            return sharp::large_hail_parameter(
                mu_pcl, lapse_rate_700_500mb, hail_growth_zone, storm_motion,
                pres.data(), hght.data(), uwin.data(), vwin.data(),
                pres.size());
        },
        nb::arg("mu_pcl"), nb::arg("lapse_rate_700_500mb"),
        nb::arg("hail_growth_zone"), nb::arg("storm_motion"),
        nb::arg("pressure"), nb::arg("height"), nb::arg("u_wind"),
        nb::arg("v_wind"),
        R"pbdoc(
Computes the Large Hail Parameter (LHP). The LHP is a multi-ingredient, 
composite index that includes thermodynamics and kinematics to attempt 
to detect environments that support very large hail. LHP has shown skill 
when differentiationg environments that support hail >= 3.5 in from those 
with < 2.0 in.

References
----------
Johnson and Sugden 2014: https://ejssm.org/archives/wp-content/uploads/2021/09/vol9-5.pdf

Parameters 
----------
mu_pcl : nwsspc.sharp.calc.parcel.Parcel 
    A previously computed most unstable parcel with CAPE and an equilibrium level
lapse_rate_700_500mb : float 
    700 - 500 hPa Lapse Rate (K)
hail_growth_zone : nwsspc.sharp.calc.layer.PressureLayer
    The layer of the atmosphere encompasing the Hail Growth Zone (Pa)
storm_motion : nwsspc.sharp.calc.winds.WindComponents 
    The storm motion vector to be used (m/s)
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure values (Pa)
height : numpy.ndarray[dtype=float32]
    1D NumPy array of height values (m)
u_wind : numpy.ndarray[dtype=float32]
    1D NumPy array of u_wind values (m/s)
v_wind : numpy.ndarray[dtype=float32]
    1D NumPy array of v_wind values (m/s)
    )pbdoc");

    m_params.def(
        "precipitable_water",
        [](sharp::PressureLayer& layer, const_prof_arr_t pres,
           const_prof_arr_t mixr) {
            if ((pres.shape(0) != mixr.shape(0))) {
                throw nb::buffer_error(
                    "All input arrays must have the same size!");
            }
            return sharp::precipitable_water(layer, pres.data(), mixr.data(),
                                             pres.size());
        },

        nb::arg("layer"), nb::arg("pres"), nb::arg("mixr"),
        R"pbdoc(
Given a PressureLayer to integrate over, compute the precipitable water 
from the given pressure and mixing ratio arrays.

Parameters 
----------
layer : nwsspc.sharp.calc.layer.PressureLayer 
    a PressureLayer over which to integrate (Pa)
pres : numpy.ndarray[dtype=float32] 
    1D NumPy array of presssure values (Pa)
mixr : numpy.ndarray[dtype=float32] 
    1D NumPy array of water vapor mixing ratio values (unitless)

Returns
-------
float
    Precipitable water content (mm)
    )pbdoc");

    m_params.def(
        "hail_growth_layer",
        [](const_prof_arr_t pres, const_prof_arr_t tmpk) {
            if ((pres.shape(0) != tmpk.shape(0))) {
                throw nb::buffer_error(
                    "All input arrays must have the same size!");
            }
            return sharp::hail_growth_layer(pres.data(), tmpk.data(),
                                            pres.size());
        },
        nb::arg("pressure"), nb::arg("temperature"),
        R"pbdoc(
Search for and return the PressureLayer of the lowest altitude 
hail growth zone. If none is found, the top and bottom of the 
PressureLayer are set to MISSING. 

Parameters 
----------
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure values (Pa)
temperature : numpy.ndarray[dtype=float32]
    1D NumPy array of temperature values (K)

Returns 
nwsspc.sharp.calc.layer.PressureLayer
    The PressureLayer containing the hail growth zone
    )pbdoc");

    m_params.def(
        "dendritic_layer",
        [](const_prof_arr_t pres, const_prof_arr_t tmpk) {
            if ((pres.shape(0) != tmpk.shape(0))) {
                throw nb::buffer_error(
                    "All input arrays must have the same size!");
            }
            return sharp::dendritic_layer(pres.data(), tmpk.data(),
                                          pres.size());
        },
        nb::arg("pressure"), nb::arg("temperature"),
        R"pbdoc(
Search for and return the PressureLayer of the lowest altitude 
dendritic growth zone. If none is found, the top and bottom of the 
PressureLayer are set to MISSING. 

Parameters 
----------
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of pressure values (Pa)
temperature : numpy.ndarray[dtype=float32]
    1D NumPy array of temperature values (K)

Returns 
-------
nwsspc.sharp.calc.layer.PressureLayer
    The PressureLayer containing the dendritic growth zone
    )pbdoc");

    m_params.def("snow_squall_parameter", &sharp::snow_squall_parameter,
                 nb::arg("wetbulb_2m"), nb::arg("mean_relh_0_2km"),
                 nb::arg("delta_thetae_0_2km"), nb::arg("mean_wind_0_2km"),
                 R"pbdoc(
The Snow Squall Parameter is a non-dimensional parameter that combines 
several ingredients believed to be beneficial for identifying snow squall 
environments by identifying the overlap of low-level potential instability, 
sufficient moisture, and strong low-level winds.


References
----------
Banacos et al. 2014: https://www.weather.gov/media/btv/research/Snow%20Squalls%20Forecasting%20and%20Hazard%20Mitigation.pdf
                
Parameters 
----------
wetbulb_2m : float 
    The surface wetbulb temperature, used to mask the parameter (K)
mean_relh_0_2km : float 
    The mean relative humidity between the surface and 2 km AGL (fraction)
delta_thetae_0_2km : float 
    The difference in equivalent potential temperature between 2 km AGL and the surface (K)
mean_wind_0_2km : float 
    The mean wind speed between the surface and 2 km AGL (m/s)

Returns 
-------
float 
    The Snow Squall Parameter
    )pbdoc");

    m_params.def("equilibrium_moisture_content",
                 &sharp::equilibrium_moisture_content, nb::arg("temperature"),
                 nb::arg("rel_humidity"),
                 R"pbdoc(
Compute the equilibrium moisture content for fuel 
as in Simard (1968).

Parameters
----------
temperature : float 
    The air temperature (K)
rel_humidity : float 
    Relative Humidity (fraction)

Returns 
-------
float 
    Equilibrium Moisture Content (fraction)
    )pbdoc");

    m_params.def(
        "equilibrium_moisture_content",
        [](const_prof_arr_t temperature, const_prof_arr_t rel_humidity) {
            auto tmpk = temperature.view();
            auto relh = rel_humidity.view();
            if ((tmpk.shape(0) != relh.shape(0))) {
                throw nb::buffer_error(
                    "temperature and rel_humidity must have the same size!");
            }

            float* emc_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                emc_arr[k] =
                    sharp::equilibrium_moisture_content(tmpk(k), relh(k));
            }

            nb::capsule owner(emc_arr,
                              [](void* p) noexcept { delete[] (float*)p; });
            return out_arr_t(emc_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("temperature"), nb::arg("rel_humidity"),
        R"pbdoc(
Compute the equilibrium moisture content for fuel 
as in Simard (1968).

Parameters
----------
temperature : numpy.ndarray[dtype=float32] 
    The air temperature (K)
rel_humidity : numpy.ndarray[dtype=float32]
    Relative Humidity (fraction)

Returns 
-------
numpy.ndarray[dtype=float32]
    Equilibrium Moisture Content (fraction)
    )pbdoc");

    m_params.def("fosberg_fire_index", &sharp::fosberg_fire_index,
                 nb::arg("temperature"), nb::arg("rel_humidity"),
                 nb::arg("wind_speed"),
                 R"pbdoc(
Compute the Fosberg Fire-Weather Index (FWWI) as in Fosberg (1978).

Parameters 
----------
temperature : float 
    The air temperature (K)
rel_humidity : float 
    The relative humidity (fraction)
wind_speed : float 
    Wind speed (m/s)

Returns 
-------
float 
    Fosberg Fire-Weather Index
    )pbdoc");

    m_params.def(
        "fosberg_fire_index",
        [](const_prof_arr_t temperature, const_prof_arr_t rel_humidity,
           const_prof_arr_t wind_speed) {
            auto tmpk = temperature.view();
            auto relh = rel_humidity.view();
            auto wspd = wind_speed.view();
            if ((tmpk.shape(0) != relh.shape(0)) ||
                (tmpk.shape(0) != wspd.shape(0))) {
                throw nb::buffer_error(
                    "temperature, rel_humidity, and wind_speed must have the "
                    "same sizes!");
            }

            float* fwwi = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                fwwi[k] = sharp::fosberg_fire_index(tmpk(k), relh(k), wspd(k));
            }

            nb::capsule owner(fwwi,
                              [](void* p) noexcept { delete[] (float*)p; });
            return out_arr_t(fwwi, {tmpk.shape(0)}, owner);
        },
        nb::arg("temperature"), nb::arg("rel_humidity"), nb::arg("wind_speed"),
        R"pbdoc(
Compute the Fosberg Fire-Weather Index (FWWI) as in Fosberg (1978).

Parameters 
----------
temperature : numpy.ndarray[dtype=float32]
    The air temperature (K)
rel_humidity : numpy.ndarray[dtype=float32]
    The relative humidity (fraction)
wind_speed : numpy.ndarray[dtype=float32]
    Wind speed (m/s)

Returns 
-------
numpy.ndarray[dtype=float32]
    Fosberg Fire-Weather Index
    )pbdoc");
}

#endif
