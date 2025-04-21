#ifndef SHARPLIB_PARAMS_BINDINGS_H
#define SHARPLIB_PARAMS_BINDINGS_H

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

// clang-format on
#include <SHARPlib/layer.h>
#include <SHARPlib/params/convective.h>
#include <SHARPlib/parcel.h>

#include <tuple>

#include "SHARPlib/winds.h"
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

            return std::make_tuple(storm_mtn.u, storm_mtn.v);
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
           sharp::PressureLayer eff_infl_lyr, sharp::Parcel* mupcl,
           const bool leftMover) {
            sharp::WindComponents storm_mtn = sharp::storm_motion_bunkers(
                pressure.data(), height.data(), u_wind.data(), v_wind.data(),
                height.size(), eff_infl_lyr, mupcl, leftMover);

            return std::make_tuple(storm_mtn.u, storm_mtn.v);
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

    m_params.def(
        "significant_tornado_parameter",
        [](sharp::Parcel pcl, float lcl_hght_agl, float storm_relative_helicity,
           float bulk_wind_difference) {
            return sharp::significant_tornado_parameter(pcl, lcl_hght_agl,
                                                        storm_relative_helicity,
                                                        bulk_wind_difference);
        },
        nb::arg("pcl"), nb::arg("lcl_hght_agl"),
        nb::arg("storm_relative_helicity"), nb::arg("bulk_wind_difference"),
        R"pbdoc(
Computes the Significant Tornado Parameter
    )pbdoc");
}

#endif
