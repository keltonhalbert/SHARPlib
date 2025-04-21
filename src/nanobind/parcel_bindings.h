#ifndef SHARPLIB_PARCEL_BINDINGS
#define SHARPLIB_PARCEL_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

// clang-format on
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>

#include "sharplib_types.h"

namespace nb = nanobind;

inline void make_parcel_bindings(nb::module_ m) {
    nb::module_ m_parcel = m.def_submodule(
        "parcel",
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Parcel Lifting Routines");

    // Bind the constructors, named fields, and default arguments
    nb::class_<sharp::lifter_wobus>(m_parcel, "lifter_wobus",
                                    R"pbdoc(
A functor that calls the Wobus Wetlift for computation of moist adiabats.

This is used to wrap the Wobus Wetlift function for parcel lifting 
routines. Functors -- classes with their operator() overloaded -- 
are used so that functions can be passed to templates in a way that 
the compiler can still optimize, rather than using function pointers
or lambdas. 

Specifically, this functor is designed to be passed as a template 
argument to Parcel::lift_parcel, so that the method of computing 
moist adiabats can be changed without changing the overall parcel 
lifting code. The reason this is awesome is that the compiler 
can still optimize and inline this code, while the user can 
configure the parcel lifting algorithm to their specifications. 

)pbdoc")
        .def(nb::init<>())
        .def_ro_static("lift_from_lcl", &sharp::lifter_wobus::lift_from_lcl,
                       R"pbdoc(
A static flag that helps the parcel lifting functions know where to lift from.
)pbdoc")
        .def("setup", &sharp::lifter_wobus::setup, nb::arg("lcl_pres"),
             nb::arg("lcl_tmpk"), R"pbdoc(
Some parcel lifters require setup in order to handle adiabatic ascent
when tracking the water vapor, liquid, and ice mixing ratios. The Wobus 
lifter does not require this, however, so this function does nothing.
             )pbdoc")
        .def("__call__", &sharp::lifter_wobus::operator(), nb::is_operator(),
             nb::arg("pres"), nb::arg("tmpk"), nb::arg("new_pres"),
             R"pbdoc(
Overloads the call operator in order to facilitate parcel lifting. 

Parameters:
    pres: Parcel pressure (Pa)
    tmpk: Parcel temperature (K)
    new_pres: Final level of parcel after lift (Pa)

Returns:
    The temperature of the lifted parcel (K)
             )pbdoc")
        .def("parcel_virtual_temperature",
             &sharp::lifter_wobus::parcel_virtual_temperature, nb::arg("pres"),
             nb::arg("tmpk"), R"pbdoc(
Computes the virtual temperature of the parcel (after saturation).

Parameters:
    pres: Parcel pressure (Pa)
    tmpk: Parcel temperature (K)

Returns:
    The virtual temperature of the parcel (K)
             )pbdoc");

    nb::class_<sharp::lifter_cm1>(m_parcel, "lifter_cm1", R"pbdoc(
Use the CM1 moist lift calculations for adiabatic and pseudoadiabatic
parcel ascent.
            )pbdoc")
        .def(nb::init<>())
        .def_ro_static("lift_from_lcl", &sharp::lifter_cm1::lift_from_lcl,
                       R"pbdoc(
A static flag that helps the parcel lifting functions know where to lift from.
The lifter_cm1 lifts from the last lifted level, rather than the LCL, because
it is an iterative solver. This results in major performance improvements while 
maintaining accuracy.
                   )pbdoc")
        .def_rw("ma_type", &sharp::lifter_cm1::ma_type, R"pbdoc(
The type of moist adiabat to use, defined by sharp::adiabat
            )pbdoc")
        .def_rw("pressure_incr", &sharp::lifter_cm1::pressure_incr, R"pbdoc(
The pressure increment (Pa) to use for the iterative solver. 
             )pbdoc")
        .def_rw("converge", &sharp::lifter_cm1::converge, R"pbdoc(
The iterative convergence criteria (K)
             )pbdoc")
        .def("setup", &sharp::lifter_cm1::setup, nb::arg("lcl_pres"),
             nb::arg("lcl_tmpk"), R"pbdoc(
This function sets the total water mixing ratio for 
adiabatic parcel ascent, and zeroes out the vapor,
liquid, and ice mixing ratio from previous parcel 
ascents. 

Parameters:
    lcl_pres: The LCL pressure (Pa)
    lcl_tmpk: The LCL temperature (K)
             )pbdoc")
        .def("__call__", &sharp::lifter_cm1::operator(), nb::is_operator(),
             nb::arg("pres"), nb::arg("tmpk"), nb::arg("new_pres"),
             R"pbdoc(
Lifts a parcel moist adiabatically/pseudoadiabatically using
sharp::moist_adiabat_cm1.

Parameters:
    pres: Parcel pressure (Pa)
    tmpk: Parcel temperature (K)
    new_pres: Final level of parcel after lift (Pa)

Returns:
    The temperature of the lifted parcel (K)
             )pbdoc")

        .def("parcel_virtual_temperature",
             &sharp::lifter_cm1::parcel_virtual_temperature, nb::arg("pres"),
             nb::arg("tmpk"), R"pbdoc(
Computes the virtual temperature of the parcel (after saturation).

Parameters:
    pres: Parcel pressure (Pa)
    tmpk: Parcel temperature (K)

Returns:
    The virtual temperature of the parcel (K)
             )pbdoc");

    nb::enum_<sharp::LPL>(m_parcel, "LPL")
        .value("SFC", sharp::LPL::SFC, "A Surface Based Parcel")
        .value("FCST", sharp::LPL::FCST, "A Forecast Surface Parcel")
        .value("MU", sharp::LPL::MU, "A Most Unstable Parcel")
        .value("ML", sharp::LPL::ML, "A Mixed-Layer Parcel")
        .value("USR", sharp::LPL::USR, "A User Defined Parcel");

    nb::class_<sharp::Parcel>(m_parcel, "Parcel", R"pbdoc(
Contains information about a Parcel's starting level and
thermodynamic attributes, as well as derived computations,
methods for constructing a parcel, and parcel ascent routines. 
                              )pbdoc")
        .def_rw("pres", &sharp::Parcel::pres, R"pbdoc(
Parcel starting pressure (Pa)
                )pbdoc")
        .def_rw("tmpk", &sharp::Parcel::tmpk, R"pbdoc(
Parcel starting temperature (K)
                )pbdoc")
        .def_rw("dwpk", &sharp::Parcel::dwpk, R"pbdoc(
Parcel starting dewpoint (K)
                )pbdoc")
        .def_rw("lcl_pressure", &sharp::Parcel::lcl_pressure, R"pbdoc(
Pressure at the Lifted Condensation Level (Pa)
                )pbdoc")
        .def_rw("lfc_pressure", &sharp::Parcel::lfc_pressure, R"pbdoc(
Pressure at the Level of Free Convection (Pa)
                )pbdoc")
        .def_rw("eql_pressure", &sharp::Parcel::eql_pressure, R"pbdoc(
Pressure at the parcel Equilibrium Level
                )pbdoc")
        .def_rw("cape", &sharp::Parcel::cape, R"pbdoc(
Parcel Convective Available Potential Energy (J/kg) between the LFC and EL
                )pbdoc")
        .def_rw("cinh", &sharp::Parcel::cinh, R"pbdoc(
Parcel Convective Inhibition (J/kg) between the LFC and EL
                )pbdoc")
        .def(nb::init<>())
        .def(
            nb::init<const float, const float, const float, const sharp::LPL>(),
            nb::arg("pressure"), nb::arg("temperature"), nb::arg("dewpoint"),
            nb::arg("lpl"), R"pbdoc(
Constructor for a Parcel

Parameters:
    pressure: Parcel initial pressure (Pa)
    temperature: Parcel initial temperature (K)
    dewpoint: Parcel initial dewpoint (K)
    lpl: Parcel Lifted Parcel Level (LPL) definition
        )pbdoc")
        .def(
            "lift_parcel",
            [](sharp::Parcel& pcl, sharp::lifter_wobus& lifter,
               const_prof_arr_t pressure) {
                const std::ptrdiff_t NZ = pressure.size();
                float* virtemp_arr = new float[NZ];

                pcl.lift_parcel(lifter, pressure.data(), virtemp_arr, NZ);

                nb::capsule owner(virtemp_arr,
                                  [](void* p) noexcept { delete[] (float*)p; });
                return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                    virtemp_arr, {pressure.shape(0)}, owner);
            },
            nb::arg("lifter"), nb::arg("pressure"),
            R"pbdoc(
Lifts a Parcel dry adiabatically from its LPL to its LCL dry
adiabatically, and then moist adiabatically from the LCL to 
the top of the profile. The moist adiabat used is determined
by the type of lifting functor passed to the function (i.e.
lifter_wobus or lifter_cm1).

Parameters:
    lifter: An instantiated lifter_wobus functor
    pressure: 1D NumPy array of Pressure levels for lifting (Pa)

Returns:
    A 1D NumPy array of parcel virtual temperature values (K)
        )pbdoc")
        .def(
            "lift_parcel",
            [](sharp::Parcel& pcl, sharp::lifter_cm1& lifter,
               const_prof_arr_t pressure) {
                const size_t NZ = pressure.size();
                float* virtemp_arr = new float[NZ];

                pcl.lift_parcel(lifter, pressure.data(), virtemp_arr, NZ);

                nb::capsule owner(virtemp_arr,
                                  [](void* p) noexcept { delete[] (float*)p; });
                return nb::ndarray<nb::numpy, float, nb::ndim<1>>(virtemp_arr,
                                                                  {NZ}, owner);
            },
            nb::arg("lifter"), nb::arg("pressure"),
            R"pbdoc(
Lifts a Parcel dry adiabatically from its LPL to its LCL dry
adiabatically, and then moist adiabatically from the LCL to 
the top of the profile. The moist adiabat used is determined
by the type of lifting functor passed to the function (i.e.
lifter_wobus or lifter_cm1).

Parameters:
    lifter: An instantiated lifter_cm1 functor
    pressure: 1D NumPy array of Pressure levels for lifting (Pa)

Returns:
    A 1D NumPy array of parcel virtual temperature values (K)
        )pbdoc")
        .def(
            "find_lfc_el",
            [](sharp::Parcel& pcl, const_prof_arr_t pres, const_prof_arr_t hght,
               const_prof_arr_t buoy) {
                pcl.find_lfc_el(pres.data(), hght.data(), buoy.data(),
                                buoy.size());
            },
            nb::arg("pressure"), nb::arg("height"), nb::arg("buoyancy"),
            R"pbdoc(
Searches the buoyancy array for the LFC and EL combination that results
in the maximum amount of CAPE in the given profile. The buoyancy array 
is typically computed by calling Parcel.lift_parcel. Once the LFC and 
EL are found, the value are set in Parcel.lfc_pres and Parcel.eql_pres.

Parameters:
    pressure: 1D NumPy array of pressure values (Pa)
    height: 1D NumPy array of height values (meters)
    Buoyancy: 1D NumPy array of buoyancy values (m/s^2)
        )pbdoc")
        .def(
            "cape_cinh",
            [](sharp::Parcel& pcl, const_prof_arr_t pres, const_prof_arr_t hght,
               const_prof_arr_t buoy) {
                pcl.cape_cinh(pres.data(), hght.data(), buoy.data(),
                              buoy.size());
                return std::make_tuple(pcl.cape, pcl.cinh);
            },
            nb::arg("pres"), nb::arg("hght"), nb::arg("buoy"),
            R"pbdoc(
Assuming that Parcel.lift_parcel has been called, cape_cinh
will integrate the area between the LFC and EL to compute CAPE,
and integrate the area between the LPL and LCL to compute CINH.

The results are stored in Parcel.cape and Parcel.cinh

Parameters:
    pres: 1D NumPy array of pressure values (Pa)
    hght: 1D NumPy array of height values (Pa)
    buoy: 1D NumPy array of buoyancy values (m/s^2)
        )pbdoc")
        .def_static("surface_parcel", &sharp::Parcel::surface_parcel,
                    nb::arg("pressure"), nb::arg("temperature"),
                    nb::arg("dewpoint"), R"pbdoc(
Given input values of surface pressure, temperature, and dewpoint 
temperature, construct and return a Surface Based Parcel.

Parameters:
    pressure: Surface pressure (Pa)
    temperature: Surface temperature (K)
    dewpoint: Surface dewpoint (K)

Returns:
    Parcel with surface values
                    )pbdoc")
        .def_static(
            "mixed_layer_parcel",
            [](sharp::PressureLayer& mix_layer, const_prof_arr_t pressure,
               const_prof_arr_t potential_temperature,
               const_prof_arr_t mixing_ratio) {
                return sharp::Parcel::mixed_layer_parcel(
                    mix_layer, pressure.data(), nullptr,
                    potential_temperature.data(), mixing_ratio.data(),
                    pressure.size());
            },
            nb::arg("mix_layer"), nb::arg("pressure"),
            nb::arg("potential_temperature"), nb::arg("mixing_ratio"),
            R"pbdoc(
Given input arrays of pressure, potential temperature, and water vapor mixing ratio, as well as a 
defined PressureLayer, compute and return a mixed-layer Parcel.

Parameters:
    mix_layer: PressureLayer over which to compute a mixed-layer parcel 
    pressure: 1D NumPy array of profile pressure values (Pa)
    potential_temperature: 1D NumPy array of profile potential temperature (K)
    mixing_ratio: 1D NumPy array of profile water vapor mixing ratio (unitless)

Returns:
    Parcel with mixed layer values

)pbdoc")
        .def_static(
            "mixed_layer_parcel",
            [](sharp::HeightLayer& mix_layer, const_prof_arr_t pressure,
               const_prof_arr_t height, const_prof_arr_t potential_temperature,
               const_prof_arr_t mixing_ratio) {
                return sharp::Parcel::mixed_layer_parcel(
                    mix_layer, pressure.data(), height.data(),
                    potential_temperature.data(), mixing_ratio.data(),
                    height.size());
            },
            nb::arg("mix_layer"), nb::arg("pressure"), nb::arg("height"),
            nb::arg("potential_temperature"), nb::arg("mixing_ratio"),
            R"pbdoc(
Given input arrays of pressure, potential temperature, and water vapor mixing ratio, as well as a 
defined PressureLayer, compute and return a mixed-layer Parcel.

Parameters:
    mix_layer: HeightLayer over which to compute a mixed-layer parcel 
    pressure: 1D NumPy array of profile pressure values (Pa)
    height: 1D NumPy array of profile height values (meters)
    potential_temperature: 1D NumPy array of profile potential temperature (K)
    mixing_ratio: 1D NumPy array of profile water vapor mixing ratio (unitless)

Returns:
    Parcel with mixed layer values
        )pbdoc")
        .def_static(
            "most_unstable_parcel",
            [](sharp::PressureLayer& layer, sharp::lifter_cm1& lifter,
               const_prof_arr_t pressure, const_prof_arr_t height,
               const_prof_arr_t temperature, const_prof_arr_t virtemp,
               const_prof_arr_t dewpoint) {
                const std::ptrdiff_t NZ = height.size();
                float* pcl_virtemp_arr = new float[NZ];
                float* pcl_buoy_arr = new float[NZ];

                sharp::Parcel max_pcl = sharp::Parcel::most_unstable_parcel(
                    layer, lifter, pressure.data(), height.data(),
                    temperature.data(), virtemp.data(), dewpoint.data(),
                    pcl_virtemp_arr, pcl_buoy_arr, NZ);

                delete[] pcl_virtemp_arr;
                delete[] pcl_buoy_arr;

                return max_pcl;
            },
            nb::arg("layer"), nb::arg("lifter"), nb::arg("pressure"),
            nb::arg("height"), nb::arg("temperature"),
            nb::arg("virtual_temperature"), nb::arg("dewpoint"),
            R"pbdoc(
Given input arrays of pressure, height, temperature, virtual temperature,
and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

Parameters:
    layer: PressureLayer for which to search for the Most Unstable Parcel
    lifter: Parcel lifting routine to use for moist ascent
    pressure: 1D NumPy array of profile pressure values (Pa)
    height: 1D NumPy array of profile height values (meters)
    temperature: 1D NumPy array of profile temperature values (K)
    virtual_temperature: 1D NumPy array of profile virtual temperature values (K)
    dewpoint: 1D NumPy array of profile dewpoint values (K)

Returns:
    Parcel with most-unstable values
        )pbdoc")
        .def_static(
            "most_unstable_parcel",
            [](sharp::HeightLayer& layer, sharp::lifter_cm1& lifter,
               const_prof_arr_t pressure, const_prof_arr_t height,
               const_prof_arr_t temperature, const_prof_arr_t virtemp,
               const_prof_arr_t dewpoint) {
                const std::ptrdiff_t NZ = height.size();

                float* pcl_virtemp_arr = new float[NZ];
                float* pcl_buoy_arr = new float[NZ];

                sharp::Parcel max_pcl = sharp::Parcel::most_unstable_parcel(
                    layer, lifter, pressure.data(), height.data(),
                    temperature.data(), virtemp.data(), dewpoint.data(),
                    pcl_virtemp_arr, pcl_buoy_arr, NZ);

                delete[] pcl_virtemp_arr;
                delete[] pcl_buoy_arr;

                return max_pcl;
            },
            nb::arg("layer"), nb::arg("lifter"), nb::arg("pressure"),
            nb::arg("height"), nb::arg("temperature"),
            nb::arg("virtual_temperature"), nb::arg("dewpoint"),
            R"pbdoc(
Given input arrays of pressure, height, temperature, virtual temperature,
and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

Parameters:
    layer: HeightLayer for which to search for the Most Unstable Parcel
    lifter: Parcel lifting routine to use for moist ascent
    pressure: 1D NumPy array of profile pressure values (Pa)
    height: 1D NumPy array of profile height values (meters)
    temperature: 1D NumPy array of profile temperature values (K)
    virtual_temperature: 1D NumPy array of profile virtual temperature values (K)
    dewpoint: 1D NumPy array of profile dewpoint values (K)

Returns:
    Parcel with most-unstable values
        )pbdoc")
        .def_static(
            "most_unstable_parcel",
            [](sharp::PressureLayer& layer, sharp::lifter_wobus& lifter,
               const_prof_arr_t pressure, const_prof_arr_t height,
               const_prof_arr_t temperature, const_prof_arr_t virtemp,
               const_prof_arr_t dewpoint) {
                const std::ptrdiff_t NZ = height.size();

                float* pcl_virtemp_arr = new float[NZ];
                float* pcl_buoy_arr = new float[NZ];

                sharp::Parcel max_pcl = sharp::Parcel::most_unstable_parcel(
                    layer, lifter, pressure.data(), height.data(),
                    temperature.data(), virtemp.data(), dewpoint.data(),
                    pcl_virtemp_arr, pcl_buoy_arr, NZ);

                delete[] pcl_virtemp_arr;
                delete[] pcl_buoy_arr;

                return max_pcl;
            },
            nb::arg("layer"), nb::arg("lifter"), nb::arg("pressure"),
            nb::arg("height"), nb::arg("temperature"),
            nb::arg("virtual_temperature"), nb::arg("dewpoint"),
            R"pbdoc(
Given input arrays of pressure, height, temperature, virtual temperature,
and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

Parameters:
    layer: PressureLayer for which to search for the Most Unstable Parcel
    lifter: Parcel lifting routine to use for moist ascent
    pressure: 1D NumPy array of profile pressure values (Pa)
    height: 1D NumPy array of profile height values (meters)
    temperature: 1D NumPy array of profile temperature values (K)
    virtual_temperature: 1D NumPy array of profile virtual temperature values (K)
    dewpoint: 1D NumPy array of profile dewpoint values (K)

Returns:
    Parcel with most-unstable values
        )pbdoc")
        .def_static(
            "most_unstable_parcel",
            [](sharp::HeightLayer& layer, sharp::lifter_wobus& lifter,
               const_prof_arr_t pressure, const_prof_arr_t height,
               const_prof_arr_t temperature, const_prof_arr_t virtemp,
               const_prof_arr_t dewpoint) {
                const std::ptrdiff_t NZ = height.size();

                float* pcl_virtemp_arr = new float[NZ];
                float* pcl_buoy_arr = new float[NZ];

                sharp::Parcel max_pcl = sharp::Parcel::most_unstable_parcel(
                    layer, lifter, pressure.data(), height.data(),
                    temperature.data(), virtemp.data(), dewpoint.data(),
                    pcl_virtemp_arr, pcl_buoy_arr, NZ);

                delete[] pcl_virtemp_arr;
                delete[] pcl_buoy_arr;

                return max_pcl;
            },
            nb::arg("layer"), nb::arg("lifter"), nb::arg("pressure"),
            nb::arg("height"), nb::arg("temperature"),
            nb::arg("virtual_temperature"), nb::arg("dewpoint"),
            R"pbdoc(
Given input arrays of pressure, height, temperature, virtual temperature,
and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

Parameters:
    layer: HeightLayer for which to search for the Most Unstable Parcel
    lifter: Parcel lifting routine to use for moist ascent
    pressure: 1D NumPy array of profile pressure values (Pa)
    height: 1D NumPy array of profile height values (meters)
    temperature: 1D NumPy array of profile temperature values (K)
    virtual_temperature: 1D NumPy array of profile virtual temperature values (K)
    dewpoint: 1D NumPy array of profile dewpoint values (K)

Returns:
    Parcel with most-unstable values
        )pbdoc");
}

#endif
