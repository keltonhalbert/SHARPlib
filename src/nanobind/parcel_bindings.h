#ifndef SHARPLIB_PARCEL_BINDINGS
#define SHARPLIB_PARCEL_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/shared_ptr.h>

// clang-format on
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <fmt/core.h>

#include <memory>

#include "binding_utils.h"
#include "sharplib_types.h"

namespace nb = nanobind;

template <typename Lft>
void bind_lift_parcel(nb::class_<sharp::Parcel>& cls, const char* doc_template,
                      const char* lifter_name) {
    std::string doc = fmt::format(doc_template, lifter_name);
    cls.def(
        "lift_parcel",
        [](sharp::Parcel& pcl, Lft& lifter, const_prof_arr_t pressure) {
            return make_output_array(pressure.size(), [&](float* out) {
                pcl.lift_parcel(lifter, pressure.data(), out, pressure.size());
            });
        },
        nb::arg("lifter"), nb::arg("pressure"), doc.c_str());
}

template <typename Layer, typename Lft>
void bind_mu_parcel(nb::class_<sharp::Parcel>& cls, const char* doc_template,
                    const char* layer_name, const char* lifter_name) {
    std::string doc = fmt::format(doc_template, layer_name, lifter_name);
    cls.def_static(
        "most_unstable_parcel",
        [](Layer& layer, Lft& lifter, const_prof_arr_t pressure,
           const_prof_arr_t height, const_prof_arr_t temperature,
           const_prof_arr_t virtemp, const_prof_arr_t dewpoint) {
            check_equal_sizes(pressure, temperature, dewpoint, height, virtemp);

            const std::size_t NZ = height.size();
            auto pcl_virtemp_arr = std::make_unique<float[]>(NZ);
            auto pcl_buoy_arr = std::make_unique<float[]>(NZ);

            sharp::Parcel max_pcl = sharp::Parcel::most_unstable_parcel(
                layer, lifter, pressure.data(), height.data(),
                temperature.data(), virtemp.data(), dewpoint.data(),
                pcl_virtemp_arr.get(), pcl_buoy_arr.get(), NZ);

            return max_pcl;
        },
        nb::arg("layer"), nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("height"), nb::arg("temperature"),
        nb::arg("virtual_temperature"), nb::arg("dewpoint"), doc.c_str());
}

template <typename Layer>
void bind_ml_parcel(nb::class_<sharp::Parcel>& cls, const char* doc) {
    if constexpr (Layer::coord == sharp::LayerCoordinate::height) {
        cls.def_static(
            "mixed_layer_parcel",
            [](Layer& mix_layer, const_prof_arr_t pressure,
               const_prof_arr_t height, const_prof_arr_t potential_temperature,
               const_prof_arr_t mixing_ratio) {
                check_equal_sizes(pressure, potential_temperature, mixing_ratio,
                                  height);
                return sharp::Parcel::mixed_layer_parcel(
                    mix_layer, pressure.data(), height.data(),
                    potential_temperature.data(), mixing_ratio.data(),
                    pressure.size());
            },
            nb::arg("mix_layer"), nb::arg("pressure"), nb::arg("height"),
            nb::arg("potential_temperature"), nb::arg("mixing_ratio"), doc);
    } else {
        cls.def_static(
            "mixed_layer_parcel",
            [](Layer& mix_layer, const_prof_arr_t pressure,
               const_prof_arr_t potential_temperature,
               const_prof_arr_t mixing_ratio) {
                check_equal_sizes(pressure, potential_temperature,
                                  mixing_ratio);
                return sharp::Parcel::mixed_layer_parcel(
                    mix_layer, pressure.data(), nullptr,
                    potential_temperature.data(), mixing_ratio.data(),
                    pressure.size());
            },
            nb::arg("mix_layer"), nb::arg("pressure"),
            nb::arg("potential_temperature"), nb::arg("mixing_ratio"), doc);
    }
}

template <typename Lft>
void bind_lut_data_class(nb::class_<sharp::lut_data<Lft>>& cls) {
    cls.def(nb::init<Lft, float, float, float, float, std::size_t,
                     std::size_t>(),
            nb::arg("lifter"), nb::arg("pmin") = 5000.0f,
            nb::arg("pmax") = 110000.0f, nb::arg("thte_min") = 210.0f,
            nb::arg("thte_max") = 430.0f, nb::arg("n_logp") = 201,
            nb::arg("n_thetae") = 221)
        .def_ro("pres_min", &sharp::lut_data<Lft>::pres_min)
        .def_ro("pres_max", &sharp::lut_data<Lft>::pres_max)
        .def_ro("thte_min", &sharp::lut_data<Lft>::thetae_min)
        .def_ro("thte_max", &sharp::lut_data<Lft>::thetae_max)
        .def_ro("num_logp", &sharp::lut_data<Lft>::num_logp)
        .def_ro("num_thetae", &sharp::lut_data<Lft>::num_thetae);
}

template <typename Lft>
void bind_lut_data(nb::module_& mod, const char* doc_template,
                   const char* lifter_name, const char* return_name) {
    std::string doc = fmt::format(doc_template, lifter_name, return_name);
    mod.def(
        "lut_data",
        [](Lft& lifter, float pmin, float pmax, float thte_min, float thte_max,
           std::size_t n_logp, std::size_t n_thetae) {
            return sharp::lut_data(lifter, pmin, pmax, thte_min, thte_max,
                                   n_logp, n_thetae);
        },
        nb::arg("lifter"), nb::arg("pmin") = 5000.0f,
        nb::arg("pmax") = 110000.0f, nb::arg("thte_min") = 210.0f,
        nb::arg("thte_max") = 430.0f, nb::arg("n_logp") = 201,
        nb::arg("n_thetae") = 221, doc.c_str());
}

template <typename Lft>
void bind_lifter_lut_class(nb::class_<sharp::lifter_lut<Lft>>& cls) {
    cls.def(nb::init<std::shared_ptr<const sharp::lut_data<Lft>>>(),
            nb::arg("data"))
        .def_ro_static("lift_from_lcl", &sharp::lifter_lut<Lft>::lift_from_lcl,
                       R"pbdoc(
A static flag that helps the parcel lifting functions know where to lift from.
)pbdoc")
        .def("setup", &sharp::lifter_lut<Lft>::setup, nb::arg("lcl_pres"),
             nb::arg("lcl_tmpk"), R"pbdoc(
Performs a setup step based on the LCL attributes.

Computes the fractional index needed to select the correct pseudoadiabat
for lookup. If the LCL is outside the table bounds, it falls back to the
direct solver.

Parameters
----------
lcl_pres : float
    The LCL pressure (Pa)
lcl_tmpk : float
    The LCL temperature (K)

Returns
-------
None
)pbdoc")
        .def("__call__", &sharp::lifter_lut<Lft>::operator(), nb::is_operator(),
             nb::arg("pres"), nb::arg("tmpk"), nb::arg("new_pres"), R"pbdoc(
Performs LUT interpolation to lift a parcel moist adiabatically.

Parameters
----------
pres : float
    Parcel pressure (Pa)
tmpk : float
    Parcel temperature (K)
new_pres : float
    Final level of parcel after lift (Pa)

Returns
-------
float
    The temperature of the lifted parcel (K)
)pbdoc")
        .def("parcel_virtual_temperature",
             &sharp::lifter_lut<Lft>::parcel_virtual_temperature,
             nb::arg("pres"), nb::arg("tmpk"), R"pbdoc(
Computes the virtual temperature of the parcel (after saturation).

Parameters
----------
pres : float
    Parcel pressure (Pa)
tmpk : float
    Parcel temperature (K)

Returns
-------
float
    The virtual temperature of the parcel (K)
)pbdoc");
}

template <typename Lft>
void bind_lifter_lut(nb::module_& mod, const char* doc_template,
                     const char* lut_name, const char* return_name) {
    std::string doc = fmt::format(doc_template, lut_name, return_name);
    mod.def(
        "lifter_lut",
        [](std::shared_ptr<const sharp::lut_data<Lft>>& lut) {
            return sharp::lifter_lut<Lft>(std::move(lut));
        },
        nb::arg("lut"), doc.c_str());
}

template <typename Lft>
void bind_lower_parcel(nb::class_<sharp::DowndraftParcel>& cls,
                       const char* doc_template, const char* lifter_name) {
    std::string doc = fmt::format(doc_template, lifter_name);
    cls.def(
        "lower_parcel",
        [](sharp::DowndraftParcel& pcl, Lft& lifter,
           const_prof_arr_t pressure) {
            return make_output_array(pressure.size(), [&](float* out) {
                pcl.lower_parcel(lifter, pressure.data(), out, pressure.size());
            });
        },
        nb::arg("lifter"), nb::arg("pressure"), doc.c_str());
}

inline void make_parcel_bindings(nb::module_ m) {
    nb::module_ m_parcel = m.def_submodule(
        "parcel",
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Parcel Lifting Routines");

    const char* lift_parcel_template_doc =
        R"pbdoc(
Lifts a Parcel dry adiabatically from its LPL to its LCL dry
adiabatically, and then moist adiabatically from the LCL to 
the top of the profile. The moist adiabat used is determined
by the type of lifting functor passed to the function (i.e.
lifter_wobus or lifter_cm1).

Parameters
----------
lifter : {}
    An instantiated lifter functor
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of Pressure levels for lifting (Pa)

Returns
-------
numpy.ndarray[dtype=float32]
    A 1D NumPy array of parcel virtual temperature values (K)
        )pbdoc";

    const char* mu_pcl_template_doc =
        R"pbdoc(
Given input arrays of pressure, height, temperature, virtual temperature,
and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

Parameters
----------
layer : {}
    Layer for which to search for the Most Unstable Parcel
lifter : {}
    Parcel lifting routine to use for moist ascent
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile pressure values (Pa)
height : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile height values (meters)
temperature : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile temperature values (K)
virtual_temperature : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile virtual temperature values (K)
dewpoint : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile dewpoint values (K)

Returns
-------
nwsspc.sharp.calc.parcel.Parcel
    Parcel with most-unstable values
        )pbdoc";

    const char* ml_pcl_pres_doc =
        R"pbdoc(
Given input arrays of pressure, potential temperature, 
and water vapor mixing ratio, as well as a defined PressureLayer, 
compute and return a mixed-layer Parcel.

Parameters
----------
mix_layer : nwsspc.sharp.calc.layer.PressureLayer 
    PressureLayer over which to compute a mixed-layer parcel 
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile pressure values (Pa)
potential_temperature : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile potential temperature (K)
mixing_ratio : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile water vapor mixing ratio (unitless)

Returns
-------
nwsspc.sharp.calc.parcel.Parcel
    Parcel with mixed layer values

)pbdoc";

    const char* ml_pcl_hght_doc =
        R"pbdoc(
Given input arrays of pressure, potential temperature, and water 
vapor mixing ratio, as well as a defined PressureLayer, compute 
and return a mixed-layer Parcel.

Parameters
----------
mix_layer : nwsspc.sharp.calc.layer.HeightLayer 
    HeightLayer over which to compute a mixed-layer parcel 
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile pressure values (Pa)
height : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile height values (meters)
potential_temperature : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile potential temperature (K)
mixing_ratio : numpy.ndarray[dtype=float32] 
    1D NumPy array of profile water vapor mixing ratio (unitless)

Returns
-------
nwsspc.sharp.calc.parcel.Parcel
    Parcel with mixed layer values
        )pbdoc";

    const char* lut_data_template_doc =
        R"pbdoc(
Constructs the lookup table (LUT) data for a parcel lifter.

Parameters
----------
lifter : {}
pmin : float 
    The minimum pressure of the lookup table (Pa)
pmax : float 
    The maximum pressure of the lookup table (Pa)
thte_min : float 
    The minimum thetae of the lookup table (K)
thte_max : float 
    The maximum thetae of the lookup table (K)
n_logp : uint
    The number of logp levels for the lookup table
n_thetae : uint
    The number of thetae levels for the lookup table

Returns
-------
{}
   A lookup table for the given lifter
)pbdoc";

    const char* lifter_lut_class_template_doc =
        R"pbdoc(
A parcel lifter functor that uses a pseudoadiabatic lookup table (LUT)
for fast moist adiabatic ascent calculations.

Instead of directly solving the moist ascent ODEs, this lifter uses
bilinear interpolation of a precomputed lookup table to determine the
parcel temperature. If the parcel's LCL falls outside the table bounds,
it falls back to the direct {} solver automatically.

LUT based parcel ascent only works for pseudoadiabats. Constructing the
LUT with a reversible adiabat type will result in an error being thrown.

Parameters
----------
data : {}
    A shared lookup table constructed with a lifter_{} instance.
)pbdoc";

    const char* lifter_lut_template_doc =
        R"pbdoc(
Constructs the parcel lifter from a LUT.

Parameters
----------
lut : {}

Returns
-------
{}
)pbdoc";

    const char* lower_parcel_template_doc =
        R"pbdoc(
Lowers a saturated nwsspc.sharp.calc.parcel.DowndraftParcel moist 
adiabatically from its LPL to the surface. The moist adiabat used 
is determined by the type of lifting functor passed to the function 
(i.e. lifter_wobus or lifter_cm1).

Unlike nwsspc.sharp.calc.parcel.Parcel.lift_parcel, the virtual 
temperature correction is not used for downdraft parcels.

Parameters
----------
lifter : {}
    An instantiated lifter functor
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of Pressure levels for lifting (Pa)

Returns
-------
numpy.ndarray[dtype=float32]
    A 1D NumPy array of parcel temperature values (K)
        )pbdoc";

    nb::class_<sharp::lifter_wobus>(m_parcel, "lifter_wobus",
                                    R"pbdoc(
A functor that calls the Wobus Wetlift for computation of moist adiabats.

This is used to wrap the Wobus Wetlift function for parcel lifting 
routines. Functors -- classes with their operator() overloaded -- 
are used so that functions can be passed to templates in a way that 
the compiler can still optimize, rather than using function pointers
or lambdas. 

Specifically, this functor is designed to be passed as a template 
argument to nwsspc.sharp.calc.parcel.Parcel.lift_parcel, so that 
the method of computing moist adiabats can be changed without 
changing the overall parcel lifting code. The reason this is awesome 
is that the compiler can still optimize and inline this code, while 
the user can configure the parcel lifting algorithm to their 
specifications. 

)pbdoc")
        .def(nb::init<>())
        .def_ro_static("lift_from_lcl", &sharp::lifter_wobus::lift_from_lcl,
                       R"pbdoc(
A static flag that helps the parcel lifting functions know where to lift from.
)pbdoc")
        .def_ro_static("ma_type", &sharp::lifter_wobus::ma_type,
                       R"pbdoc(
A static flag for the type of solver this is 
)pbdoc")
        .def_rw("converge", &sharp::lifter_wobus::converge, R"pbdoc(
The iterative convergence criteria (K)
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

Parameters
----------
pres : float 
    Parcel pressure (Pa)
tmpk : float 
    Parcel temperature (K)
new_pres : float 
    Final level of parcel after lift (Pa)

Returns
-------
float
    The temperature of the lifted parcel (K)
             )pbdoc")
        .def("parcel_virtual_temperature",
             &sharp::lifter_wobus::parcel_virtual_temperature, nb::arg("pres"),
             nb::arg("tmpk"), R"pbdoc(
Computes the virtual temperature of the parcel (after saturation).

Parameters
----------
pres : float 
    Parcel pressure (Pa)
tmpk : float 
    Parcel temperature (K)

Returns
-------
float
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
The type of moist adiabat to use, defined by nwsspc.sharp.calc.thermo.adiabat
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

Parameters
----------
lcl_pres : float 
    The LCL pressure (Pa)
lcl_tmpk : float 
    The LCL temperature (K)

Returns
-------
None
             )pbdoc")
        .def("__call__", &sharp::lifter_cm1::operator(), nb::is_operator(),
             nb::arg("pres"), nb::arg("tmpk"), nb::arg("new_pres"),
             R"pbdoc(
Lifts a parcel moist adiabatically/pseudoadiabatically using
nwsspc.sharp.calc.thermo.moist_adiabat_cm1.

Parameters
----------
pres : float 
    Parcel pressure (Pa)
tmpk : float 
    Parcel temperature (K)
new_pres : float 
    Final level of parcel after lift (Pa)

Returns
-------
float
    The temperature of the lifted parcel (K)
             )pbdoc")

        .def("parcel_virtual_temperature",
             &sharp::lifter_cm1::parcel_virtual_temperature, nb::arg("pres"),
             nb::arg("tmpk"), R"pbdoc(
Computes the virtual temperature of the parcel (after saturation).

Parameters
----------
pres : float 
    Parcel pressure (Pa)
tmpk : float 
    Parcel temperature (K)

Returns
-------
float
    The virtual temperature of the parcel (K)
             )pbdoc");

    for_each_base_lifter([&](auto tag, const char* suffix,
                             const char* display_name, const char* lifter_fqn,
                             const char* lut_data_fqn,
                             const char* lifter_lut_fqn) {
        using Lft = typename decltype(tag)::type;

        std::string lut_data_cls_name = fmt::format("lut_data_{}", suffix);
        auto lut_data_cls = nb::class_<sharp::lut_data<Lft>>(
            m_parcel, lut_data_cls_name.c_str());
        bind_lut_data_class(lut_data_cls);

        bind_lut_data<Lft>(m_parcel, lut_data_template_doc, lifter_fqn,
                           lut_data_fqn);

        std::string lifter_lut_cls_name = fmt::format("lifter_lut_{}", suffix);
        std::string lifter_lut_cls_doc = fmt::format(
            lifter_lut_class_template_doc, display_name, lut_data_fqn, suffix);
        auto lifter_lut_cls = nb::class_<sharp::lifter_lut<Lft>>(
            m_parcel, lifter_lut_cls_name.c_str(), lifter_lut_cls_doc.c_str());
        bind_lifter_lut_class(lifter_lut_cls);

        bind_lifter_lut<Lft>(m_parcel, lifter_lut_template_doc, lut_data_fqn,
                             lifter_lut_fqn);
    });

    nb::enum_<sharp::LPL>(m_parcel, "LPL")
        .value("SFC", sharp::LPL::SFC, "A Surface Based Parcel")
        .value("FCST", sharp::LPL::FCST, "A Forecast Surface Parcel")
        .value("MU", sharp::LPL::MU, "A Most Unstable Parcel")
        .value("ML", sharp::LPL::ML, "A Mixed-Layer Parcel")
        .value("USR", sharp::LPL::USR, "A User Defined Parcel");

    auto parcel_class =
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
            .def_rw("mpl_pressure", &sharp::Parcel::mpl_pressure, R"pbdoc(
Pressure at the Maximum Parcel Level
                )pbdoc")
            .def_rw("cape", &sharp::Parcel::cape, R"pbdoc(
Parcel Convective Available Potential Energy (J/kg) between the LFC and EL
                )pbdoc")
            .def_rw("cinh", &sharp::Parcel::cinh, R"pbdoc(
Parcel Convective Inhibition (J/kg) between the LFC and EL
                )pbdoc")
            .def(nb::init<>())
            .def(nb::init<const float, const float, const float,
                          const sharp::LPL>(),
                 nb::arg("pressure"), nb::arg("temperature"),
                 nb::arg("dewpoint"), nb::arg("lpl"), R"pbdoc(
Constructor for a Parcel

Parameters
----------
pressure : float 
    Parcel initial pressure (Pa)
temperature : float 
    Parcel initial temperature (K)
dewpoint : float 
    Parcel initial dewpoint (K)
lpl : nwsspc.sharp.calc.parcel.LPL 
    Parcel Lifted Parcel Level (LPL) definition
        )pbdoc")
            .def(
                "find_lfc_el",
                [](sharp::Parcel& pcl, const_prof_arr_t pres,
                   const_prof_arr_t hght, const_prof_arr_t buoy) {
                    check_equal_sizes(pres, hght, buoy);
                    pcl.find_lfc_el(pres.data(), hght.data(), buoy.data(),
                                    buoy.size());
                    return std::make_tuple(pcl.lfc_pressure, pcl.eql_pressure);
                },
                nb::arg("pressure"), nb::arg("height"), nb::arg("buoyancy"),
                R"pbdoc(
Searches the buoyancy array for the LFC and EL combination that results
in the maximum amount of CAPE in the given profile. The buoyancy array 
is typically computed by calling nwsspc.sharp.calc.parcel.Parcel.lift_parcel. 
Once the LFC and EL are found, the value are set in 
nwsspc.sharp.calc.parcel.Parcel.lfc_pres and 
nwsspc.sharp.calc.parcel.Parcel.eql_pres via the provided parcel.

The value of eql_pres is MISSING if there is no qualifying level 
found within the data bounds (e.g. incomplete data, or EL above 
the available data). Any calls to nwsspc.sharp.calc.parcel.Parcel.cape_cinh 
will still compute CAPE without the presence of an EL, using the best-available
data. 

Parameters
----------
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of pressure values (Pa)
height : numpy.ndarray[dtype=float32]
    1D NumPy array of height values (meters)
buoyancy : numpy.ndarray[dtype=float32] 
    1D NumPy array of buoyancy values (m/s^2)

Returns
-------
tuple[float, float]
    (LFC_PRES, EL_PRES)
        )pbdoc")
            .def(
                "maximum_parcel_level",
                [](sharp::Parcel& pcl, const_prof_arr_t pres,
                   const_prof_arr_t hght, const_prof_arr_t buoy) {
                    check_equal_sizes(pres, hght, buoy);
                    return pcl.maximum_parcel_level(pres.data(), hght.data(),
                                                    buoy.data(), pres.size());
                },
                nb::arg("pressure"), nb::arg("height"), nb::arg("buoyancy"),
                R"pbdoc(
Find the pressure of the Maximum Parcel Level (MPL).

The Maximum Parcel Level (MPL) is the level a parcel would reach 
if it expended all of its integrated positive buoyancy past the 
Equilibrium Level. It is found by integrating negatively buoyant 
area above the Equilibrium Level until the integrated negative 
buoyancy is equal in magnitude to the Convective Available 
Potential Energy between the Level of Free Convection and the 
Equilibrium Level. 

For valid calculations, nwsspc.sharp.calc.parcel.Parcel.cape_cinh 
must be called first, or nwsspc.sharp.calc.parcel.Parcel.cape and 
nwsspc.sharp.calc.parcel.Parcel.eql_pressure must be set. 

A value of nwsspc.sharp.calc.constants.MISSING is returned if:
  * CAPE is 0 
  * nwsspc.sharp.calc.parcel.Parcel.eql_pressure is MISSING
  * No valid MPL candidate is found within the profile
    In this scenario, it likely exceeds the top of the available data

In addition to being returned, the result is stored inside of 
nwsspc.sharp.calc.parcel.Parcel.mpl_pressure.

Parameters
----------
pres : numpy.ndarray[dtype=float32] 
    1D NumPy array of pressure values (Pa)
hght : numpy.ndarray[dtype=float32] 
    1D NumPy array of height values (meters)
buoyancy : numpy.ndarray[dtype=float32] 
    1D NumPy array of buoyancy values (m/s^2)

Returns 
-------
float 
    The pressure of the Maximum Parcel Level (Pa)
        )pbdoc")
            .def(
                "cape_cinh",
                [](sharp::Parcel& pcl, const_prof_arr_t pres,
                   const_prof_arr_t hght, const_prof_arr_t buoy) {
                    check_equal_sizes(pres, hght, buoy);
                    pcl.cape_cinh(pres.data(), hght.data(), buoy.data(),
                                  buoy.size());
                    return std::make_tuple(pcl.cape, pcl.cinh);
                },
                nb::arg("pres"), nb::arg("hght"), nb::arg("buoy"),
                R"pbdoc(
Assuming that nwsspc.sharp.calc.parcel.Parcel.lift_parcel has 
been called, cape_cinh will integrate the area between the LFC 
and EL to compute CAPE, and integrate the area between the LPL 
and LCL to compute CINH.

If eql_pressure is MISSING, but lfc_pressure is defined, the 
routine will compute CAPE with the available data despite the 
lack of a defined equilibrium level. This is useful for 
incomplete profile data, or pressure-level data where the EL 
is above the top pressure value.

The results are stored in nwsspc.sharp.calc.parcel.Parcel.cape 
and nwsspc.sharp.calc.parcel.Parcel.cinh via the provided parcel.

Parameters
----------
pres : numpy.ndarray[dtype=float32] 
    1D NumPy array of pressure values (Pa)
hght : numpy.ndarray[dtype=float32] 
    1D NumPy array of height values (meters)
buoyancy : numpy.ndarray[dtype=float32] 
    1D NumPy array of buoyancy values (m/s^2)

Returns
-------
tuple[float, float]
    (CAPE, CINH)
        )pbdoc")
            .def(
                "lifted_index",
                [](sharp::Parcel& pcl, const float plev, const_prof_arr_t pres,
                   const_prof_arr_t vtmpk, const_prof_arr_t pcl_vtmpk) {
                    check_equal_sizes(pres, vtmpk, pcl_vtmpk);
                    return pcl.lifted_index(plev, pres.data(), vtmpk.data(),
                                            pcl_vtmpk.data(), pres.size());
                },
                nb::arg("pres_lev"), nb::arg("pressure"),
                nb::arg("virtual_temperature"),
                nb::arg("parcel_virtual_temperature"),
                R"pbdoc(
Using the parcel and environment virtual temperature, compute the lifted index 
at the given level. The pressure level typically used is 500 hPa. The lifted index 
is the difference between the environment and parcel virtual temperatures at the 
requested level. 

Parameters
----------
pres_lev : float 
    The pressure level (Pa) to compute the lifted index at
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure values (Pa)
virtual_temperature : numpy.ndarray[dtype=float32]
    1D NumPy array of environment virtual temperature values (K)
parcel_virtual_temperature : numpy.ndarray[dtype=float32]
    1D NumPy array of parcel virtual temperature values (K)

Returns 
-------
float 
    The lifted index
        )pbdoc")
            .def_static("surface_parcel", &sharp::Parcel::surface_parcel,
                        nb::arg("pressure"), nb::arg("temperature"),
                        nb::arg("dewpoint"), R"pbdoc(
Given input values of surface pressure, temperature, and dewpoint 
temperature, construct and return a Surface Based Parcel.

Parameters
----------
pressure : float 
    Surface pressure (Pa)
temperature : float 
    Surface temperature (K)
dewpoint : float 
    Surface dewpoint (K)

Returns
-------
nwsspc.sharp.calc.parcel.Parcel
    Parcel with surface values
                    )pbdoc");

    for_each_lifter([&](auto tag, const char* lifter_name) {
        using Lft = typename decltype(tag)::type;
        bind_lift_parcel<Lft>(parcel_class, lift_parcel_template_doc,
                              lifter_name);
    });

    for_each_lifter([&](auto tag, const char* lifter_name) {
        using Lft = typename decltype(tag)::type;
        bind_mu_parcel<sharp::PressureLayer, Lft>(
            parcel_class, mu_pcl_template_doc,
            "nwsspc.sharp.calc.layer.PressureLayer", lifter_name);
        bind_mu_parcel<sharp::HeightLayer, Lft>(
            parcel_class, mu_pcl_template_doc,
            "nwsspc.sharp.calc.layer.HeightLayer", lifter_name);
    });

    bind_ml_parcel<sharp::PressureLayer>(parcel_class, ml_pcl_pres_doc);
    bind_ml_parcel<sharp::HeightLayer>(parcel_class, ml_pcl_hght_doc);

    auto dpcl_class =
        nb::class_<sharp::DowndraftParcel>(m_parcel, "DowndraftParcel", R"pbdoc(
Contains information about a DowndraftParcel's starting level and
thermodynamic attributes, as well as derived computations,
methods for constructing a parcel, and parcel descent routines. 
                              )pbdoc")
            .def_rw("pres", &sharp::DowndraftParcel::pres, R"pbdoc(
DowndraftParcel starting pressure (Pa)
                )pbdoc")
            .def_rw("tmpk", &sharp::DowndraftParcel::tmpk, R"pbdoc(
DowndraftParcel starting temperature (K)
                )pbdoc")
            .def_rw("dwpk", &sharp::DowndraftParcel::dwpk, R"pbdoc(
DowndraftParcel starting dewpoint (K)
                )pbdoc")
            .def_rw("cape", &sharp::DowndraftParcel::cape, R"pbdoc(
DowndraftParcel Convective Available Potential Energy (J/kg)
                )pbdoc")
            .def_rw("cinh", &sharp::DowndraftParcel::cinh, R"pbdoc(
DowndraftParcel Convective Inhibition (J/kg)
                )pbdoc")
            .def(nb::init<>())
            .def(nb::init<const float, const float, const float>(),
                 nb::arg("pressure"), nb::arg("temperature"),
                 nb::arg("dewpoint"),
                 R"pbdoc(
Constructor for a DowndraftParcel

Parameters
----------
pressure : float 
    DowndraftParcel initial pressure (Pa)
temperature : float 
    DowndraftParcel initial temperature (K)
dewpoint : float 
    DowndraftParcel initial dewpoint (K)
        )pbdoc")
            .def(
                "cape_cinh",
                [](sharp::DowndraftParcel& pcl, const_prof_arr_t pres,
                   const_prof_arr_t hght, const_prof_arr_t buoy) {
                    check_equal_sizes(pres, hght, buoy);
                    pcl.cape_cinh(pres.data(), hght.data(), buoy.data(),
                                  buoy.size());
                    return std::make_tuple(pcl.cape, pcl.cinh);
                },
                nb::arg("pres"), nb::arg("hght"), nb::arg("buoy"),
                R"pbdoc(
Assuming that nwsspc.sharp.calc.parcel.DowndraftParcel.lower_parcel 
has been called, cape_cinh will integrate the area between the LPL
and the surface to compute downdraft CAPE and downdraft CINH.

The results are stored in nwsspc.sharp.calc.parcel.DowndraftParcel.cape 
and nwsspc.sharp.calc.parcel.DowndraftParcel.cinh via the provided parcel.

Parameters
----------
pres : numpy.ndarray[dtype=float32] 
    1D NumPy array of pressure values (Pa)
hght : numpy.ndarray[dtype=float32] 
    1D NumPy array of height values (meters)
buoyancy : numpy.ndarray[dtype=float32] 
    1D NumPy array of buoyancy values (m/s^2)

Returns
-------
tuple[float, float]
    (DCAPE, DCINH)
        )pbdoc")
            .def_static(
                "min_thetae",
                [](sharp::PressureLayer& search_layer,
                   const_prof_arr_t pressure, const_prof_arr_t temperature,
                   const_prof_arr_t dewpoint, const_prof_arr_t thetae,
                   const float mean_depth) {
                    check_equal_sizes(pressure, temperature, dewpoint, thetae);
                    return sharp::DowndraftParcel::min_thetae(
                        search_layer, pressure.data(), temperature.data(),
                        dewpoint.data(), thetae.data(), pressure.size(),
                        mean_depth);
                },
                nb::arg("search_layer"), nb::arg("pressure"),
                nb::arg("temperature"), nb::arg("dewpoint"), nb::arg("thetae"),
                nb::arg("mean_depth") = 10000.0f,
                R"pbdoc(
Define a downdraft parcel. 

Defines a downdraft parcel within a given search layer. 
The downdraft parcel is defined as the minimum layer-mean 
equivalent potential temperature (Theta-E) within the 
search layer. Typical values are to search within the lowest
400 hPa of the profile, and a mean depth of 100 hPa. 

Parameters 
----------
search_layer : nwsspc.sharp.calc.layer.PressureLayer 
    The layer over which to search for the downdraft parcel 
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure (Pa)
temperature : numpy.ndarray[dtype=float32]
    1D NumPy array of temperature (K)
dewpoint : numpy.ndarray[dtype=float32]
    1D NumPy array of dewpoint (K)
thetae : numpy.ndarray[dtype=float32]
    1D NumPy array of thetae (K)
mean_depth : float
    The layer depth for calculating mean thetae.

Returns 
-------
nwsspc.sharp.calc.parcel.DowndraftParcel 
    Downdraft Parcel
    )pbdoc");

    for_each_lifter([&](auto tag, const char* lifter_name) {
        using Lft = typename decltype(tag)::type;
        bind_lower_parcel<Lft>(dpcl_class, lower_parcel_template_doc,
                               lifter_name);
    });
}

#endif
