#ifndef SHARPLIB_PARCEL_BINDINGS
#define SHARPLIB_PARCEL_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

// clang-format on
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
}

#endif
