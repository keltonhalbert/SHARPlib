#ifndef SHARPLIB_CONSTANTS_BINDINGS
#define SHARPLIB_CONSTANTS_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>

// clang-format on
#include <SHARPlib/constants.h>

namespace nb = nanobind;

void make_constants_bindings(nb::module_ m) {
    nb::module_ m_const =
        m.def_submodule("constants",
                        "Sounding and Hodograph Analysis and Research Program "
                        "Library (SHARPlib) :: Constants");

    m_const.attr("RGAS") = sharp::RGAS;
    m_const.attr("MDRY") = sharp::MDRY;
    m_const.attr("MWATER") = sharp::MWATER;
    m_const.attr("DRYAIR_SPEC_HEAT_RATIO") = sharp::DRYAIR_SPEC_HEAT_RATIO;
    m_const.attr("ZEROCNK") = sharp::ZEROCNK;
    m_const.attr("HPA_TO_PA") = sharp::HPA_TO_PA;
    m_const.attr("PA_TO_HPA") = sharp::PA_TO_HPA;
    m_const.attr("THETA_REF_PRESSURE") = sharp::THETA_REF_PRESSURE;
    m_const.attr("MISSING") = sharp::MISSING;
    m_const.attr("GRAVITY") = sharp::GRAVITY;
    m_const.attr("TOL") = sharp::TOL;
    m_const.attr("PI") = sharp::PI;
    m_const.attr("RDGAS") = sharp::RDGAS;
    m_const.attr("RVGAS") = sharp::RVGAS;
    m_const.attr("EPSILON") = sharp::EPSILON;
    m_const.attr("CP_DRYAIR") = sharp::CP_DRYAIR;
    m_const.attr("CP_VAPOR") = sharp::CP_VAPOR;
    m_const.attr("CP_LIQUID") = sharp::CP_LIQUID;
    m_const.attr("CP_ICE") = sharp::CP_ICE;
    m_const.attr("EXP_LV") = sharp::EXP_LV;
    m_const.attr("EXP_LS") = sharp::EXP_LS;
    m_const.attr("LV1") = sharp::LV1;
    m_const.attr("LV2") = sharp::LV2;
    m_const.attr("LS1") = sharp::LS1;
    m_const.attr("LS2") = sharp::LS2;
    m_const.attr("ROCP") = sharp::ROCP;
    m_const.attr("GAMMA_D") = sharp::GAMMA_D;
    m_const.attr("PRANDTL") = sharp::PRANDTL;
    m_const.attr("VKSQ") = sharp::VKSQ;
    m_const.attr("VAPPRES_REF") = sharp::VAPPRES_REF;
}

#endif
