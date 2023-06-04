#include <SHARPlib/CWrap/params/convective_wrap.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/params/convective.h>

#include <stdlib.h>

void sharp_effective_inflow_layer(const float* pressure, const float* height,
                                  const float* temperature,
                                  const float* dewpoint,
                                  const float* virtemp_arr,
                                  const float* buoy_arr, const int N,
                                  sharp_PressureLayer_t* elyr,
                                  float cape_thresh, float cinh_thresh) {
    if (elyr == NULL) return;
    sharp::PressureLayer* l = static_cast<sharp::PressureLayer*>(elyr->obj);
    sharp::PressureLayer out = sharp::effective_inflow_layer(
        pressure, height, temperature, dewpoint, virtemp_arr, buoy_arr, N,
        cape_thresh, cinh_thresh);

    l->bottom = out.bottom;
    l->top = out.top;
}

void sharp_storm_motion_bunkers_np(const float* pressure, const float* height,
                                   const float* u_wind, const float* v_wind,
                                   const int N,
                                   sharp_HeightLayer_t* mn_wind_lyr_agl,
                                   sharp_HeightLayer_t* wind_shr_lyr_agl,
                                   sharp_WindComponents_t* storm_motion,
                                   int leftMover, int pressureWeighted) {
    if ((mn_wind_lyr_agl == NULL) || (wind_shr_lyr_agl == NULL) ||
        (storm_motion == NULL))
        return;
    sharp::HeightLayer* mnl =
        static_cast<sharp::HeightLayer*>(mn_wind_lyr_agl->obj);
    sharp::HeightLayer* shl =
        static_cast<sharp::HeightLayer*>(wind_shr_lyr_agl->obj);
    sharp::WindComponents* stm =
        static_cast<sharp::WindComponents*>(storm_motion->obj);

    sharp::WindComponents out =
        sharp::storm_motion_bunkers(pressure, height, u_wind, v_wind, N, *mnl,
                                    *shl, leftMover, pressureWeighted);
    stm->u = out.u;
    stm->v = out.v;
}

void sharp_storm_motion_bunkers(
    const float* pressure, const float* height, const float* u_wind,
    const float* v_wind, const int N, sharp_PressureLayer_t* eff_infl_lyr,
    sharp_Parcel_t* pcl, sharp_WindComponents_t* storm_motion, int leftMover) {

    if ((storm_motion == NULL) || (pcl == NULL) || (eff_infl_lyr == NULL)) return;

    sharp::PressureLayer* eil =
        static_cast<sharp::PressureLayer*>(eff_infl_lyr->obj);
	sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    sharp::WindComponents* stm =
        static_cast<sharp::WindComponents*>(storm_motion->obj);

    sharp::WindComponents out = sharp::storm_motion_bunkers(
        pressure, height, u_wind, v_wind, N, eil, pcl, leftMover);
    stm->u = out.u;
    stm->v = out.v;
}

float sharp_entrainment_cape(const float* pressure, const float* height,
                             const float* temperature, const float* mse_arr,
                             const float* u_wind, const float* v_wind,
                             const int N, sharp_Parcel_t* pcl) {
    if ((pcl == NULL)) return sharp::MISSING;
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    return sharp::entrainment_cape(pressure, height, temperature, mse_arr,
                                   u_wind, v_wind, N, pc);
}

float sharp_energy_helicity_index(float cape, float helicity) {
    return sharp::energy_helicity_index(cape, helicity);
}

float sharp_supercell_composite_parameter(float mu_cape, float eff_srh,
                                          float eff_shear) {
    return sharp::supercell_composite_parameter(mu_cape, eff_srh, eff_shear);
}

float sharp_significant_tornado_parameter(sharp_Parcel_t* pcl,
                                          float lcl_hght_agl, float helicity,
                                          float bulk_shear) {
    if ((pcl == NULL)) return sharp::MISSING;
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    return sharp::significant_tornado_parameter(*pc, lcl_hght_agl, helicity,
                                                bulk_shear);
}
