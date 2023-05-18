#include <SHARPlib/CWrap/params_wrap.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/params.h>

void sharp_effective_inflow_layer(sharp_Profile_t* prof,
                                  sharp_PressureLayer_t* elyr,
                                  float cape_thresh, float cinh_thresh) {
    if ((prof == NULL) || (elyr == NULL)) return;
    sharp::Profile* p = static_cast<sharp::Profile*>(prof->obj);
    sharp::PressureLayer* l = static_cast<sharp::PressureLayer*>(elyr->obj);
    sharp::PressureLayer out =
        sharp::effective_inflow_layer(p, cape_thresh, cinh_thresh);

    l->bottom = out.bottom;
    l->top = out.top;
}

void sharp_storm_motion_bunkers_np(sharp_Profile_t* prof,
                                   sharp_HeightLayer_t* mn_wind_lyr_agl,
                                   sharp_HeightLayer_t* wind_shr_lyr_agl,
                                   sharp_WindComponents_t* storm_motion,
                                   int leftMover, int pressureWeighted) {
    if ((prof == NULL) || (mn_wind_lyr_agl == NULL) ||
        (wind_shr_lyr_agl == NULL) || (storm_motion == NULL))
        return;
    sharp::Profile* p = static_cast<sharp::Profile*>(prof->obj);
    sharp::HeightLayer* mnl =
        static_cast<sharp::HeightLayer*>(mn_wind_lyr_agl->obj);
    sharp::HeightLayer* shl = 
        static_cast<sharp::HeightLayer*>(wind_shr_lyr_agl->obj);
    sharp::WindComponents* stm =
        static_cast<sharp::WindComponents*>(storm_motion->obj);

    sharp::WindComponents out = sharp::storm_motion_bunkers(
        p, *mnl, *shl, leftMover, pressureWeighted);
    stm->u = out.u;
    stm->v = out.v;
}

void sharp_storm_motion_bunkers(sharp_Profile_t* prof,
                                sharp_WindComponents_t* storm_motion,
                                int leftMover) {
    if ((prof == NULL) || (storm_motion == NULL)) return;
    sharp::Profile* p = static_cast<sharp::Profile*>(prof->obj);
    sharp::WindComponents* stm =
        static_cast<sharp::WindComponents*>(storm_motion->obj);
    sharp::WindComponents out = sharp::storm_motion_bunkers(p, leftMover);
    stm->u = out.u;
    stm->v = out.v;
}

float sharp_entrainment_cape(sharp_Profile_t* prof, sharp_Parcel_t* pcl) {
    if ((prof == NULL) || (pcl == NULL)) return sharp::MISSING;
    sharp::Profile* pf = static_cast<sharp::Profile*>(prof->obj);
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    return sharp::entrainment_cape(pf, pc);
}

float sharp_energy_helicity_index(float cape, float helicity) {
    return sharp::energy_helicity_index(cape, helicity);
}

float sharp_supercell_composite_parameter(float mu_cape, float eff_srh,
                                          float eff_shear) {
    return sharp::supercell_composite_parameter(mu_cape, eff_srh, eff_shear);
}

float sharp_significant_tornado_parameter(sharp_Profile_t* prof,
                                          sharp_Parcel_t* pcl, float helicity,
                                          float bulk_shear) {
    if ((prof == NULL) || (pcl == NULL)) return sharp::MISSING;
    sharp::Profile* pf = static_cast<sharp::Profile*>(prof->obj);
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    return sharp::significant_tornado_parameter(pf, *pc, helicity, bulk_shear);
}
