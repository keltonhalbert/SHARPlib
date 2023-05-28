#include <SHARPlib/CWrap/layer_wrap.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <stdlib.h>

sharp_PressureLayer_t* sharp_PressureLayer_create(float bottom, float top) {
    sharp_PressureLayer_t* l;
    sharp::PressureLayer* lyr;

    l = (sharp_PressureLayer_t*)malloc(sizeof(*l));
    lyr = new sharp::PressureLayer(bottom, top);
    l->obj = lyr;
    return l;
}

sharp_HeightLayer_t* sharp_HeightLayer_create(float bottom, float top) {
    sharp_HeightLayer_t* l;
    sharp::HeightLayer* lyr;

    l = (sharp_HeightLayer_t*)malloc(sizeof(*l));
    lyr = new sharp::HeightLayer(bottom, top);
    l->obj = lyr;
    return l;
}

sharp_LayerIndex_t* sharp_LayerIndex_create() {
    sharp_LayerIndex_t* i;
    sharp::LayerIndex* idx;

    i = (sharp_LayerIndex_t*)malloc(sizeof(*i));
    idx = new sharp::LayerIndex();
    i->obj = idx;
    return i;
}

void sharp_PressureLayer_delete(sharp_PressureLayer_t* l) {
    if (l == NULL) return;
    delete static_cast<sharp::PressureLayer*>(l->obj);
    free(l);
}

void sharp_HeightLayer_delete(sharp_HeightLayer_t* l) {
    if (l == NULL) return;
    delete static_cast<sharp::HeightLayer*>(l->obj);
    free(l);
}

void sharp_LayerIndex_delete(sharp_LayerIndex_t* i) {
    if (i == NULL) return;
    delete static_cast<sharp::LayerIndex*>(i->obj);
    free(i);
}

float sharp_PressureLayer_get_bottom(sharp_PressureLayer_t* l) {
    if (l == NULL) return sharp::MISSING;
    return static_cast<sharp::PressureLayer*>(l->obj)->bottom;
}

float sharp_PressureLayer_get_top(sharp_PressureLayer_t* l) {
    if (l == NULL) return sharp::MISSING;
    return static_cast<sharp::PressureLayer*>(l->obj)->top;
}

float sharp_HeightLayer_get_bottom(sharp_HeightLayer_t* l) {
    if (l == NULL) return sharp::MISSING;
    return static_cast<sharp::HeightLayer*>(l->obj)->bottom;
}

float sharp_HeightLayer_get_top(sharp_HeightLayer_t* l) {
    if (l == NULL) return sharp::MISSING;
    return static_cast<sharp::HeightLayer*>(l->obj)->top;
}

int sharp_LayerIndex_get_bottom(sharp_LayerIndex_t* i) {
    if (i == NULL) return 0;
    return static_cast<sharp::LayerIndex*>(i->obj)->kbot;
}

int sharp_LayerIndex_get_top(sharp_LayerIndex_t* i) {
    if (i == NULL) return 0;
    return static_cast<sharp::LayerIndex*>(i->obj)->ktop;
}

void sharp_get_PressureLayer_index(sharp_PressureLayer_t* lyr,
                                   sharp_LayerIndex_t* idx,
                                   const float* pressure, const int NZ) {
    if ((lyr == NULL) || (idx == NULL)) return;
    sharp::PressureLayer* pl = static_cast<sharp::PressureLayer*>(lyr->obj);
    sharp::LayerIndex* i = static_cast<sharp::LayerIndex*>(idx->obj);
    sharp::LayerIndex out = sharp::get_layer_index(*pl, pressure, NZ);
    i->kbot = out.kbot;
    i->ktop = out.ktop;
}

void sharp_get_HeightLayer_index(sharp_HeightLayer_t* lyr,
                                   sharp_LayerIndex_t* idx,
                                   const float* height, const int NZ) {
    if ((lyr == NULL) || (idx == NULL)) return;
    sharp::HeightLayer* hl = static_cast<sharp::HeightLayer*>(lyr->obj);
    sharp::LayerIndex* i = static_cast<sharp::LayerIndex*>(idx->obj);
    sharp::LayerIndex out = sharp::get_layer_index(*hl, height, NZ);
    i->kbot = out.kbot;
    i->ktop = out.ktop;
}


void sharp_HeightLayer_to_PressureLayer(sharp_HeightLayer_t* hlyr,
                                        sharp_PressureLayer_t* plyr,
                                        const float* pressure,
                                        const float* height, const int NZ,
                                        bool isAGL) {
    if ((hlyr == NULL) || (plyr == NULL)) return;
    sharp::HeightLayer* hl = static_cast<sharp::HeightLayer*>(hlyr->obj);
    sharp::PressureLayer* pl = static_cast<sharp::PressureLayer*>(plyr->obj);

    sharp::PressureLayer out =
        sharp::height_layer_to_pressure(*hl, pressure, height, NZ, isAGL);
    hl->bottom = out.bottom;
    hl->top = out.top;
}

void sharp_PressureLayer_to_HeightLayer(sharp_PressureLayer_t* plyr,
                                        sharp_HeightLayer_t* hlyr,
                                        const float* pressure,
                                        const float* height, const int NZ,
                                        bool toAGL) {
    if ((hlyr == NULL) || (plyr == NULL)) return;
    sharp::PressureLayer* pl = static_cast<sharp::PressureLayer*>(plyr->obj);
    sharp::HeightLayer* hl = static_cast<sharp::HeightLayer*>(hlyr->obj);

    sharp::HeightLayer out =
        sharp::pressure_layer_to_height(*pl, pressure, height, NZ, toAGL);
    pl->bottom = out.bottom;
    pl->top = out.top;
}

float sharp_PressureLayer_min(sharp_PressureLayer_t* plyr, const float* pres,
                              const float* data, const int NZ, float* lvl_of_min) {
    if (plyr == NULL) return sharp::MISSING;
    sharp::PressureLayer* pl = static_cast<sharp::PressureLayer*>(plyr->obj);
    return sharp::layer_min(*pl, pres, data, NZ, lvl_of_min);
}

float sharp_PressureLayer_max(sharp_PressureLayer_t* plyr, const float* pres,
                              const float* data, const int NZ, float* lvl_of_max) {
    if (plyr == NULL) return sharp::MISSING;
    sharp::PressureLayer* pl = static_cast<sharp::PressureLayer*>(plyr->obj);
    return sharp::layer_max(*pl, pres, data, NZ, lvl_of_max);
}

float sharp_HeightLayer_min(sharp_HeightLayer_t* hlyr, const float* hght,
                            const float* data, const int NZ, float* lvl_of_min) {
    if (hlyr == NULL) return sharp::MISSING;
    sharp::HeightLayer* hl = static_cast<sharp::HeightLayer*>(hlyr->obj);
    return sharp::layer_min(*hl, hght, data, NZ, lvl_of_min);
}

float sharp_HeightLayer_max(sharp_HeightLayer_t* hlyr, const float* hght,
                            const float* data, const int NZ, float* lvl_of_max) {
    if (hlyr == NULL) return sharp::MISSING;
    sharp::HeightLayer* hl = static_cast<sharp::HeightLayer*>(hlyr->obj);
    return sharp::layer_max(*hl, hght, data, NZ, lvl_of_max);
}

float sharp_PressureLayer_mean(sharp_PressureLayer_t* plyr,
                               const float* pressure, const float* data,
                               const int NZ) {
    sharp::PressureLayer* pl = static_cast<sharp::PressureLayer*>(plyr->obj);
    return sharp::layer_mean(*pl, pressure, data, NZ);
}

float sharp_HeightLayer_mean(sharp_HeightLayer_t* hlyr, const float* height,
                             const float* pressure, const float* data, const int NZ,
                             int isAGL) {
    sharp::HeightLayer* hl = static_cast<sharp::HeightLayer*>(hlyr->obj);
    return sharp::layer_mean(*hl, height, pressure, data, NZ, isAGL);
}

float sharp_PressureLayer_integrate(sharp_PressureLayer_t* plyr,
                                    const float* data, const float* pressure,
                                    const int NZ, int integ_sign, int weighted) {
    sharp::PressureLayer* pl = static_cast<sharp::PressureLayer*>(plyr->obj);
    return sharp::integrate_layer_trapz(*pl, data, pressure, NZ, integ_sign,
                                        weighted);
}

float sharp_HeightLayer_integrate(sharp_HeightLayer_t* hlyr, const float* data,
                                  const float* height, const int NZ, int integ_sign,
                                  int weighted) {
    sharp::HeightLayer* hl = static_cast<sharp::HeightLayer*>(hlyr->obj);
    return sharp::integrate_layer_trapz(*hl, data, height, NZ, integ_sign,
                                        weighted);
}
