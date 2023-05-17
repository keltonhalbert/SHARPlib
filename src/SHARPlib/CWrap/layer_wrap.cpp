#include <SHARPlib/CWrap/layer_wrap.h>
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
    sharp::PressureLayer* lyr;
    lyr = static_cast<sharp::PressureLayer*>(l->obj);
    return lyr->bottom;
}

float sharp_PressureLayer_get_top(sharp_PressureLayer_t* l) {
    sharp::PressureLayer* lyr;
    lyr = static_cast<sharp::PressureLayer*>(l->obj);
    return lyr->top;
}

float sharp_HeightLayer_get_bottom(sharp_HeightLayer_t* l) {
    sharp::HeightLayer* lyr;
    lyr = static_cast<sharp::HeightLayer*>(l->obj);
    return lyr->bottom;
}

float sharp_HeightLayer_get_top(sharp_HeightLayer_t* l) {
    sharp::HeightLayer* lyr;
    lyr = static_cast<sharp::HeightLayer*>(l->obj);
    return lyr->top;
}
