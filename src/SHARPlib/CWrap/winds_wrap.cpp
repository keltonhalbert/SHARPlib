#include <SHARPlib/CWrap/winds_wrap.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/winds.h>
#include <stdlib.h>

sharp_WindVector_t* sharp_WindVector_create(float speed, float direction) {
    sharp_WindVector_t* v;
    sharp::WindVector* vec;

    v = (sharp_WindVector_t*)malloc(sizeof(*v));
    vec = new sharp::WindVector();
    vec->speed = speed;
    vec->direction = direction;
    v->obj = vec;
    return v;
}

sharp_WindComponents_t* sharp_WindComponents_create(float u, float v) {
    sharp_WindComponents_t* c;
    sharp::WindComponents* comp;

    c = (sharp_WindComponents_t*)malloc(sizeof(*c));
    comp = new sharp::WindComponents();
    comp->u = u;
    comp->v = v;
    c->obj = comp;
    return c;
}

void sharp_WindVector_delete(sharp_WindVector_t* v) {
    if (v == NULL) return;
    delete static_cast<sharp::WindVector*>(v->obj);
    free(v);
}

void sharp_WindComponents_delete(sharp_WindComponents_t* c) {
    if (c == NULL) return;
    delete static_cast<sharp::WindComponents*>(c->obj);
    free(c);
}

float sharp_WindVector_get_speed(sharp_WindVector_t* v) {
    if (v == NULL) return sharp::MISSING;
    return static_cast<sharp::WindVector*>(v->obj)->speed;
}

float sharp_WindVector_get_direction(sharp_WindVector_t* v) {
    if (v == NULL) return sharp::MISSING;
    return static_cast<sharp::WindVector*>(v->obj)->direction;
}

float sharp_WindComponents_get_u(sharp_WindComponents_t* c) {
    if (c == NULL) return sharp::MISSING;
    return static_cast<sharp::WindComponents*>(c->obj)->u;
}

float sharp_WindComponents_get_v(sharp_WindComponents_t* c) {
    if (c == NULL) return sharp::MISSING;
    return static_cast<sharp::WindComponents*>(c->obj)->v;
}

float sharp_u_component(float wind_speed, float wind_direction) {
    return sharp::u_component(wind_speed, wind_direction);
}

float sharp_v_component(float wind_speed, float wind_direction) {
    return sharp::v_component(wind_speed, wind_direction);
}

float sharp_vector_angle(float u_comp, float v_comp) {
    return sharp::vector_angle(u_comp, v_comp);
}

float sharp_vector_magnitude(float u_comp, float v_comp) {
    return sharp::vector_magnitude(u_comp, v_comp);
}

void sharp_components_to_vector(sharp_WindComponents_t* comp,
                                sharp_WindVector_t* vec) {
    if ((comp == NULL) || (vec == NULL)) return;
    sharp::WindComponents* c = static_cast<sharp::WindComponents*>(comp->obj);
    sharp::WindVector* v = static_cast<sharp::WindVector*>(vec->obj);
    sharp::WindVector out = sharp::components_to_vector(*c);
    v->speed = out.speed;
    v->direction = out.direction;
}

void sharp_vector_to_components(sharp_WindVector_t* vec,
                                sharp_WindComponents_t* comp) {
    if ((vec == NULL) || (comp == NULL)) return;
    sharp::WindVector* v = static_cast<sharp::WindVector*>(vec->obj);
    sharp::WindComponents* c = static_cast<sharp::WindComponents*>(comp->obj);
    sharp::WindComponents out = sharp::vector_to_components(*v);
    c->u = out.u;
    c->v = out.v;
}

void sharp_mean_wind(sharp_PressureLayer_t* plyr, sharp_WindComponents_t* cmp,
                     const float* pres, const float* u_wind,
                     const float* v_wind, const int N, const int weighted) {
    if ((plyr == NULL) || (cmp == NULL)) return;
    sharp::PressureLayer* p = static_cast<sharp::PressureLayer*>(plyr->obj);
    sharp::WindComponents* c = static_cast<sharp::WindComponents*>(cmp->obj);
    sharp::WindComponents out =
        sharp::mean_wind(*p, pres, u_wind, v_wind, N, weighted);
    c->u = out.u;
    c->v = out.v;
}

void sharp_PressureLayer_wind_shear(sharp_PressureLayer_t* plyr,
                                    sharp_WindComponents_t* cmp,
                                    const float* pres, const float* u_wind,
                                    const float* v_wind, const int N) {
    if ((plyr == NULL) || (cmp == NULL)) return;
    sharp::PressureLayer* p = static_cast<sharp::PressureLayer*>(plyr->obj);
    sharp::WindComponents* c = static_cast<sharp::WindComponents*>(cmp->obj);
    sharp::WindComponents out = sharp::wind_shear(*p, pres, u_wind, v_wind, N);
    c->u = out.u;
    c->v = out.v;
}

void sharp_HeightLayer_wind_shear(sharp_HeightLayer_t* hlyr,
                                  sharp_WindComponents_t* cmp,
                                  const float* hght, const float* u_wind,
                                  const float* v_wind, const int N) {
    if ((hlyr == NULL) || (cmp == NULL)) return;
    sharp::HeightLayer* h = static_cast<sharp::HeightLayer*>(hlyr->obj);
    sharp::WindComponents* c = static_cast<sharp::WindComponents*>(cmp->obj);
    sharp::WindComponents out = sharp::wind_shear(*h, hght, u_wind, v_wind, N);
    c->u = out.u;
    c->v = out.v;
}

float sharp_HeightLayer_helicity(sharp_HeightLayer_t* hlyr,
                                 sharp_WindComponents_t* storm_motion,
                                 const float* height, const float* u_wind,
                                 const float* v_wind, const int N) {
    if ((hlyr == NULL) || (storm_motion == NULL)) return sharp::MISSING;
    sharp::HeightLayer* h = static_cast<sharp::HeightLayer*>(hlyr->obj);
    sharp::WindComponents* stm =
        static_cast<sharp::WindComponents*>(storm_motion->obj);
    return sharp::helicity(*h, *stm, height, u_wind, v_wind, N);
}

float sharp_PressureLayer_helicity(sharp_PressureLayer_t* plyr,
                                   sharp_WindComponents_t* storm_motion,
                                   const float* pressure, const float* u_wind,
                                   const float* v_wind, const int N) {
    if ((plyr == NULL) || (storm_motion == NULL)) return sharp::MISSING;
    sharp::PressureLayer* p = static_cast<sharp::PressureLayer*>(plyr->obj);
    sharp::WindComponents* stm =
        static_cast<sharp::WindComponents*>(storm_motion->obj);
    return sharp::helicity(*p, *stm, pressure, u_wind, v_wind, N);
}
