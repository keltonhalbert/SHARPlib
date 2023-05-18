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

