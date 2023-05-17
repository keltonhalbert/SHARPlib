#include <SHARPlib/CWrap/winds_wrap.h>
#include <SHARPlib/winds.h>
#include <stdlib.h>

sharp_WindVector_t* sharp_WindVector_create() {
    sharp_WindVector_t* v;
    sharp::WindVector* vec;

    v = (sharp_WindVector_t*)malloc(sizeof(*v));
    vec = new sharp::WindVector();
    v->obj = vec;
    return v;
}

sharp_WindComponents_t* sharp_WindComponents_create() {
    sharp_WindComponents_t* c;
    sharp::WindComponents* comp;

    c = (sharp_WindComponents_t*)malloc(sizeof(*c));
    comp = new sharp::WindComponents();
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
