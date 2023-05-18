#ifndef __SHARP_WINDS_WRAP_H__
#define __SHARP_WINDS_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif


struct sharp_WindVector {
    void *obj;
};

struct sharp_WindComponents {
    void *obj;
};

typedef struct sharp_WindVector sharp_WindVector_t;
typedef struct sharp_WindComponents sharp_WindComponents_t;

sharp_WindVector_t* sharp_WindVector_create(float speed, float direction);
sharp_WindComponents_t* sharp_WindComponents_create(float u, float v);

void sharp_WindVector_delete(sharp_WindVector_t* v);
void sharp_WindComponents_delete(sharp_WindComponents_t* c);

float sharp_WindVector_get_speed(sharp_WindVector_t* v);
float sharp_WindVector_get_direction(sharp_WindVector_t* v);

float sharp_WindComponents_get_u(sharp_WindComponents_t* c);
float sharp_WindComponents_get_v(sharp_WindComponents_t* c);

#ifdef __cplusplus
}
#endif

#endif // __SHARP_WINDS_WRAP_H__
