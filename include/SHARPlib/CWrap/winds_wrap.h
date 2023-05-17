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

sharp_WindVector_t* sharp_WindVector_create();
sharp_WindComponents_t* sharp_WindComponents_create();

void sharp_WindVector_delete(sharp_WindVector_t* v);
void sharp_WindComponents_delete(sharp_WindComponents_t* c);

#ifdef __cplusplus
}
#endif

#endif // __SHARP_WINDS_WRAP_H__
