#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

typedef struct {
  size_t data[2];
} ind2;

typedef struct {
  size_t data[3];
} ind3;

typedef struct {
  size_t data[4];
} ind4;

int ind2l(ivec2 shape, ivec2 ind);
int ind2lc(ivec2 shape, ivec2 ind);
int indc2l(ivec2 shape, ivec2 indc);
int indc2lc(ivec2 shape, ivec2 indc);
ivec2 l2ind(ivec2 shape, int l);
ivec2 l2indc(ivec2 shape, int l);
ivec2 lc2ind(ivec2 shape, int lc);
ivec2 lc2indc(ivec2 shape, int lc);
int l2lc(ivec2 shape, int l);
int lc2l(ivec2 shape, int lc);
int xy_to_lc_and_cc(ivec2 shape, dvec2 xymin, dbl h, dvec2 xy, dvec2 *cc);

int ind2l3(ivec3 shape, ivec3 ind);
ivec3 l2ind3(ivec3 shape, int l);

#ifdef __cplusplus
}
#endif
