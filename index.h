#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

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

#ifdef __cplusplus
}
#endif
