#pragma once

struct sjs;

#include "hermite.h"
#include "vec.h"

typedef dbl (*sfield)(dvec2);
typedef dvec2 (*vfield)(dvec2);

typedef struct sjs sjs_s;

void sjs_alloc(sjs_s **sjs);
void sjs_dealloc(sjs_s **sjs);
void sjs_init(sjs_s *sjs, ivec2 shape, dbl h, sfield s, vfield grad_s);
void sjs_deinit(sjs_s *sjs);
void sjs_add_fac_pt_src(sjs_s *sjs, ivec2 ind, dbl r, int *nf, int *nfc);
void sjs_solve(sjs_s *sjs);
dbl sjs_T(sjs_s *sjs, dvec2 xy);
dbl sjs_Tx(sjs_s *sjs, dvec2 xy);
dbl sjs_Ty(sjs_s *sjs, dvec2 xy);
dbl sjs_Txy(sjs_s *sjs, dvec2 xy);
bicubic *sjs_bicubic(sjs_s *sjs, ivec2 cind);
