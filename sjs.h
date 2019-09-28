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
int sjs_add_fac_pt_src(sjs_s *sjs, ivec2 ind, dbl r);
void sjs_solve(sjs_s *sjs);
bicubic *sjs_bicubic(sjs_s *sjs, ivec2 ind);
