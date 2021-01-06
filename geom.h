#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct {
  dbl min[3], max[3];
} rect3;

void lin_comb_unit_vec_3(dbl const t_in[3][3], dbl const b[3], dbl t_out[3]);

#ifdef __cplusplus
}
#endif
