#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct {
  dbl min[3], max[3];
} rect3;

/**
 * Compute an approximate convex combination of three unit vectors. The rows
 * Are
 *
 * http://mathforum.org/library/drmath/view/65316.html
 * https://en.wikipedia.org/wiki/Slerp
 */
void lin_comb_unit_vec_3(dbl const t_in[3][3], dbl const b[3], dbl t_out[3]);

#ifdef __cplusplus
}
#endif
