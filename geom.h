#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct mesh3 mesh3_s;

typedef struct {
  dbl min[3], max[3];
} rect3;

/**
 * Check if four points are coplanar.
 */
bool points_are_coplanar(dbl const **x);

bool ray_and_face_are_coplanar(mesh3_s const *mesh, size_t l0, size_t l1,
                               size_t l2, dbl const *ray);

#ifdef __cplusplus
}
#endif
