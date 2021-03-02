#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "geom.h"

typedef enum camera_type {
  CAMERA_TYPE_NONE,
  CAMERA_TYPE_ORTHOGRAPHIC,
  CAMERA_TYPE_PERSPECTIVE
} camera_type_e;

typedef struct camera {
  camera_type_e type;
  dbl pos[3], look[3], left[3], up[3];
  dbl width, height; // parameters for orthographic camera
  dbl near, fovy, aspect; // parameters for perspective camera
  size_t dim[2];
} camera_s;

void camera_init(camera_s *camera);
ray3 camera_get_ray_for_index(camera_s const *camera, int i, int j);

#ifdef __cplusplus
}
#endif
