#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "geom.h"

typedef enum camera_type {
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

ray3 camera_get_ray(camera_s const *camera, dbl x[2]);

#ifdef __cplusplus
}
#endif
