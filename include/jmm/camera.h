#pragma once

#include "def.h"
#include "geom.h"

typedef enum camera_type {
  CAMERA_TYPE_UNINITIALIZED,
  CAMERA_TYPE_ORTHOGRAPHIC,
  CAMERA_TYPE_PERSPECTIVE
} camera_type_e;

typedef struct camera {
  camera_type_e type;
  dbl pos[3], look[3], left[3], up[3];
  dbl width, height; // parameters for orthographic camera
  dbl fovy, aspect; // parameters for perspective camera
  size_t dim[2];
} camera_s;

void camera_reset(camera_s *camera);
ray3 camera_get_ray_for_index(camera_s const *camera, int i, int j);
