#include <cgreen/cgreen.h>

#include "camera.h"

Describe(camera);

BeforeEach(camera) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(camera) {}

Ensure(camera, get_ray_works_for_orthographic_camera) {
  camera_s camera = {
    .type = CAMERA_TYPE_ORTHOGRAPHIC,
    .pos = {0, 0, -1},
    .look = {0, 0, 1},
    .left = {0, 1, 0},
    .up = {-1, 0, 0},
    .near = 1,
    .width = 2,
    .height = 1,
    .dim = {2, 4}
  };

  dbl x[2];
  ray3 ray;

  x[0] = 0;
  x[1] = 0;
  ray = camera_get_ray(&camera, x);
  assert_that_double(ray.org[0], is_nearly_double(-0.25));
  assert_that_double(ray.org[1], is_nearly_double(0.75));
  assert_that_double(ray.org[2], is_nearly_double(0));
  assert_that_double(ray.dir[0], is_nearly_double(0));
  assert_that_double(ray.dir[1], is_nearly_double(0));
  assert_that_double(ray.dir[2], is_nearly_double(1));
}
