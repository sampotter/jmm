#include <cgreen/cgreen.h>

#include "camera.h"

Describe(camera);

BeforeEach(camera) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(camera) {}

Ensure(camera, get_ray_for_index_works_for_orthographic_camera) {
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

  ray3 ray;

  ray = camera_get_ray_for_index(&camera, 0, 0);
  assert_that_double(ray.org[0], is_nearly_double(-0.25));
  assert_that_double(ray.org[1], is_nearly_double(0.75));
  assert_that_double(ray.org[2], is_nearly_double(-1));
  assert_that_double(ray.dir[0], is_nearly_double(0));
  assert_that_double(ray.dir[1], is_nearly_double(0));
  assert_that_double(ray.dir[2], is_nearly_double(1));
}

Ensure(camera, get_ray_for_index_works_for_perspective_camera) {
  camera_s camera = {
    .type = CAMERA_TYPE_PERSPECTIVE,
    .pos = {0, -1, 0},
    .look = {0, 1, 0},
    .left = {-1, 0, 0},
    .up = {0, 0, 1},
    .near = 1,
    .fovy = 90,
    .aspect = 2,
    .dim = {2, 4}
  };

  ray3 ray;

  ray = camera_get_ray_for_index(&camera, 0, 0);
  assert_that_double(ray.org[0], is_nearly_double(0));
  assert_that_double(ray.org[1], is_nearly_double(-1));
  assert_that_double(ray.org[2], is_nearly_double(0));
  assert_that_double(ray.dir[0], is_nearly_double(-0.8017837257372732));
  assert_that_double(ray.dir[1], is_nearly_double(0.5345224838248488));
  assert_that_double(ray.dir[2], is_nearly_double(0.2672612419124244));
}
