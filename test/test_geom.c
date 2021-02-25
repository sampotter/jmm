#include <cgreen/cgreen.h>

#include "geom.h"

Describe(geom);

BeforeEach(geom) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(geom) {}

Ensure(geom, ray3_intersects_tetra3_works) {
  ray3 ray;
  tetra3 tetra;
  bool hit;
  dbl t;

  tetra = (tetra3) {
    .v = {
      {0, 0, 0},
      {1, 0, 0},
      {0, 1, 0},
      {0, 0, 1}
    }
  };

  ray = (ray3) {
    .org = {1, 1, 1},
    .dir = {-1/SQRT3, -1/SQRT3, -1/SQRT3}
  };
  hit = ray3_intersects_tetra3(&ray, &tetra, &t);
  assert_that(hit);
  assert_that_double(t, is_nearly_double(SQRT3 - 1/SQRT3));

  ray = (ray3) {
    .org = {0, 0, 0},
    .dir = {1, 0, 0}
  };
  hit = ray3_intersects_tetra3(&ray, &tetra, &t);
  assert_that(hit);
  assert_that_double(t, is_nearly_double(0));

  ray = (ray3) {
    .org = {2, 0, 0},
    .dir = {-1, 0, 0}
  };
  hit = ray3_intersects_tetra3(&ray, &tetra, &t);
  assert_that(hit);
  assert_that_double(t, is_nearly_double(1));

  ray = (ray3) {
    .org = {0, 2, 0},
    .dir = {0, -1, 0}
  };
  hit = ray3_intersects_tetra3(&ray, &tetra, &t);
  assert_that(hit);
  assert_that_double(t, is_nearly_double(1));

  ray = (ray3) {
    .org = {0, 0, 2},
    .dir = {0, 0, -1}
  };
  hit = ray3_intersects_tetra3(&ray, &tetra, &t);
  assert_that(hit);
  assert_that_double(t, is_nearly_double(1));
}

Ensure(geom, ray3_intersects_tri3_works) {
  ray3 ray;
  tri3 tri;
  bool hit;
  dbl t;

  tri = (tri3) {
    .v = {
      {0, 0, 0},
      {1, 0, 0},
      {0, 1, 0}
    }
  };

  ray = (ray3) {
    .org = {0, 0, 1},
    .dir = {0, 0, -1}
  };
  hit = ray3_intersects_tri3(&ray, &tri, &t);
  assert_that(hit);
  assert_that_double(t, is_nearly_double(1));

  ray = (ray3) {
    .org = {1./3, 1./3, 0},
    .dir = {0, 0, 1},
  };
  hit = ray3_intersects_tri3(&ray, &tri, &t);
  assert_that(hit);
  assert_that_double(t, is_nearly_double(0));
}
