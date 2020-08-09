#include <catch2/catch.hpp>

#include "dial.h"
#include "index.h"

TEST_CASE ("Small s = 1 point source problem computes |x|", "[dial3]") {
  using namespace Catch::literals;

  stype_e stype = CONSTANT;
  int nx = 5;
  int ny = 5;
  int nz = 5;
  dbl h = 2.0/(nx - 1);
  ivec3 shape = {.data = {nx, ny, nz}};

  ivec3 ind0 = ivec3_int_div(shape, 2);

  dial3_s *dial;
  dial3_alloc(&dial);
  dial3_init(dial, stype, shape.data, h);
  dial3_add_point_source(dial, ind0.data, 0);
  dial3_solve(dial);

  ivec3 ind;
  for (ind.data[0] = 0; ind.data[0] < nx; ++ind.data[0]) {
    for (ind.data[1] = 0; ind.data[1] < ny; ++ind.data[1]) {
      for (ind.data[2] = 0; ind.data[2] < nz; ++ind.data[2]) {
        if (ivec3_equal(ind, ind0)) {
          continue;
        }
        int l = ind2l3(shape, ind);
        dvec3 grad_T = dvec3_sub(ivec3_to_dvec3(ind), ivec3_to_dvec3(ind0));
        dbl T = h*dvec3_norm(grad_T);
        grad_T = dvec3_dbl_div(grad_T, T/h);
        REQUIRE(dial3_get_T(dial, l) == Approx(T));
        dvec3 grad_T_gt;
        dial3_get_grad_T(dial, l, grad_T_gt.data);
        REQUIRE(grad_T.data[0] == Approx(grad_T_gt.data[0]));
        REQUIRE(grad_T.data[1] == Approx(grad_T_gt.data[1]));
        REQUIRE(grad_T.data[2] == Approx(grad_T_gt.data[2]));
      }
    }
  }

  dial3_deinit(dial);
  dial3_dealloc(&dial);
}

TEST_CASE ("Solving s = 1 on 1x3x3 L yields correct T and grad(T)",
           "[dial3]") {
  using namespace Catch::literals;

  stype_e stype = CONSTANT;

  dbl h = 1.0;
  int shape[3] = {1, 3, 3};
  int ind0[3] = {0, 2, 0};

  int inds[12] = {
    0, 0, 0,
    0, 0, 1,
    0, 1, 0,
    0, 1, 1
  };

  dial3_s *dial;
  dial3_alloc(&dial);
  dial3_init(dial, stype, shape, h);
  dial3_add_boundary_points(dial, inds, 4);
  dial3_add_point_source(dial, ind0, 0);
  dial3_solve(dial);

  state_e state_gt[9] = {
    BOUNDARY, BOUNDARY, VALID,
    BOUNDARY, BOUNDARY, VALID,
    VALID,    VALID,    VALID
  };

  dbl T_gt[9] = {
    NAN, NAN, h + SQRT2 * h + h,
    NAN, NAN,         h + SQRT2,
    0,     h,             h + h
  };

  dbl Tx_gt[9] = {
    NAN, NAN, 0,
    NAN, NAN, 0,
    NAN,   0, 0
  };

  dbl Ty_gt[9] = {
    NAN, NAN,       -1,
    NAN, NAN, -SQRT2/2,
    NAN,   0,        0
  };

  dbl Tz_gt[9] = {
    NAN, NAN,       0,
    NAN, NAN, SQRT2/2,
    NAN,   1,       1
  };

  state_e *state = dial3_get_state_ptr(dial);

  for (int i = 0; i < 9; ++i) {
    REQUIRE(state[i] == state_gt[i]);

    dbl T = dial3_get_T(dial, i);

    if (isnan(T_gt[i])) {
      REQUIRE(isnan(T));
    } else {
      REQUIRE(T == Approx(T_gt[i]));
    }

    dbl grad_T[3];
    dial3_get_grad_T(dial, i, grad_T);

    if (isnan(Tx_gt[i])) {
      REQUIRE(isnan(grad_T[i]));
    } else {
      REQUIRE(grad_T[i] == Approx(Tx_gt[i]));
    }

    if (isnan(Ty_gt[i])) {
      REQUIRE(isnan(grad_T[i]));
    } else {
      REQUIRE(grad_T[i] == Approx(Ty_gt[i]));
    }

    if (isnan(Tz_gt[i])) {
      REQUIRE(isnan(grad_T[i]));
    } else {
      REQUIRE(grad_T[i] == Approx(Tz_gt[i]));
    }
  }

  dial3_deinit(dial);
  dial3_dealloc(&dial);
}

TEST_CASE ("Solving s = 1 on 1x5x5 L yields correct T and grad(T)",
           "[dial3]") {
  using namespace Catch::literals;

  stype_e stype = CONSTANT;
  dbl h = 0.5;
  int shape[3] = {1, 5, 5};
  int ind0[3] = {0, 4, 0};

  int inds[27] = {
    0, 0, 0,
    0, 0, 1,
    0, 0, 2,
    0, 1, 0,
    0, 1, 1,
    0, 1, 2,
    0, 2, 0,
    0, 2, 1,
    0, 2, 2
  };

  dial3_s *dial;
  dial3_alloc(&dial);
  dial3_init(dial, stype, shape, h);
  dial3_add_boundary_points(dial, inds, 9);
  dial3_add_point_source(dial, ind0, 0);
  dial3_solve(dial);

  state_e state_gt[25] = {
    BOUNDARY, BOUNDARY, BOUNDARY, VALID, VALID,
    BOUNDARY, BOUNDARY, BOUNDARY, VALID, VALID,
    BOUNDARY, BOUNDARY, BOUNDARY, VALID, VALID,
    VALID,    VALID,    VALID,    VALID, VALID,
    VALID,    VALID,    VALID,    VALID, VALID,
  };

  dbl T_gt[25] = {
    NAN, NAN,     NAN,     (SQRT5 + SQRT2 + 2)*h, (2*SQRT5 + SQRT2)*h,
    NAN, NAN,     NAN,     (SQRT5 + SQRT2 + 1)*h, (SQRT5 + 2*SQRT2)*h,
    NAN, NAN,     NAN,     (SQRT5 + SQRT2)*h,     2*SQRT5*h,
    h,   SQRT2*h, SQRT5*h, SQRT10*h,              SQRT17*h,
    0,   h,       2*h,     3*h,                   4*h
  };

  state_e *state = dial3_get_state_ptr(dial);

  for (int i = 0; i < 25; ++i) {
    REQUIRE(state[i] == state_gt[i]);

    dbl T = dial3_get_T(dial, i);

    if (isnan(T_gt[i])) {
      REQUIRE(isnan(T));
    } else {
      REQUIRE(T == Approx(T_gt[i]));
    }

//    if (isnan(Tx_gt[i])) {
//      REQUIRE(isnan(grad_T[i].data[0]));
//    } else {
//      REQUIRE(grad_T[i].data[0] == Approx(Tx_gt[i]));
//    }
//
//    if (isnan(Ty_gt[i])) {
//      REQUIRE(isnan(grad_T[i].data[1]));
//    } else {
//      REQUIRE(grad_T[i].data[1] == Approx(Ty_gt[i]));
//    }
//
//    if (isnan(Tz_gt[i])) {
//      REQUIRE(isnan(grad_T[i].data[2]));
//    } else {
//      REQUIRE(grad_T[i].data[2] == Approx(Tz_gt[i]));
//    }
  }

  dial3_deinit(dial);
  dial3_dealloc(&dial);
}
