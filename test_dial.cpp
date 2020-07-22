#include <catch2/catch.hpp>

#include "dial.h"
#include "index.h"

TEST_CASE ("Small s = 1 point source problem computes |x|", "[dial3]") {
  using namespace Catch::literals;

  stype_e stype = CONSTANT;
  int nx = 11;
  int ny = 11;
  int nz = 11;
  dbl h = 2.0/(nx - 1);
  ivec3 shape = {.data = {nx, ny, nz}};

  ivec3 ind0 = ivec3_int_div(shape, 2);

  dial3_s *dial;
  dial3_alloc(&dial);
  dial3_init(dial, stype, shape, h);
  dial3_add_point_source_with_trial_nbs(dial, ind0, 0);
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
        dvec3 grad_T_gt = dial3_get_grad_T(dial, l);
        REQUIRE(grad_T.data[0] == Approx(grad_T_gt.data[0]));
        REQUIRE(grad_T.data[1] == Approx(grad_T_gt.data[1]));
        REQUIRE(grad_T.data[2] == Approx(grad_T_gt.data[2]));
      }
    }
  }

  dial3_deinit(dial);
  dial3_dealloc(&dial);
}
