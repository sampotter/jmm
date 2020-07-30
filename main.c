#include <stdio.h>

#include "dial.h"

int main() {
  stype_e stype = CONSTANT;
  int nx = 111;
  int ny = 111;
  int nz = 111;
  dbl h = 2.0/(nx - 1);
  ivec3 shape = {.data = {nx, ny, nz}};

  ivec3 ind0 = ivec3_int_div(shape, 2);

  dial3_s *dial;
  dial3_alloc(&dial);
  dial3_init(dial, stype, shape, h);
  dial3_add_point_source_with_trial_nbs(dial, ind0, 0);
  dial3_solve(dial);
}
