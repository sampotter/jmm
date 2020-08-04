#include <stdio.h>
#include <stdlib.h>

#include "dial.h"

int main() {
  stype_e stype = CONSTANT;
  int n = 65;
  int nx = n;
  int ny = n;
  int nz = n;
  dbl h = 2.0/(nx - 1);
  int shape[3] = {nx, ny, nz};
  int ind0[3] = {nx, ny/2, 0};

  int m = 3*(n/4);
  int *ind = malloc(3*sizeof(int)*m*m*m);
  int l = 0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      for (int k = 0; k < m; ++k) {
        ind[3*l] = i;
        ind[3*l + 1] = j;
        ind[3*l + 2] = k;
        ++l;
      }
    }
  }

  dial3_s *dial;
  dial3_alloc(&dial);
  dial3_init(dial, stype, shape, h);
  dial3_add_boundary_points(dial, ind, m*m*m);
  dial3_add_point_source_with_trial_nbs(dial, ind0, 0);
  dial3_solve(dial);

  free(ind);
}
