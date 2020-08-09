#include <stdio.h>
#include <stdlib.h>

#include "dial.h"
#include "index.h"

int main() {
  stype_e stype = CONSTANT;
  int n = 5;
  int nx = n;
  int ny = n;
  int nz = n;
  dbl h = 2.0/(nx - 1);
  int shape[3] = {nx, ny, nz};
  int ind0[3] = {n/2, n - 1, 0};

  int mx = 5, my = 3, mz = 3;

  int num_bd = mx*my*mz;
  int *ind = malloc(3*sizeof(int)*num_bd);
  int l = 0;
  for (int i = 0; i < mx; ++i) {
    for (int j = 0; j < my; ++j) {
      for (int k = 0; k < mz; ++k) {
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
  dial3_add_boundary_points(dial, ind, num_bd);
  dial3_add_point_source(dial, ind0, 0);
  
  printf("%d\n", dial3_get_state_ptr(dial)[33]);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        printf("%d ", dial3_get_state_ptr(dial)[ind2l3((ivec3) {.data = {nx, ny, nz}}, (ivec3) {.data = {i, j, k}})]);
      }
      printf("\n");
    }
    printf("\n");
  }
  
  dial3_solve(dial);
  
  printf("%d\n", dial3_get_state_ptr(dial)[33]);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        printf("%1.2f ", dial3_get_T(dial, ind2l3((ivec3) {.data = {nx, ny, nz}}, (ivec3) {.data = {i, j, k}})));
      }
      printf("\n");
    }
    printf("\n");
  }
  
  int l_ = ind2l3((ivec3) {.data = {nx, ny, nz}}, (ivec3) {.data = {1, 1, 3}});
  printf("T[%d] = %1.2f\n", l_, dial3_get_T(dial, l_));

  free(ind);
}
