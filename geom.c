#include "geom.h"

#include <math.h>

#include "mat.h"
#include "mesh3.h"

void lin_comb_unit_vec_3(dbl const t_in[3][3], dbl const b[3], dbl t_out[3]) {
  for (int i = 0; i < 3; ++i) {
    t_out[i] = 0;
    for (int j = 0; j < 3; ++j) {
      t_out[i] += b[i]*t_in[i][j];
    }
  }
  dbl3_normalize(t_out);

bool points_are_coplanar(dbl const **x) {
  dbl dx[3][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);
  dbl3_sub(x[3], x[2], dx[2]);
  dbl det = dbl33_det(dx);
  return fabs(det) < 1e-15;
}

bool ray_and_face_are_coplanar(mesh3_s const *mesh,
                               size_t l0, size_t l1, size_t l2, dbl const *ray) {
  dbl const *x[3] = {
    mesh3_get_vert_ptr(mesh, l0),
    mesh3_get_vert_ptr(mesh, l1),
    mesh3_get_vert_ptr(mesh, l2)
  };

  dbl dx[2][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);

  dbl n[3];
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);

  return fabs(dbl3_dot(n, ray)) < 1e-14;
}
