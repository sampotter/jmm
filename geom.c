#include "geom.h"

#include <math.h>

#include "mat.h"
#include "mesh3.h"

bool points_are_coplanar(dbl const **x) {
  dbl dx[3][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);
  dbl3_sub(x[3], x[2], dx[2]);
  dbl det = dbl33_det(dx);
  return fabs(det) < 1e-15;
}

bool ray_and_face_are_coplanar(mesh3_s const *mesh, size_t l0, size_t l1,
                               size_t l2, dbl const *ray) {
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

dbl tri_area(dbl const x[3], dbl const y[3], dbl const z[3]) {
  dbl xy[3], xz[3];
  dbl3_sub(y, x, xy);
  dbl3_sub(z, x, xz);
  dbl cp[3];
  dbl3_cross(xy, xz, cp);
  return dbl3_norm(cp)/2;
}

void get_bary_coords_3d(dbl const *x[3], dbl const y[3], dbl b[3]) {
  dbl dx[2][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);

  dbl n[3];
  dbl3_cross(dx[0], dx[1], n);
  dbl area = dbl3_normalize(n)/2;

  b[0] = tri_area(y, x[1], x[2])/area;
  b[1] = tri_area(x[0], y, x[2])/area;
  b[2] = tri_area(x[0], x[1], y)/area;
}

bool adj_tetra_pair_is_convex(mesh3_s const *mesh, size_t l0,
                              size_t const lf[3], size_t l1) {
  // Get vertex data from mesh
  // TODO: probably faster to just get pointers
  dbl x0[3], x1[3], x[3][3];
  mesh3_copy_vert(mesh, l0, x0);
  mesh3_copy_vert(mesh, l1, x1);
  for (int i = 0; i < 3; ++i)
    mesh3_copy_vert(mesh, lf[i], x[i]);

  dbl dx[2][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);

  // Compute normal vector for plane spanned by (x1 - x0, x2 - x0).
  dbl n[3];
  dbl3_cross(dx[0], dx[1], n);
  dbl area = dbl3_normalize(n)/2;

  // Find intersection between line [x0, x1] and plane p(x) = n'*(x - x0).
  dbl d[3];
  dbl3_sub(x1, x0, d);
  dbl s = (dbl3_dot(n, x[0]) - dbl3_dot(n, x0))/dbl3_dot(n, d);
  dbl xs[3];
  dbl3_saxpy(s, x0, d, xs);

  // Compute barycentric coordinates of xs
  dbl const b[3] = {
    tri_area(xs, x[1], x[2])/area,
    tri_area(x[0], xs, x[2])/area,
    tri_area(x[0], x[1], xs)/area
  };

  return b[0] > 0 && b[1] > 0 && b[2] > 0;
}

dbl min_tetra_altitude(dbl const x[4][3]) {
  dbl n[3], dx[3][3];

  dbl h = INFINITY, newh;

  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);
  dbl3_sub(x[3], x[0], dx[2]);
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
  newh = fabs(dbl3_dot(n, dx[2]));
  h = fmin(h, newh);

  dbl3_sub(x[0], x[1], dx[0]);
  dbl3_sub(x[2], x[1], dx[1]);
  dbl3_sub(x[3], x[1], dx[2]);
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
  newh = fabs(dbl3_dot(n, dx[2]));
  h = fmin(h, newh);

  dbl3_sub(x[0], x[2], dx[0]);
  dbl3_sub(x[1], x[2], dx[1]);
  dbl3_sub(x[3], x[2], dx[2]);
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
  newh = fabs(dbl3_dot(n, dx[2]));
  h = fmin(h, newh);

  dbl3_sub(x[0], x[3], dx[0]);
  dbl3_sub(x[1], x[3], dx[1]);
  dbl3_sub(x[2], x[3], dx[2]);
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
  newh = fabs(dbl3_dot(n, dx[2]));
  h = fmin(h, newh);

  return h;
}
