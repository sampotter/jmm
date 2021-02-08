#include "geom.h"

#include <assert.h>
#include <math.h>

#include "mat.h"
#include "mesh3.h"
#include "vec.h"

rect3 rect3_make_empty() {
  rect3 rect;
  dbl3_inf(rect.min);
  dbl3_neginf(rect.max);
  return rect;
}

void rect3_get_extent(rect3 const *rect, dbl extent[3]) {
  dbl3_sub(rect->max, rect->min, extent);
}

void rect3_insert_point(rect3 *rect, dbl x[3]) {
  dbl3_min(rect->min, x, rect->min);
  dbl3_max(rect->max, x, rect->max);
}

dbl rect3_surface_area(rect3 const *rect) {
  dbl extent[3];
  rect3_get_extent(rect, extent);
  return 2*dbl3_normsq(extent);
}

bool rect3_overlaps(rect3 const *r1, rect3 const *r2) {
  if (r1->max[0] < r2->min[0] || r1->min[0] > r2->max[0]) return false;
  if (r1->max[1] < r2->min[1] || r1->min[1] > r2->max[1]) return false;
  if (r1->max[2] < r2->min[2] || r1->min[2] > r2->max[2]) return false;
  return true;
}

bool ray3_intersects_rect3(ray3 const *ray, rect3 const *rect) {
  dbl const atol = 1e-15;

  dbl const *p = ray->org, *d = ray->dir;
  dbl const *m = rect->min, *M = rect->max;

  // TODO: handle this case
  assert(fabs(d[0]) > atol || fabs(d[1]) > atol || fabs(d[2]) > atol);

  dbl t, pt[3];

  t = (m[0] - p[0])/d[0];
  dbl3_saxpy(t, d, p, pt);
  if (m[1] <= pt[1] && pt[1] <= M[1] && m[2] <= pt[2] && pt[2] <= M[2])
    return true;

  t = (M[0] - p[0])/d[0];
  dbl3_saxpy(t, d, p, pt);
  if (m[1] <= pt[1] && pt[1] <= M[1] && m[2] <= pt[2] && pt[2] <= M[2])
    return true;

  t = (m[1] - p[1])/d[1];
  dbl3_saxpy(t, d, p, pt);
  if (m[0] <= pt[0] && pt[0] <= M[0] && m[2] <= pt[2] && pt[2] <= M[2])
    return true;

  t = (M[1] - p[1])/d[1];
  dbl3_saxpy(t, d, p, pt);
  if (m[0] <= pt[0] && pt[0] <= M[0] && m[2] <= pt[2] && pt[2] <= M[2])
    return true;

  t = (m[2] - p[2])/d[2];
  dbl3_saxpy(t, d, p, pt);
  if (m[0] <= pt[0] && pt[0] <= M[0] && m[1] <= pt[1] && pt[1] <= M[1])
    return true;

  t = (M[2] - p[2])/d[2];
  dbl3_saxpy(t, d, p, pt);
  if (m[0] <= pt[0] && pt[0] <= M[0] && m[1] <= pt[1] && pt[1] <= M[1])
    return true;

  return false;
}

bool ray3_intersects_tri3(ray3 const *ray, tri3 const *tri, dbl *t) {
  dbl const atol = 1e-15;

  dbl dv[2][3], n[3];
  dbl3_sub(tri->v[1], tri->v[0], dv[0]);
  dbl3_sub(tri->v[2], tri->v[0], dv[1]);
  dbl3_cross(dv[0], dv[1], n);
  dbl3_normalize(n); // TODO: probably unnecessary

  if (fabs(dbl3_dot(n, ray->dir)) < atol)
    assert(false); // TODO: handle this case

  *t = dbl3_dot(n, tri->v[0]) - dbl3_dot(n, ray->org);
  *t /= dbl3_dot(n, ray->dir);

  return *t >= 0;
}

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
