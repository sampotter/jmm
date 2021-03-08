#include "geom.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "mat.h"
#include "mesh2.h"
#include "mesh3.h"
#include "vec.h"

bool mesh3_tetra_contains_point(mesh3_tetra_s const *tetra, dbl const x[3]) {
  tetra3 tetra_ = mesh3_get_tetra(tetra->mesh, tetra->l);
  return tetra3_contains_point(&tetra_, x);
}

void mesh3_tetra_get_bary_coords(mesh3_tetra_s const *tetra, dbl const x[3], dbl b[4]) {
  tetra3 tetra_ = mesh3_get_tetra(tetra->mesh, tetra->l);
  tetra3_get_bary_coords(&tetra_, x, b);
}

void mesh3_tetra_get_point(mesh3_tetra_s const *tetra, dbl const b[4], dbl x[3]) {
  tetra3 tetra_ = mesh3_get_tetra(tetra->mesh, tetra->l);
  tetra3_get_point(&tetra_, b, x);
}

void tri3_get_normal(tri3 const *tri, dbl normal[3]) {
  dbl t[2][3];
  dbl3_sub(tri->v[1], tri->v[0], t[0]);
  dbl3_sub(tri->v[2], tri->v[0], t[1]);
  dbl3_cross(t[0], t[1], normal);
  dbl3_normalize(normal);
}

void tri3_get_centroid(tri3 const *tri, dbl c[3]) {
  for (int j = 0; j < 3; ++j) {
    c[j] = 0;
    for (int i = 0; i < 3; ++i)
      c[j] += tri->v[i][j];
    c[j] /= 3;
  }
}

bool tri3_contains_point(tri3 const *tri, dbl x[3]) {
  dbl const atol = 1e-13;
  dbl normal[3];
  tri3_get_normal(tri, normal);
  if (fabs(dbl3_dot(normal, x) - dbl3_dot(normal, tri->v[0])) > atol)
    return false;
  dbl b[3];
  tri3_get_bary_coords(tri, x, b);
  return dbl3_valid_bary_coord(b);
}

void tri3_get_bary_coords(tri3 const *tri, dbl const x[3], dbl b[4]) {
  dbl dx[2][3];
  dbl3_sub(tri->v[1], tri->v[0], dx[0]);
  dbl3_sub(tri->v[2], tri->v[0], dx[1]);

  dbl n[3];
  dbl3_cross(dx[0], dx[1], n);
  dbl area = dbl3_normalize(n)/2;

  b[0] = tri_area(x, tri->v[1], tri->v[2])/area;
  b[1] = tri_area(tri->v[0], x, tri->v[2])/area;
  b[2] = tri_area(tri->v[0], tri->v[1], x)/area;
}

bool tetra3_contains_point(tetra3 const *tetra, dbl const x[3]) {
  dbl b[4];
  tetra3_get_bary_coords(tetra, x, b);
  return dbl4_valid_bary_coord(b);
}

void tetra3_get_bary_coords(tetra3 const *tetra, dbl const x[3], dbl b[4]) {
  dbl lhs[4][4];
  for (int j = 0; j < 4; ++j) {
    lhs[0][j] = 1;
    for (int i = 1; i < 4; ++i)
      lhs[i][j] = tetra->v[j][i - 1];
  }

  dbl rhs[4] = {1, x[0], x[1], x[2]}, vol = dbl44_det(lhs);

  dbl tmp[4];
  for (int j = 0; j < 4; ++j) {
    dbl44_get_col(lhs, j, tmp);
    dbl44_set_col(lhs, j, rhs);
    b[j] = dbl44_det(lhs)/vol;
    dbl44_set_col(lhs, j, tmp);
  }
}

void tetra3_get_centroid(tetra3 const *tetra, dbl c[3]) {
  for (int j = 0; j < 3; ++j) {
    c[j] = 0;
    for (int i = 0; i < 4; ++i)
      c[j] += tetra->v[i][j];
    c[j] /= 4;
  }
}

void tetra3_get_point(tetra3 const *tetra, dbl const b[4], dbl x[3]) {
  x[0] = x[1] = x[2] = 0;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j)
      x[j] += b[i]*tetra->v[i][j];
}

rect3 rect3_make_empty(void) {
  rect3 rect;
  dbl3_inf(rect.min);
  dbl3_neginf(rect.max);
  return rect;
}

void rect3_get_extent(rect3 const *rect, dbl extent[3]) {
  dbl3_sub(rect->max, rect->min, extent);
}

void rect3_get_centroid(rect3 const *rect, dbl centroid[3]) {
  centroid[0] = (rect->min[0] + rect->max[0])/2;
  centroid[1] = (rect->min[1] + rect->max[1])/2;
  centroid[2] = (rect->min[2] + rect->max[2])/2;
}

void rect3_get_half_extent(rect3 const *rect, dbl half_extent[3]) {
  half_extent[0] = (rect->max[0] - rect->min[0])/2;
  half_extent[1] = (rect->max[1] - rect->min[1])/2;
  half_extent[2] = (rect->max[2] - rect->min[2])/2;
}

void rect3_insert_point(rect3 *rect, dbl const x[3]) {
  dbl3_min(rect->min, x, rect->min);
  dbl3_max(rect->max, x, rect->max);
}

void rect3_insert_tri3(rect3 *rect, tri3 const *tri) {
  for (int i = 0; i < 3; ++i) {
    dbl3_min(rect->min, tri->v[i], rect->min);
    dbl3_max(rect->max, tri->v[i], rect->max);
  }
}

void rect3_insert_tetra3(rect3 *rect, tetra3 const *tetra) {
  for (int i = 0; i < 4; ++i) {
    dbl3_min(rect->min, tetra->v[i], rect->min);
    dbl3_max(rect->max, tetra->v[i], rect->max);
  }
}

void rect3_insert_mesh2_tri(rect3 *rect, mesh2_tri_s const *tri) {
  tri3 tri_ = mesh2_get_tri(tri->mesh, tri->l);
  rect3_insert_tri3(rect, &tri_);
}

void rect3_insert_mesh3_tetra(rect3 *rect, mesh3_tetra_s const *tetra) {
  tetra3 tetra_ = mesh3_get_tetra(tetra->mesh, tetra->l);
  rect3_insert_tetra3(rect, &tetra_);
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

bool rect3_occludes_ray3(rect3 const *rect, ray3 const *ray) {
  dbl const atol = 1e-15;

  dbl const *p = ray->org, *d = ray->dir;
  dbl const *m = rect->min, *M = rect->max;

  // TODO: handle this case
  assert(fabs(d[0]) > atol || fabs(d[1]) > atol || fabs(d[2]) > atol);

  dbl t, pt[3];

  t = (m[0] - p[0])/d[0];
  dbl3_saxpy(t, d, p, pt);
  if (t >= 0 && m[1] <= pt[1] && pt[1] <= M[1] && m[2] <= pt[2] && pt[2] <= M[2])
    return true;

  t = (M[0] - p[0])/d[0];
  dbl3_saxpy(t, d, p, pt);
  if (t >= 0 && m[1] <= pt[1] && pt[1] <= M[1] && m[2] <= pt[2] && pt[2] <= M[2])
    return true;

  t = (m[1] - p[1])/d[1];
  dbl3_saxpy(t, d, p, pt);
  if (t >= 0 && m[0] <= pt[0] && pt[0] <= M[0] && m[2] <= pt[2] && pt[2] <= M[2])
    return true;

  t = (M[1] - p[1])/d[1];
  dbl3_saxpy(t, d, p, pt);
  if (t >= 0 && m[0] <= pt[0] && pt[0] <= M[0] && m[2] <= pt[2] && pt[2] <= M[2])
    return true;

  t = (m[2] - p[2])/d[2];
  dbl3_saxpy(t, d, p, pt);
  if (t >= 0 && m[0] <= pt[0] && pt[0] <= M[0] && m[1] <= pt[1] && pt[1] <= M[1])
    return true;

  t = (M[2] - p[2])/d[2];
  dbl3_saxpy(t, d, p, pt);
  if (t >= 0 && m[0] <= pt[0] && pt[0] <= M[0] && m[1] <= pt[1] && pt[1] <= M[1])
    return true;

  return false;
}

void ray3_get_point(ray3 const *ray, dbl t, dbl x[3]) {
  dbl3_saxpy(t, ray->dir, ray->org, x);
}

bool ray3_intersects_rect3(ray3 const *ray, rect3 const *rect, dbl *t) {
  dbl const atol = 1e-15;

  dbl const *p = ray->org, *d = ray->dir;
  dbl const *m = rect->min, *M = rect->max;

  // TODO: handle this case
  assert(fabs(d[0]) > atol || fabs(d[1]) > atol || fabs(d[2]) > atol);

  dbl s, pt[3];

  *t = INFINITY;

  s = (m[0] - p[0])/d[0];
  dbl3_saxpy(s, d, p, pt);
  if (s >= 0 && m[1] <= pt[1] && pt[1] <= M[1] && m[2] <= pt[2] && pt[2] <= M[2])
    *t = fmin(*t, s);

  s = (M[0] - p[0])/d[0];
  dbl3_saxpy(s, d, p, pt);
  if (s >= 0 && m[1] <= pt[1] && pt[1] <= M[1] && m[2] <= pt[2] && pt[2] <= M[2])
    *t = fmin(*t, s);

  s = (m[1] - p[1])/d[1];
  dbl3_saxpy(s, d, p, pt);
  if (s >= 0 && m[0] <= pt[0] && pt[0] <= M[0] && m[2] <= pt[2] && pt[2] <= M[2])
    *t = fmin(*t, s);

  s = (M[1] - p[1])/d[1];
  dbl3_saxpy(s, d, p, pt);
  if (s >= 0 && m[0] <= pt[0] && pt[0] <= M[0] && m[2] <= pt[2] && pt[2] <= M[2])
    *t = fmin(*t, s);

  s = (m[2] - p[2])/d[2];
  dbl3_saxpy(s, d, p, pt);
  if (s >= 0 && m[0] <= pt[0] && pt[0] <= M[0] && m[1] <= pt[1] && pt[1] <= M[1])
    *t = fmin(*t, s);

  s = (M[2] - p[2])/d[2];
  dbl3_saxpy(s, d, p, pt);
  if (s >= 0 && m[0] <= pt[0] && pt[0] <= M[0] && m[1] <= pt[1] && pt[1] <= M[1])
    *t = fmin(*t, s);

  return *t >= 0;
}

bool ray3_intersects_mesh3_tetra(ray3 const *ray, mesh3_tetra_s const *tetra, dbl *t) {
  tetra3 tetra_ = mesh3_get_tetra(tetra->mesh, tetra->l);
  return ray3_intersects_tetra3(ray, &tetra_, t);
}

bool ray3_intersects_tri3(ray3 const *ray, tri3 const *tri, dbl *t) {
  dbl const atol = 1e-15;

  dbl const *v[3] = {tri->v[0], tri->v[1], tri->v[2]};
  dbl dv[2][3], n[3], b[3], xt[3];
  dbl3_sub(tri->v[1], tri->v[0], dv[0]);
  dbl3_sub(tri->v[2], tri->v[0], dv[1]);
  dbl3_cross(dv[0], dv[1], n);
  dbl3_normalize(n); // TODO: probably unnecessary

  // First, check if the origin of the ray lies in the plane of the
  // triangle.
  *t = dbl3_dot(n, tri->v[0]) - dbl3_dot(n, ray->org);
  if (fabs(*t) < atol) {
    get_bary_coords_3d(v, ray->org, b);
    return dbl3_valid_bary_coord(b);
  }

  // If the ray direction and triangle normal are orthogonal, then the
  // triangle would appear infinitely thin to the ray, so return false
  if (fabs(dbl3_dot(n, ray->dir)) < atol)
    return false;

  // Compute the ray parameter---if it's negative, the triangle is
  // behind the start of the ray and we can return early.
  *t /= dbl3_dot(n, ray->dir);
  if (*t < 0)
    return false;

  // Finally, check whether the ray passes through the triangle
  ray3_get_point(ray, *t, xt);
  get_bary_coords_3d(v, xt, b);
  return dbl3_valid_bary_coord(b);
}

static void tetra3_get_face(tetra3 const *tetra, int p[3], tri3 *tri) {
  assert(p[0] != p[1] && p[1] != p[2]);
  for (int i = 0; i < 3; ++i) {
    assert(0 <= p[i] && p[i] < 4);
    for (int j = 0; j < 3; ++j)
      tri->v[i][j] = tetra->v[p[i]][j];
  }
}

bool ray3_intersects_tetra3(ray3 const *ray, tetra3 const *tetra, dbl *t) {
  tri3 tri;
  dbl s;

  *t = INFINITY;

  /**
   * Get each of the faces of tetra using `tetra3_get_face`, intersect
   * `ray` with these faces, and set `*t` to the minimum intersection
   * parameter.
   */

  tetra3_get_face(tetra, (int[]) {0, 1, 2}, &tri);
  if (ray3_intersects_tri3(ray, &tri, &s))
    *t = fmin(*t, s);

  tetra3_get_face(tetra, (int[]) {0, 1, 3}, &tri);
  if (ray3_intersects_tri3(ray, &tri, &s))
    *t = fmin(*t, s);

  tetra3_get_face(tetra, (int[]) {0, 2, 3}, &tri);
  if (ray3_intersects_tri3(ray, &tri, &s))
    *t = fmin(*t, s);

  tetra3_get_face(tetra, (int[]) {1, 2, 3}, &tri);
  if (ray3_intersects_tri3(ray, &tri, &s))
    *t = fmin(*t, s);

  // We've found an intersection if `*t < INFINITY` (its initial value
  // at the start of `ray3_intersects_tetra`).
  return isfinite(*t);
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
