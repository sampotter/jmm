#include <jmm/geom.h>

#include <assert.h>
#include <math.h>
#include <string.h>

#include <jmm/mat.h>
#include <jmm/mesh2.h>
#include <jmm/mesh3.h>
#include <jmm/util.h>
#include <jmm/vec.h>

void line3_get_closest_point(line3 const *line, dbl const x[3], dbl y[3]) {
  dbl xy[3]; dbl3_sub(line->y, line->x, xy);
  dbl t = (dbl3_dot(xy, x) - dbl3_dot(xy, line->x))/dbl3_normsq(xy);
  dbl3_saxpy(t, xy, line->x, y);
}

bool line3_point_colinear(line3 const *line, dbl const x[3], dbl atol) {
  dbl y[3];
  line3_get_closest_point(line, x, y);
  return dbl3_dist(x, y) < atol;
}

bool line3_point_in_interval(line3 const *line, dbl const x[3], dbl atol) {
  dbl xy[3];
  dbl3_sub(line->y, line->x, xy);

  dbl t = (dbl3_dot(xy, x) - dbl3_dot(xy, line->x))/dbl3_normsq(xy);
  if (t < -atol || t > 1 + atol)
    return false;

  dbl y[3];
  dbl3_saxpy(t, xy, line->x, y);
  return dbl3_dist(x, y) < atol;
}

bool mesh3_tetra_contains_point(mesh3_tetra_s const *tetra, dbl const x[3], dbl const *eps) {
  tetra3 tetra_ = mesh3_get_tetra(tetra->mesh, tetra->l);
  return tetra3_contains_point(&tetra_, x, eps);
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

bool tri3_contains_point(tri3 const *tri, dbl const x[3]) {
  dbl const atol = 1e-13;
  dbl normal[3];
  tri3_get_normal(tri, normal);
  if (fabs(dbl3_dot(normal, x) - dbl3_dot(normal, tri->v[0])) > atol)
    return false;
  dbl b[3];
  tri3_get_bary_coords(tri, x, b);
  return dbl3_valid_bary_coord(b);
}

void tri3_get_bary_coords(tri3 const *tri, dbl const x[3], dbl b[3]) {
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

/* Compute the closest point to `x` that lies in `tri` and store the
 * result in `y`. */
void tri3_get_closest_point(tri3 const *tri, dbl const x[3], dbl y[3]) {
  /* This implementation is a transliteration of the "final version of
   * the triangle closest point query" starting on page 141 of
   * "Real-time Collision Detection" (2005). Not totally sure how
   * robust or good this code is, but *do not* want to try
   * implementing this myself at the moment. */

  // Rename variables to match RTCD implementation
  dbl const *p = x, *a = tri->v[0], *b = tri->v[1], *c = tri->v[2];

  // Check if `p` is in the vertex region outside `a`
  dbl ab[3]; dbl3_sub(b, a, ab);
  dbl ac[3]; dbl3_sub(c, a, ac);
  dbl ap[3]; dbl3_sub(p, a, ap);
  dbl d1 = dbl3_dot(ab, ap), d2 = dbl3_dot(ac, ap);
  if (d1 <= 0 && d2 <= 0) {
    dbl3_copy(a, y); // bary coords = (1, 0, 0)
    return;
  }

  // Check if `p` is in vertex region outside `b`
  dbl bp[3]; dbl3_sub(p, b, bp);
  dbl d3 = dbl3_dot(ab, bp);
  dbl d4 = dbl3_dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) {
    dbl3_copy(b, y); // bary coords = (0, 1, 0)
    return;
  }

  // Check if `p` in edge region of `ab`, and if so return projection
  // of `p` onto `ab`
  dbl vc = d1*d4 - d3*d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) {
    dbl v = d1/(d1 - d3);
    return dbl3_saxpy(v, ab, a, y); // bary coords = (1 - v, v, 0)
  }

  // Check if `p` is in the vertex region outside `c`
  dbl cp[3]; dbl3_sub(p, c, cp);
  dbl d5 = dbl3_dot(ab, cp);
  dbl d6 = dbl3_dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) {
    dbl3_copy(c, y); // bary coords = (0, 0, 1)
    return;
  }

  // Check if `p` is in the edge region of `ac`, and if so, return the
  // projection of `p` onto `ac`
  dbl vb = d5*d2 - d1*d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) {
    dbl w = d2/(d2 - d6);
    dbl3_saxpy(w, ac, a, y); // bary coords = (1 - w, 0, w)
    return;
  }

  // Check if `p` is in the edge region of `bc`, and if so, return
  // the projection of `p` onto `bc`
  dbl va = d3*d6 - d5*d4;
  if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) {
    dbl w = (d4 - d3)/(d4 - d3 + d5 - d6);
    dbl bc[3]; dbl3_sub(c, b, bc);
    dbl3_saxpy(w, bc, b, y); // bary coords = (0, 1 - w, w)
    return;
  }

  // Conclude that `p` is inside the face region. Compute `q` through
  // its barycentric coordinates `(u, v, w)`.
  dbl denom = 1/(va + vb + vc);
  dbl v = vb*denom, w = vc*denom;
  dbl tmp[3]; dbl3_saxpy(v, ab, a, tmp);
  dbl3_saxpy(w, ac, tmp, y);
}

dbl tri3_dist(tri3 const *tri, dbl const x[3]) {
  dbl y[3];
  tri3_get_closest_point(tri, x, y);
  return dbl3_dist(x, y);
}

bool tri3_coplanar(tri3 const *tri, tri3 const *other_tri, dbl const *atol) {
  dbl const atol_ = atol ? *atol : 1e-14;
  dbl n[3], other_n[3];
  tri3_get_normal(tri, n);
  tri3_get_normal(other_tri, other_n);
  dbl dot = dbl3_dot(n, other_n);
  return dot > 1 - atol_ || dot < atol_ - 1;
}

bool tri3_equal(tri3 const *tri1, tri3 const *tri2) {
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (tri1->v[i][j] != tri2->v[i][j])
        return false;
  return true;
}

static bool same_side(dbl const a[3], dbl const b[3], dbl const c[3],
                      dbl const p[3], dbl const q[3],
                      dbl const atol) {
  dbl ab[3]; dbl3_sub(b, a, ab);
  dbl ac[3]; dbl3_sub(c, a, ac);
  dbl ap[3]; dbl3_sub(p, a, ap);
  dbl aq[3]; dbl3_sub(q, a, aq);

  dbl n[3]; dbl3_cross(ab, ac, n); dbl3_normalize(n);

  dbl n_dot_ap = shrink(dbl3_ndot(n, ap), atol);
  dbl n_dot_aq = shrink(dbl3_ndot(n, aq), atol);

  return !(signum(n_dot_ap) == -signum(n_dot_aq));
}

rect3 tetra3_get_bounding_box(tetra3 const *tetra) {
  rect3 rect = {
    .min = {INFINITY, INFINITY, INFINITY},
    .max = {-INFINITY, -INFINITY, -INFINITY}
  };

  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      rect.min[j] = fmin(rect.min[j], tetra->v[i][j]);
      rect.max[j] = fmax(rect.max[j], tetra->v[i][j]);
    }
  }

  return rect;
}

bool tetra3_contains_point(tetra3 const *tetra, dbl const x[3], dbl const *eps) {
  dbl const atol = eps ? *eps : 1e-14;
  dbl const (*v)[3] = tetra->v;
  return same_side(v[0], v[1], v[2], v[3], x, atol) &&
         same_side(v[1], v[2], v[3], v[0], x, atol) &&
         same_side(v[2], v[3], v[0], v[1], x, atol) &&
         same_side(v[3], v[0], v[1], v[2], x, atol);
}

void tetra3_get_bary_coords(tetra3 const *tetra, dbl const x[3], dbl b[4]) {
  // Set up the LHS of the problem. The first three rows consist of
  // the components of the tetrahedron vertices, and the last row is
  // all ones.
  dbl lhs[4][4];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 4; ++j)
      lhs[i][j] = tetra->v[j][i];
  for (int j = 0; j < 4; ++j)
    lhs[3][j] = 1;

  // Set up the RHS of the problem. The first three entries are just
  // `x`, and the last is one.
  dbl rhs[4];
  for (int i = 0; i < 3; ++i)
    rhs[i] = x[i];
  rhs[3] = 1;

  // Solve the system. After solving, we should have sum(b) == 1,
  // close to machine precision.
  dbl44_dbl4_solve(lhs, rhs, b);

  dbl4_normalize1(b);
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

grid2_s
tetra3_get_covering_xy_subgrid(tetra3 const *tetra, grid2_s const *grid,
                               size_t offset[2]) {
  rect3 bbox = tetra3_get_bounding_box(tetra);

  dbl2 dx0;
  dbl2_sub(bbox.min, grid->xymin, dx0);

  dbl2 dx1;
  dbl2_sub(bbox.max, grid->xymin, dx1);

  dbl h = grid->h;

  size_t ijmin[2] = {floor(dx0[0]/h), floor(dx0[1]/h)};
  size_t ijmax[2] = {floor(dx1[0]/h), floor(dx1[1]/h)};

  offset[0] = ijmin[0];
  offset[1] = ijmin[1];

  return (grid2_s) {
    .shape = {
      [0] = ijmax[0] - ijmin[0] + 1,
      [1] = ijmax[1] - ijmin[1] + 1
    },
    .xymin = {h*ijmin[0], h*ijmin[1]},
    .h = h,
    .order = grid->order
  };
}

bool tetra3_equal(tetra3 const *tetra1, tetra3 const *tetra2) {
  for (size_t i = 0; i < 4; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (tetra1->v[i][j] != tetra2->v[i][j])
        return false;
  return true;
}

rect3 rect3_make_empty(void) {
  rect3 rect;
  dbl3_inf(rect.min);
  dbl3_neginf(rect.max);
  return rect;
}

rect3 rect3_get_bounding_box_for_points(size_t n, dbl3 const *x) {
  rect3 rect = {
    .min = {INFINITY, INFINITY, INFINITY},
    .max = {-INFINITY, -INFINITY, -INFINITY}
  };

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      rect.min[j] = fmin(rect.min[j], x[i][j]);
      rect.max[j] = fmax(rect.max[j], x[i][j]);
    }
  }

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

bool rect3_is_empty(rect3 const *rect) {
  return rect->min[0] >= rect->max[0] || rect->min[1] >= rect->max[1]
    || rect->min[2] >= rect->max[2];
}

bool rect3_overlaps(rect3 const *r1, rect3 const *r2) {
  if (r1->max[0] < r2->min[0] || r1->min[0] > r2->max[0]) return false;
  if (r1->max[1] < r2->min[1] || r1->min[1] > r2->max[1]) return false;
  if (r1->max[2] < r2->min[2] || r1->min[2] > r2->max[2]) return false;
  return true;
}

bool rect3_occludes_ray3(rect3 const *rect, ray3 const *ray) {
  dbl const *p = ray->org, *d = ray->dir;
  dbl const *m = rect->min, *M = rect->max;

#if JMM_DEBUG
  // TODO: handle this case
  dbl const atol = 1e-14;
  assert(fabs(d[0]) > atol || fabs(d[1]) > atol || fabs(d[2]) > atol);
#endif

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

ray3 ray3_make_empty() {
  return (ray3) {.dir = {NAN, NAN, NAN}, .org = {NAN, NAN, NAN}};
}

void ray3_get_point(ray3 const *ray, dbl t, dbl x[3]) {
  dbl3_saxpy(t, ray->dir, ray->org, x);
}

dbl ray3_intersect_rect3(ray3 const *ray, rect3 const *rect) {
  dbl x = ray->org[0], y = ray->org[1], z = ray->org[2];
  dbl u = ray->dir[0], v = ray->dir[1], w = ray->dir[2];
  dbl x0 = rect->min[0], y0 = rect->min[1], z0 = rect->min[2];
  dbl x1 = rect->max[0], y1 = rect->max[1], z1 = rect->max[2];

  /* ray parameters */
  dbl s, t = INFINITY;
  dbl xs, ys, zs;

  dbl const atol = 1e-15;

  if (fabs(u) > atol) {
    /* Check for intersection with {x0} x [y0, y1] x [z0, z1] */
    s = (x0 - x)/u;
    ys = y + s*v;
    zs = z + s*w;
    if (s >= 0 && y0 <= ys && ys <= y1 && z0 <= zs && zs <= z1)
      t = fmin(s, t);

    /* Check for intersection with {x1} x [y0, y1] x [z0, z1] */
    s = (x1 - x)/u;
    ys = y + s*v;
    zs = z + s*w;
    if (s >= 0 && y0 <= ys && ys <= y1 && z0 <= zs && zs <= z1)
      t = fmin(s, t);
  }

  if (fabs(v) > atol) {
    /* Check for intersection with [x0, x1] x {y0} x [z0, z1] */
    s = (y0 - y)/v;
    xs = x + s*u;
    zs = z + s*w;
    if (s >= 0 && x0 <= xs && xs <= x1 && z0 <= zs && zs <= z1)
      t = fmin(s, t);

    /* Check for intersection with [x0, x1] x {y1} x [z0, z1] */
    s = (y1 - y)/v;
    xs = x + s*u;
    zs = z + s*w;
    if (s >= 0 && x0 <= xs && xs <= x1 && z0 <= zs && zs <= z1)
      t = fmin(s, t);
  }

  if (fabs(w) > atol) {
    /* Check for intersection with [x0, x1] x [y0, y1] x {z0} */
    s = (z0 - z)/w;
    xs = x + s*u;
    ys = y + s*v;
    if (s >= 0 && x0 <= xs && xs <= x1 && y0 <= ys && ys <= y1)
      t = fmin(s, t);

    /* Check for intersection with [x0, x1] x [y0, y1] x {z1} */
    s = (z1 - z)/w;
    xs = x + s*u;
    ys = y + s*v;
    if (s >= 0 && x0 <= xs && xs <= x1 && y0 <= ys && ys <= y1)
      t = fmin(s, t);
  }

  return t;
}

bool ray3_intersects_mesh3_tetra(ray3 const *ray, mesh3_tetra_s const *tetra, dbl *t) {
  tetra3 tetra_ = mesh3_get_tetra(tetra->mesh, tetra->l);
  return ray3_intersects_tetra3(ray, &tetra_, t);
}

/* Try to find `t` such that `ray->org + t*ray->dir` intersects
 * `tri`. This will return `true` if an intersection with `t >= 0`
 * exists, and `false` otherwise. Each argument is assumed to be
 * non-`NULL`. */
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

dbl ray3_closest_point_on_line(ray3 const *ray, line3 const *line,
                               dbl *t_ray, dbl *t_line) {
  dbl dy[3]; dbl3_sub(line->y, line->x, dy);

  dbl tmp1[2][2];
  tmp1[0][0] = dbl3_normsq(ray->dir);
  tmp1[1][0] = tmp1[0][1] = dbl3_dot(ray->dir, dy);
  tmp1[1][1] = dbl3_normsq(dy);

  dbl tmp2[3]; dbl3_sub(ray->org, line->x, tmp2);

  dbl tmp3[2];
  tmp3[0] = dbl3_dot(ray->dir, tmp2);
  tmp3[1] = dbl3_dot(dy, tmp2);

  dbl tmp4[2];
  dbl22_dbl2_solve(tmp1, tmp3, tmp4);
  tmp4[0] *= -1;

  if (t_ray != NULL)
    *t_ray = tmp4[0];

  if (t_line != NULL)
    *t_line = tmp4[1];

  dbl ray_cp[3]; dbl3_saxpy(tmp4[0], ray->dir, ray->org, ray_cp);
  dbl line_cp[3]; dbl3_saxpy(tmp4[1], dy, line->x, line_cp);

  return dbl3_dist(ray_cp, line_cp);
}

bool rect3_contains_point(rect3 const *rect, dbl3 const x) {
  return rect->min[0] <= x[0] && x[0] <= rect->max[0]
    && rect->min[1] <= x[1] && x[1] <= rect->max[1]
    && rect->min[2] <= x[2] && x[2] <= rect->max[2];
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

bool ray3_and_tri3_are_parallel(ray3 const *ray, tri3 const *tri) {
  dbl const atol = 1e-14;
  dbl n[3];
  tri3_get_normal(tri, n);
  return fabs(dbl3_dot(ray->dir, n)) < atol;
}

bool points_are_coplanar(dbl const **x) {
  dbl dx[3][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);
  dbl3_sub(x[3], x[2], dx[2]);
  dbl det = dbl33_det(dx);
  return fabs(det) < 1e-15;
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
