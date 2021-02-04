#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct mesh3 mesh3_s;

typedef struct {
  dbl min[3], max[3];
} rect3;

rect3 rect3_make_empty();
void rect3_get_extent(rect3 const *rect, dbl extent[3]);
void rect3_insert_point(rect3 *rect, dbl x[3]);
dbl rect3_surface_area(rect3 const *rect);
bool rect3_overlaps(rect3 const *r1, rect3 const *r2);

typedef struct {
  dbl v[3][3];
} tri3;

typedef struct {
  /**
   * The origin of the ray.
   */
  dbl org[3];

  /**
   * The direction of the ray.
   */
  dbl dir[3];
} ray3;

bool ray3_intersects_rect3(ray3 const *ray, rect3 const *rect);
bool ray3_intersects_tri3(ray3 const *ray, tri3 const *tri, dbl *t);

/**
 * Check if four points are coplanar.
 */
bool points_are_coplanar(dbl const **x);

bool ray_and_face_are_coplanar(mesh3_s const *mesh, size_t l0, size_t l1,
                               size_t l2, dbl const *ray);

/**
 * Compute the area of a triangle in R^3 with vertices given by `x`,
 * `y`, and `z`.
 */
dbl tri_area(dbl const x[3], dbl const y[3], dbl const z[3]);

/**
 * Assuming that x[0], x[1], and x[2] point to the vertices of a
 * triangle in 3D, compute the barycentric coordinates of y with
 * respect to x[0], x[1], and x[2], storing them in b. This assumes
 * that y is roughly coplanar with x[0], x[1], and x[2], but will
 * forge ahead and do the computation regardless of the exact
 * geometry.
 */
void get_bary_coords_3d(dbl const *x[3], dbl const y[3], dbl b[3]);

/**
 * Check if a pair of tetrahedra in `mesh` is convex. The tetrahedron
 * are required to share the face indexed by `lf`. The indices `l0`
 * and `l1` index the points on either side of the separating face
 * `lf`.
 */
bool adj_tetra_pair_is_convex(mesh3_s const *mesh, size_t l0,
                              size_t const lf[3], size_t l1);

/**
 * Compute the minimum altitude of the tetrahedron with vertices given
 * by `x`.
 */
dbl min_tetra_altitude(dbl const x[4][3]);

// Defined in triBoxOverlap.c. Code by Tomas Akenine-Moller. Checks
// for an overlap between a box and a triangle.
int triBoxOverlap(dbl boxcenter[3], dbl boxhalfsize[3], dbl triverts[3][3]);

#ifdef __cplusplus
}
#endif
