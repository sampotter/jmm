#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "def.h"

typedef struct {
  dbl x[3], y[3];
} line3;

void line3_get_closest_point(line3 const *line, dbl const x[3], dbl y[3]);
bool line3_point_colinear(line3 const *line, dbl const x[3], dbl atol);

typedef struct {
  dbl min[3], max[3];
} rect3;

bool mesh3_tetra_contains_point(mesh3_tetra_s const *tetra, dbl const x[3], dbl const *eps);
void mesh3_tetra_get_bary_coords(mesh3_tetra_s const *tetra, dbl const x[3], dbl b[4]);
void mesh3_tetra_get_point(mesh3_tetra_s const *tetra, dbl const b[4], dbl x[3]);

typedef struct {
  dbl v[3][3];
} tri3;

void tri3_get_normal(tri3 const *tri, dbl normal[3]);
void tri3_get_centroid(tri3 const *tri, dbl centroid[3]);
bool tri3_contains_point(tri3 const *tri, dbl const x[3]);
void tri3_get_bary_coords(tri3 const *tri, dbl const x[3], dbl b[3]);
void tri3_get_closest_point(tri3 const *tri, dbl const x[3], dbl y[3]);
dbl tri3_dist(tri3 const *tri, dbl const x[3]);
bool tri3_coplanar(tri3 const *tri, tri3 const *other_tri, dbl const *atol);

typedef struct {
  dbl v[4][3];
} tetra3;

bool tetra3_contains_point(tetra3 const *tetra, dbl const x[3], dbl const *eps);
void tetra3_get_bary_coords(tetra3 const *tetra, dbl const x[3], dbl b[4]);
void tetra3_get_centroid(tetra3 const *tetra, dbl centroid[3]);
void tetra3_get_point(tetra3 const *tetra, dbl const b[4], dbl x[3]);

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

rect3 rect3_make_empty(void);
void rect3_get_extent(rect3 const *rect, dbl extent[3]);
void rect3_get_centroid(rect3 const *rect, dbl centroid[3]);
void rect3_get_half_extent(rect3 const *rect, dbl half_extent[3]);
void rect3_insert_point(rect3 *rect, dbl const x[3]);
void rect3_insert_tri3(rect3 *rect, tri3 const *tri);
void rect3_insert_tetra3(rect3 *rect, tetra3 const *tetra);
void rect3_insert_mesh2_tri(rect3 *rect, mesh2_tri_s const *tri);
void rect3_insert_mesh3_tetra(rect3 *rect, mesh3_tetra_s const *tetra);
dbl rect3_surface_area(rect3 const *rect);
bool rect3_overlaps(rect3 const *r1, rect3 const *r2);
bool rect3_occludes_ray3(rect3 const *rect, ray3 const *ray);

ray3 ray3_make_empty();
void ray3_get_point(ray3 const *ray, dbl t, dbl x[3]);
bool ray3_intersects_rect3(ray3 const *ray, rect3 const *rect, dbl *t);
bool ray3_intersects_mesh3_tetra(ray3 const *ray, mesh3_tetra_s const *tetra, dbl *t);
bool ray3_intersects_tri3(ray3 const *ray, tri3 const *tri, dbl *t);
bool ray3_intersects_tetra3(ray3 const *ray, tetra3 const *tetra, dbl *t);
bool ray3_and_tri3_are_parallel(ray3 const *ray, tri3 const *tri);

/**
 * Check if four points are coplanar.
 */
bool points_are_coplanar(dbl const **x);

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
 * Compute the minimum altitude of the tetrahedron with vertices given
 * by `x`.
 */
dbl min_tetra_altitude(dbl const x[4][3]);

// Defined in triBoxOverlap.c. Code by Tomas Akenine-Moller. Checks
// for an overlap between a box and a triangle.
int triBoxOverlap(dbl const center[3], dbl const half[3], dbl const v[3][3]);

#ifdef __cplusplus
}
#endif
