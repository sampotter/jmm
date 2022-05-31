#include "bmesh.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "bb.h"
#include "geom.h"
#include "hybrid.h"
#include "mat.h"
#include "mesh3.h"
#include "util.h"

bool bmesh33_cell_intersect(bmesh33_cell_s const *cell, ray3 const *ray, dbl *t) {
  dbl const atol = 1e-13;

  // TODO: this a first draft... can definitely be
  // simplified/optimized/otherwise improved.

  mesh3_tetra_s tetra = {cell->mesh, cell->l};

  dbl t0, t1, x0[3], x1[3], b0[4], b1[4];

  /**
   * The first thing we need to do is define the interval which we're
   * going to try to intersect with the level set.
   *
   * We try to shoot two rays, since a ray can intersect a tetrahedron
   * in at most two points. If the first ray escapes, we missed. If
   * the second ray (shot from just beyond the first intersection)
   * misses, then we know that [0, t0] is our interval (in this case,
   * we verify that r(0) is inside the tetrahedron---if it isn't, we
   * get another case that we have to handle separately). Otherwise,
   * [t0, t1] is our interval.
   */
  if (!ray3_intersects_mesh3_tetra(ray, &tetra, &t0))
    return false;

  ray3_get_point(ray, t0, x0);

  bool overlap[4];
  int face_inds[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
  tri3 tetra_face[4];
  for (int i = 0; i < 4; ++i) {
    tetra_face[i] = mesh3_tetra_get_face(&tetra, face_inds[i]);
    overlap[i] = tri3_contains_point(&tetra_face[i], x0);
  }

  t1 = INFINITY;
  for (int i = 0; i < 4; ++i) {
    if (overlap[i]) continue;
    dbl tmp;
    if (ray3_intersects_tri3(ray, &tetra_face[i], &tmp))
      t1 = fmin(tmp, t1);
  }

  // If `t1` is infinity, we either grazed the tetrahedron or shot a
  // ray from inside the tetrahedron (or from its boundary). First we
  // set the interval of interest to [0, t0] and then try to handle
  // special cases early.
  if (isinf(t1)) {
    dbl3_copy(x0, x1);
    dbl3_copy(ray->org, x0);

    t1 = t0;
    t0 = 0;

    // If `x0` isn't in `tetra` and we only found one intersection
    // point, then the ray grazed the tetrahedron. We can immediately
    // check whether we hit the level set by evaluating the Bezier
    // tetra to check if the ray hit (which is now parametrized by t =
    // t1) lies on the level set.
    if (!mesh3_tetra_contains_point(&tetra, x0, &atol)) {
      mesh3_tetra_get_bary_coords(&tetra, x1, b1);
      bool hit = fabs(bb33_f(cell->bb, b1) - cell->level) < atol;
      if (hit)
        *t = t1;
      return hit;
    }
  } else {
    ray3_get_point(ray, t1, x1);
  }

  mesh3_tetra_get_bary_coords(&tetra, x0, b0);
  mesh3_tetra_get_bary_coords(&tetra, x1, b1);

  // We should handle this case separately---i.e., just check whether
  // the Bezier tetrahedron evaluated at the two coincident interval
  // endpoints equals `level`.
  if (dbl4_dist(b0, b1) < atol &&
      fabs(bb33_f(cell->bb, b0) - cell->level) < atol) {
    *t = t0;
    return true;
  }

  /**
   * Now we want to use the hybrid rootfinder to check whether f([t0,
   * t1]) contains `level`.
   */
  // TODO: It's probably possible to do this here without appealing to
  // a hybrid rootfinder but we can look into that later. An idea is
  // to restrict the Bezier tetra to an interval by pulling out the
  // data for a univariate cubic, but I'm not sure whether a cubic has
  // high enough degree to reproduce the 20-param Bezier tetra exactly
  // or not and too lazy to check right now.

  cubic_s cubic = bb33_restrict_along_interval(cell->bb, b0, b1);
  cubic_add_constant(&cubic, -cell->level);

  dbl root[3] = {INFINITY, INFINITY, INFINITY};
  int num_roots = cubic_get_real_roots(&cubic, root);

  int bad_roots = 0;
  for (int i = 0; i < num_roots; ++i) {
    if (root[i] < 0 || 1 < root[i]) {
      root[i] = INFINITY;
      ++bad_roots;
    }
  }
  num_roots -= bad_roots;

  dbl3_sort(root);

  /* Compute `t` as the arc-length distance from the ray origin to the
   * first intersection with the tetrahedron plus the arc length
   * distance from `x0` to `xs`, which is just s scaled by
   * `dbl3_dist(x0, x1)`.
   *
   * If we didn't find any roots in the interval [0, 1], root[0] will
   * be infinity, which will keep *t set to infinity after
   * returning. */
  *t = t0 + root[0]*dbl3_dist(x0, x1);

  return isfinite(*t);
}

void bmesh33_cell_Df(bmesh33_cell_s const *cell, dbl3 const x, dbl3 Df) {
  mesh3_tetra_s tetra = {cell->mesh, cell->l};

  /* Compute the gradient at x in barycentric coordinates */
  dbl4 b;
  mesh3_tetra_get_bary_coords(&tetra, x, b);

  /* Set up transform matrix. The first three rows of A correspond to
   * the standard directions in R^3, which we use to compute the
   * gradient. */
  static dbl44 A;
  size_t lv[4];
  mesh3_cv(cell->mesh, cell->l, lv);
  for (size_t i = 0; i < 4; ++i) {
    dbl const *xi = mesh3_get_vert_ptr(cell->mesh, lv[i]);
    for (size_t j = 0; j < 3; ++j)
      A[j][i] = xi[j];
    A[3][i] = 1;
  }
  dbl44_invert(A);
  dbl44_transpose(A);

  dbl const atol = 1e-13;
  assert(fabs(dbl4_sum(A[0])) < atol);
  assert(fabs(dbl4_sum(A[1])) < atol);
  assert(fabs(dbl4_sum(A[2])) < atol);
  assert(fabs(1 - dbl4_sum(A[3])) < atol);

  /* Convert back to Cartesian coordinates */
  for (size_t i = 0; i < 3; ++i)
    Df[i] = bb33_df(cell->bb, b, A[i]);
}

bool bmesh33_cell_equal(bmesh33_cell_s const *c1, bmesh33_cell_s const *c2) {
  return c1->mesh == c2-> mesh && c1->l == c2->l;
}

struct bmesh33 {
  mesh3_s const *mesh;
  bool mesh_owner;
  size_t num_cells;
  bb33 *bb;
  dbl level;
};

void bmesh33_alloc(bmesh33_s **bmesh) {
  *bmesh = malloc(sizeof(bmesh33_s));
}

void bmesh33_dealloc(bmesh33_s **bmesh) {
  free(*bmesh);
  *bmesh = NULL;
}

void bmesh33_init_from_mesh3_and_jets(bmesh33_s *bmesh, mesh3_s const *mesh,
                                      jet31t const *jet) {
  bmesh->mesh = mesh;
  bmesh->mesh_owner = false;
  bmesh->num_cells = mesh3_ncells(mesh);

  // Interpolate jets to create Bezier tetrahedra for each cell
  bmesh->bb = malloc(bmesh->num_cells*sizeof(bb33));
  for (size_t l = 0; l < bmesh->num_cells; ++l)
    bb33_init_from_cell_and_jets(&bmesh->bb[l], mesh, jet, l);

  bmesh->level = NAN;
}

void bmesh33_deinit(bmesh33_s *bmesh) {
  free(bmesh->bb);
  bmesh->bb = NULL;

  if (bmesh->mesh_owner) {
    mesh3_deinit((mesh3_s *)bmesh->mesh);
    mesh3_dealloc((mesh3_s **)&bmesh->mesh);
  }
}

size_t bmesh33_num_cells(bmesh33_s const *bmesh) {
  return bmesh->num_cells;
}

dbl bmesh33_get_level(bmesh33_s const *bmesh) {
  return bmesh->level;
}

mesh3_s const *bmesh33_get_mesh_ptr(bmesh33_s const *bmesh) {
  return bmesh->mesh;
}

bmesh33_s *bmesh33_restrict_to_level(bmesh33_s const *bmesh, dbl level) {
  // First, precompute which Bezier tetra *might* bracket the
  // level. Since the graph of each Bezier tetra is contained in the
  // convex hull of the defining control points, we can quickly rule
  // out points outside of this interval (just check if the level is
  // greater than the maximum or smaller than the minimum Bezier
  // ordinate), but we may still have some false positives in what we
  // compute this way.
  size_t num_brack = 0;
  bool *brack = malloc(bmesh->num_cells*sizeof(bool));
  for (size_t l = 0; l < bmesh->num_cells; ++l) {
    brack[l] = bb33_convex_hull_brackets_value(&bmesh->bb[l], level);
    num_brack += brack[l];
  }

  bmesh33_s *level_bmesh = malloc(sizeof(bmesh33_s));
  level_bmesh->bb = malloc(num_brack*sizeof(bb33));

  dbl *verts = malloc(4*num_brack*sizeof(dbl[3]));
  size_t *cells = malloc(num_brack*sizeof(size_t[4]));

  // For each bracketed cell, we go ahead and splat all of the
  // incident vertices into `verts`. There will be duplicates. For now
  // we're trading off memory use for simplicity. Later on we can come
  // back and optimize this. While we're doing this, we also copy over
  // Bezier tetra data so that we don't have to recompute it.
  size_t lv = 0, lc = 0, cv[4], j = 0;
  for (size_t l = 0; l < bmesh->num_cells; ++l) {
    if (!brack[l]) continue;
    level_bmesh->bb[lc] = bmesh->bb[l]; // Copy Bezier tetra data
    mesh3_cv(bmesh->mesh, l, cv); // Grab the inds for this cell
    for (int i = 0; i < 4; ++i) {
      mesh3_copy_vert(bmesh->mesh, cv[i], &verts[4*3*lc + 3*i]);
      cells[4*lc + i] = j++; // Splat new vert inds into `cells`
    }
    lv += 4;
    ++lc;
  }

  dbl eps = mesh3_get_eps(bmesh->mesh);

  mesh3_alloc((mesh3_s **)&level_bmesh->mesh);
  mesh3_init((mesh3_s *)level_bmesh->mesh, verts, lv, cells, lc, false, &eps);

  level_bmesh->mesh_owner = true;
  level_bmesh->num_cells = mesh3_ncells(level_bmesh->mesh);
  assert(level_bmesh->num_cells == lc);

  level_bmesh->level = level;

  free(cells);
  free(verts);
  free(brack);

  return level_bmesh;
}

bmesh33_cell_s bmesh33_get_cell(bmesh33_s const *bmesh, size_t l) {
  assert(l < bmesh->num_cells);
  return (bmesh33_cell_s) {
    .bb = &bmesh->bb[l],
    .mesh = bmesh->mesh,
    .l = l,
    .level = bmesh33_get_level(bmesh)
  };
}

/* Evaluate `bmesh` at the point `x`. If `x` lies outside the mesh,
 * return `NAN`. */
dbl bmesh33_f(bmesh33_s const *bmesh, dbl3 const x) {
  // TODO: very inefficient implementation! Optimize this using rtree.
  for (size_t l = 0; l < bmesh->num_cells; ++l) {
    if (mesh3_cell_contains_point(bmesh->mesh, l, x)) {
      tetra3 tetra = mesh3_get_tetra(bmesh->mesh, l);
      dbl4 b;
      tetra3_get_bary_coords(&tetra, x, b);
      return bb33_f(&bmesh->bb[l], b);
    }
  }
  return NAN;
}

bb33 *bmesh33_get_bb_ptr(bmesh33_s const *bmesh, size_t lc) {
  return &bmesh->bb[lc];
}
