#include "bmesh.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "bb.h"
#include "geom.h"
#include "hybrid.h"
#include "mesh3.h"

typedef struct {
  bb33 const *bb;
  dbl b0[4], db[4];
  dbl level;
} bmesh33_cell_ray_intersects_level_context;

static dbl
bmesh33_cell_ray_intersects_level_hybrid_cost_func(
  dbl t, bmesh33_cell_ray_intersects_level_context *context)
{
  dbl bt[4];
  dbl4_saxpy(t, context->db, context->b0, bt);
  return bb33_f(context->bb, bt) - context->level;
}

bool bmesh33_cell_ray_intersects_level(bmesh33_cell const *cell,
                                       ray3 const *ray, dbl level, dbl b[4]) {
  // TODO: this a first draft... can definitely be
  // simplified/optimized/otherwise improved.

  dbl const atol = 1e-13;

  dbl t0, t1, x[3], b1[4];
  mesh3_tetra_s tetra = {cell->mesh, cell->l};
  bmesh33_cell_ray_intersects_level_context context = {.bb = cell->bb, .level = level};

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

  ray3 ray_ = *ray;
  dbl3_saxpy(t0, ray->dir, ray->org, ray_.org);
  if (ray3_intersects_mesh3_tetra(&ray_, &tetra, &t1)) {
    dbl3_saxpy(t0, ray->dir, ray->org, x);
    mesh3_tetra_get_bary_coords(&tetra, x, context.b0);

    dbl3_saxpy(t1, ray_.dir, ray_.org, x);
    mesh3_tetra_get_bary_coords(&tetra, x, b1);
  } else {
    t1 = t0;
    t0 = 0;

    mesh3_tetra_get_bary_coords(&tetra, ray_.org, context.b0);

    // TODO: a special case to handle later: if the `ray_.org` isn't
    // in `tetra` and we only found one intersection point, then the
    // ray grazed the tetrahedron. We can immediately check whether we
    // hit the level set by just evaluating the Bezier tetrahedron at
    // the intersection point and check whether it's equal to the
    // level value.
    if (!mesh3_tetra_contains_point(&tetra, ray_.org)) {
      printf("t0 = %0.16g\nt1 = %g0.16\n", t0, t1);
      for (int i = 0; i < 4; ++i)
        printf("b0[%d] = %0.16g\n", i, context.b0[i]);
      printf("sum(b0) = %0.16g\n", dbl4_sum(context.b0));
      assert(false);
    }

    dbl3_saxpy(t1, ray_.dir, ray_.org, x);
    mesh3_tetra_get_bary_coords(&tetra, x, b1);
  }
  dbl4_sub(b1, context.b0, context.db);

  // We should handle this case separately---i.e., just check whether
  // the Bezier tetrahedron evaluated at the two coincident interval
  // endpoints equals `level`.
  assert(dbl4_norm(context.db) > atol);

  /**
   * Now we want to use the hybrid rootfinder to check whether f([t0,
   * t1]) contains `level`. It's probably possible to do this here
   * without appealing to a hybrid rootfinder but we can look into
   * that later. An idea is to restrict the Bezier tetra to an
   * interval by pulling out the data for a univariate cubic, but I'm
   * not sure whether a cubic has high enough degree to reproduce the
   * 20-param Bezier tetra exactly or not and too lazy to check right
   * now.
   */
  dbl s = hybrid(
    (hybrid_cost_func_t)bmesh33_cell_ray_intersects_level_hybrid_cost_func,
    0, 1, &context);

  dbl4_saxpy(s, context.db, context.b0, b);

  return true;
}

struct bmesh33 {
  mesh3_s const *mesh;
  bool mesh_owner;
  size_t num_cells;
  bb33 *bb;
};

void bmesh33_alloc(bmesh33_s **bmesh) {
  *bmesh = malloc(sizeof(bmesh33_s));
}

void bmesh33_dealloc(bmesh33_s **bmesh) {
  free(*bmesh);
  *bmesh = NULL;
}

void bmesh33_init_from_mesh3_and_jets(bmesh33_s *bmesh, mesh3_s const *mesh,
                                      jet3 const *jet) {
  bmesh->mesh = mesh;
  bmesh->mesh_owner = false;
  bmesh->num_cells = mesh3_ncells(mesh);

  // Interpolate jets to create Bezier tetrahedra for each cell
  bmesh->bb = malloc(bmesh->num_cells*sizeof(bb33));
  for (size_t l = 0; l < bmesh->num_cells; ++l)
    bb33_init_from_cell_and_jets(&bmesh->bb[l], mesh, jet, l);
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

mesh3_s const *bmesh33_get_mesh_ptr(bmesh33_s const *bmesh) {
  return bmesh->mesh;
}

bmesh33_s *bmesh33_get_level_bmesh(bmesh33_s const *bmesh, dbl level) {
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

  mesh3_alloc((mesh3_s **)&level_bmesh->mesh);
  mesh3_init((mesh3_s *)level_bmesh->mesh, verts, lv, cells, lc, false);

  level_bmesh->mesh_owner = true;
  level_bmesh->num_cells = mesh3_ncells(level_bmesh->mesh);
  assert(level_bmesh->num_cells == lc);

  free(cells);
  free(verts);
  free(brack);

  return level_bmesh;
}

bmesh33_cell bmesh33_get_cell(bmesh33_s const *bmesh, size_t l) {
  assert(l < bmesh->num_cells);
  return (bmesh33_cell) {.bb = &bmesh->bb[l], .mesh = bmesh->mesh, .l = l};
}
