#include <cgreen/cgreen.h>

#include "bmesh.h"
#include "camera.h"
#include "mesh3.h"
#include "rtree.h"

Describe(bmesh33);

BeforeEach(bmesh33) {
  double_absolute_tolerance_is(2e-14);
  double_relative_tolerance_is(1e-15);
}

AfterEach(bmesh33) {}

/**
 * This function creates a `bmesh33_s` that approximates the sphere
 * over a tetrahedron mesh that discretizes [-1, 1]^3, where each
 * octant is discretized into five tetrahedra so that the whole thing
 * is perfectly symmetric.
 */
static void
create_approximate_sphere_bmesh33(
  bmesh33_s **bmesh_handle, mesh3_s **mesh_handle, jet3 **jet_handle)
{
  /**
   * First, set up the `mesh3_s` discretizing [-1, 1]^3.
   */

  // Set up the first eight vertices as a template (we'll mirror these
  // onto the other eight octants).
  dbl verts[64][3] = {
    [0] = {0, 0, 0},
    [1] = {0, 0, 1},
    [2] = {0, 1, 0},
    [3] = {0, 1, 1},
    [4] = {1, 0, 0},
    [5] = {1, 0, 1},
    [6] = {1, 1, 0},
    [7] = {1, 1, 1}
  };

  // Explicitly lay out the signs for each reflection of the first
  // eight template vertices onto the remaining seven octans.
  int signs[8][3] = {
    { 1,  1,  1},
    { 1,  1, -1},
    { 1, -1,  1},
    { 1, -1, -1},
    {-1,  1,  1},
    {-1,  1, -1},
    {-1, -1,  1},
    {-1, -1, -1}
  };

  // Iterate through the first eight octants repeatedly, mirroring
  // them onto the remaining octants. There will be duplicate
  // vertices, which is perfectly fine.
  for (int i = 1; i < 8; ++i)
    for (int j = 0; j < 8; ++j)
      for (int k = 0; k < 3; ++k)
        verts[8*i + j][k] = signs[i][k]*verts[j][k];

  // Lay out the template cell, which should be symmetric under
  // permutation of the orthant axes.
  size_t cells[40][4] = {
    [0] = {0, 1, 3, 5}, // 000, 001, 011, 101
    [1] = {0, 2, 3, 6}, // 000, 010, 011, 110
    [2] = {0, 4, 5, 6}, // 000, 100, 101, 110
    [3] = {0, 3, 5, 6}, // 000, 011, 101, 110
    [4] = {3, 5, 6, 7}, // 011, 101, 110, 111
  };

  // Now, repeat the template cell's index list by shifting it for
  // each other orthant. (Nothing fancy needs to happen here.)
  for (int i = 1; i < 8; ++i)
    for (int j = 0; j < 5; ++j)
      for (int k = 0; k < 4; ++k)
        cells[5*i + j][k] = 8*i + cells[j][k];

  mesh3_alloc(mesh_handle);
  mesh3_init(*mesh_handle, &verts[0][0], 64, &cells[0][0], 40, false, NULL);

  /**
   * Next, compute the jets for each vertex in `verts`.
   */

  *jet_handle = malloc(64*sizeof(jet3));

  for (int i = 0; i < 8; ++i) {
    // We specially set the jet at (0, 0, 0) to zero to avoid a
    // singularity. This makes the overall approximation worse, but
    // this is just a test...
    (*jet_handle)[8*i] = (jet3) {.f = 0, .Df = {0, 0, 0}};

    dbl *x;
    for (int j = 1; j < 8; ++j) {
      x = verts[8*i + j];
      jet3 J = {.f = dbl3_norm(x)};
      dbl3_normalized(x, J.Df);
      (*jet_handle)[8*i + j] = J;
    }
  }

  /**
   * Finally, set up the `bmesh33_s` instance.
   */

  bmesh33_alloc(bmesh_handle);
  bmesh33_init_from_mesh3_and_jets(*bmesh_handle, *mesh_handle, *jet_handle);
}

static void
destroy_approximate_sphere_bmesh33(
  bmesh33_s **bmesh_handle, mesh3_s **mesh_handle, jet3 **jet_handle)
{
  bmesh33_deinit(*bmesh_handle);
  bmesh33_dealloc(bmesh_handle);

  mesh3_deinit(*mesh_handle);
  mesh3_dealloc(mesh_handle);

  free(*jet_handle);
  *jet_handle = NULL;
}

#define SET_UP_APPROXIMATE_SPHERE()                         \
  bmesh33_s *bmesh;                                         \
  mesh3_s *mesh;                                            \
  jet3 *jet;                                                \
  create_approximate_sphere_bmesh33(&bmesh, &mesh, &jet);

#define TEAR_DOWN_APPROXIMATE_SPHERE()                      \
  destroy_approximate_sphere_bmesh33(&bmesh, &mesh, &jet);

Ensure(bmesh33, approximate_sphere_setup_and_teardown_works) {
  SET_UP_APPROXIMATE_SPHERE();

  /**
   * Just check basic properties of the structs initialized by
   * SET_UP_APPROXIMATE_SPHERE.
   */
  assert_that(bmesh33_num_cells(bmesh), is_equal_to(40));

  TEAR_DOWN_APPROXIMATE_SPHERE();
}

Ensure(bmesh33, mesh3_cell_contains_point_works) {
  SET_UP_APPROXIMATE_SPHERE();

  size_t num_cells = mesh3_ncells(mesh);

  dbl3 x;

  for (size_t i = 0; i < num_cells; ++i) {
    mesh3_get_centroid(mesh, i, x);
    for (size_t j = 0; j < num_cells; ++j) {
      if (i == j) {
        assert_that(mesh3_cell_contains_point(mesh, j, x));
      } else {
        assert_false(mesh3_cell_contains_point(mesh, j, x));
      }
    }
  }

  TEAR_DOWN_APPROXIMATE_SPHERE();
}

Ensure(bmesh33, ray_intersects_level_works_on_approximate_sphere) {
  SET_UP_APPROXIMATE_SPHERE();

  dbl const level = 0.5;

  bmesh33_s *level_bmesh = bmesh33_restrict_to_level(bmesh, level);

  /* After restricting to T = 0.5, we lose the corner tetrahedra,
   * leaving 32 cells in the level's tetrahedron mesh (i.e., the ones
   * containing (+/- 1, +/- 1, +/- 1) in each octant). */
  assert_that(bmesh33_num_cells(level_bmesh), is_equal_to(32));

  rtree_s *rtree;
  rtree_alloc(&rtree);
  rtree_init(rtree, 4, RTREE_SPLIT_STRATEGY_SURFACE_AREA);
  rtree_insert_bmesh33(rtree, level_bmesh);

  camera_s camera = {
    .type = CAMERA_TYPE_ORTHOGRAPHIC,
    .pos = {0, -2, 0},
    .look = {0, 1, 0},
    .left = {-1, 0, 0},
    .up = {0, 0, 1},
    .width = 2.0,
    .height = 2.0,
    .dim = {33, 33}
  };

  FILE *fp = fopen(
    "bmesh33_ray_intersects_level_works_on_approximate_sphere.txt", "r");

  ray3 ray;
  isect isect;
  dbl t_gt;
  for (size_t i = 0; i < camera.dim[0]; ++i){
    for (size_t j = 0; j < camera.dim[1]; ++j) {
      // Shoot the (i, j)th camera ray
      ray = camera_get_ray_for_index(&camera, i, j);
      rtree_intersect(rtree, &ray, &isect);

      // Read the correct groundtruth value for the intersection
      // parameter from disk
      fscanf(fp, "%lf\n", &t_gt);

      // Check that isect.t agrees with the groundtruth value
      if (isinf(t_gt)) {
        assert_that(isinf(isect.t));
        assert_that(isect.obj, is_null);
      } else {
        assert_that_double(t_gt, is_nearly_double(isect.t));
        assert_that(robj_get_type(isect.obj), is_equal_to(ROBJ_BMESH33_CELL));
      }

      if (isect.obj == NULL)
        continue;

      // Make sure that the intersected point lies on the correct
      // level set!
      dbl3 xt; ray3_get_point(&ray, isect.t, xt);
      bmesh33_cell_s cell = *(bmesh33_cell_s *)robj_get_data(isect.obj);
      tetra3 tetra = mesh3_get_tetra(cell.mesh, cell.l);
      dbl4 b; tetra3_get_bary_coords(&tetra, xt, b);
      dbl f = bb33_f(cell.bb, b);
      assert_that_double(f, is_nearly_double(level));
    }
  }

  fclose(fp);

  rtree_deinit(rtree);
  rtree_dealloc(&rtree);

  bmesh33_deinit(level_bmesh);
  bmesh33_dealloc(&level_bmesh);

  TEAR_DOWN_APPROXIMATE_SPHERE();
}
