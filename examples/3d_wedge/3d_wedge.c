#include "3d_wedge.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <error.h>
#include <hybrid.h>
#include <mat.h>
#include <mesh2.h>

#include "mesh3_extra.h"

void jmm_3d_wedge_spec_dump(jmm_3d_wedge_spec_s const *spec, char const *path) {
  FILE *fp = fopen(path, "w");

  fprintf(fp, "verbose: %s\n", spec->verbose ? "true" : "false");
  fprintf(fp, "visualize: %s\n", spec->visualize ? "true" : "false");
  fprintf(fp, "maxvol: %g\n", spec->maxvol);
  fprintf(fp, "n: %g\n", spec->n);
  fprintf(fp, "w: %g\n", spec->w);
  fprintf(fp, "h: %g\n", spec->h);
  fprintf(fp, "R: %g\n", spec->R);

  fclose(fp);
}

void jmm_3d_wedge_problem_alloc(jmm_3d_wedge_problem_s **wedge) {
  *wedge = malloc(sizeof(jmm_3d_wedge_problem_s));
}

jmm_error_e jmm_3d_wedge_problem_init(jmm_3d_wedge_problem_s *wedge,
                                      jmm_3d_wedge_spec_s const *spec)
{
  wedge->spec = *spec; // just copy the spec over

  jmm_error_e error = JMM_ERROR_NONE;

  mesh3_alloc(&wedge->mesh);

  error = mesh3_init_from_3d_wedge_spec(wedge->mesh, spec);
  if (error != JMM_ERROR_NONE) {
    mesh3_dealloc(&wedge->mesh);
    return error;
  }

  // we allocate these here, but we'll initialize them later in solve
  eik3_alloc(&wedge->eik_direct);
  eik3_alloc(&wedge->eik_o_refl);
  eik3_alloc(&wedge->eik_n_refl);

  wedge->jet_direct_gt = malloc(mesh3_nverts(wedge->mesh)*sizeof(jet32t));

  return error;
}

void jmm_3d_wedge_problem_deinit(jmm_3d_wedge_problem_s *wedge) {
  mesh3_deinit(wedge->mesh);
  mesh3_dealloc(&wedge->mesh);

  if (eik3_is_initialized(wedge->eik_direct))
    eik3_deinit(wedge->eik_direct);

  if (eik3_is_initialized(wedge->eik_o_refl))
    eik3_deinit(wedge->eik_o_refl);

  if (eik3_is_initialized(wedge->eik_n_refl))
    eik3_deinit(wedge->eik_n_refl);

  eik3_dealloc(&wedge->eik_direct);
  eik3_dealloc(&wedge->eik_o_refl);
  eik3_dealloc(&wedge->eik_n_refl);

  free(wedge->jet_direct_gt);
  wedge->jet_direct_gt = NULL;
}

void jmm_3d_wedge_problem_dealloc(jmm_3d_wedge_problem_s **wedge) {
  free(*wedge);
  *wedge = NULL;
}

typedef struct {
  dbl3 xsrc, v1, x, t_out;
  dbl F;
} F_context;

static dbl dFdt(dbl t, F_context *context) {
  dbl3 xt_minus_xsrc;
  for (size_t i = 0; i < 3; ++i)
    xt_minus_xsrc[i] = t*context->v1[i] - context->xsrc[i];

  dbl3 x_minus_xt;
  for (size_t i = 0; i < 3; ++i)
    x_minus_xt[i] = context->x[i] - t*context->v1[i];

  dbl x_minus_xt_norm = dbl3_norm(x_minus_xt);
  dbl xt_minus_xsrc_norm = dbl3_norm(xt_minus_xsrc);

  context->F = x_minus_xt_norm + xt_minus_xsrc_norm;

  dbl3_normalized(x_minus_xt, context->t_out);

  dbl lhs = x_minus_xt_norm*dbl3_dot(context->v1, xt_minus_xsrc);
  dbl rhs = xt_minus_xsrc_norm*dbl3_dot(context->v1, x_minus_xt);

  return lhs - rhs;
}

static void set_jet_direct_gt(jmm_3d_wedge_problem_s *wedge, dbl sp, dbl phip) {
  size_t nverts = mesh3_nverts(wedge->mesh);

  F_context context = {
    .xsrc = {sp*cos(-phip), sp*sin(-phip), 0},
    .v1 = {0, 0, 1},
    .x = {NAN, NAN, NAN}
  };

  jmm_3d_wedge_spec_s const *spec = &wedge->spec;

  dbl tmax = spec->h/2, tmin = -spec->h/2, t;

  jet32t *jet = NULL;

  for (size_t i = 0; i < nverts; ++i) {
    mesh3_copy_vert(wedge->mesh, i, context.x);

    /* Get a pointer to the current jet: */
    jet = &wedge->jet_direct_gt[i];

    /* Compute the cylindrical angle of x about the wedge in order to
     * determine visibility. */
    dbl phi = atan2(context.x[1], context.x[0]);

    /* Compute the radius of x in cylindrical coordinates. */
    dbl rho = hypot(context.x[0], context.x[1]);

    /* The target point is directly visible: */
    if (rho == 0 || !(0 < phi && phi < JMM_PI - phip)) {
      /* Compute the eikonal and its gradient: */
      dbl3_sub(context.x, context.xsrc, jet->Df);
      jet->f = dbl3_norm(jet->Df);
      dbl3_dbl_div_inplace(jet->Df, jet->f);

      /* And its Hessian: */
      dbl33 outer;
      dbl3_outer(jet->Df, jet->Df, outer);
      dbl33_eye(jet->D2f);
      dbl33_sub_inplace(jet->D2f, outer);
      dbl33_dbl_div_inplace(jet->D2f, jet->f);
    }

    /* Edge diffraction has occurred: */
    else {
      if (!hybrid((hybrid_cost_func_t)dFdt, tmin, tmax, &context, &t))
        assert(false); // TODO: ?!

      /* Copy over eikonal value and gradient: */
      jet->f = context.F;
      dbl3_copy(context.t_out, jet->Df);

      /* Compute unit vector pointing from diffracting edge to x: */
      dbl3 v2 = {context.x[0]/rho, context.x[1]/rho, 0};

      /* Compute a unit vector orthogonal to v1 and v2: */
      dbl3 q1;
      dbl3_cross(context.v1, v2, q1);

      /* Compute au nit vector orthogonal to q1 and the ray
       * direction: */
      dbl3 q2;
      dbl3_cross(jet->Df, q1, q2);

      /* Compute the first curvature outer product: */
      dbl33 outer1;
      dbl3_outer(q1, q1, outer1);
      dbl33_dbl_div_inplace(outer1, rho);

      /* ... and the second curvature outer product: */
      dbl33 outer2;
      dbl3_outer(q2, q2, outer2);
      dbl33_dbl_div_inplace(outer2, jet->f);

      /* Sum them up to get the Hessian of the eikonal: */
      dbl33_add(outer1, outer2, jet->D2f);
    }
  }
}

jmm_error_e
jmm_3d_wedge_problem_solve(jmm_3d_wedge_problem_s *wedge, dbl sp, dbl phip,
                           dbl rfac, double omega) {
  (void)omega;

  dbl3 xsrc; // the point source
  dbl3 x; // a varying mesh vertex

  /* Reset the eikonal problems if we've called solve before */

  if (eik3_is_initialized(wedge->eik_direct))
    eik3_deinit(wedge->eik_direct);

  if (eik3_is_initialized(wedge->eik_n_refl))
    eik3_deinit(wedge->eik_n_refl);

  if (eik3_is_initialized(wedge->eik_o_refl))
    eik3_deinit(wedge->eik_o_refl);

  /* Compute the location of the point source */

  xsrc[0] = sp*cos(phip);
  xsrc[1] = -sp*sin(phip);
  xsrc[2] = 0;

  /* Set up and solve the direct eikonal problem */

  eik3_init(wedge->eik_direct, wedge->mesh, FTYPE_POINT_SOURCE);

  /* Find all vertices which lie inside the factoring radius, and
   * initialize the vertices of all cells incident on those vertices
   * with exact data. Note that there may be some cells which
   * intersect the factoring radius, but which aren't discovered this
   * way. */
  size_t num_initialized = 0;
  for (size_t i = 0; i < mesh3_ncells(wedge->mesh); ++i) {
    size_t cv[4];
    mesh3_cv(wedge->mesh, i, cv);

    /* first, check if any of the current cell's vertices intersect
     * the initialization ball... */
    bool found_intersection = false;
    for (size_t j = 0; j < 4; ++j) {
      mesh3_copy_vert(wedge->mesh, cv[j], x);
      if (dbl3_dist(x, xsrc) <= rfac) {
        found_intersection = true;
        break;
      }
    }
    if (!found_intersection)
      continue;

    /* ... and if they do, initialize all of the vertices in the cell
     * with exact data */
    for (size_t j = 0; j < 4; ++j) {
      if (eik3_is_trial(wedge->eik_direct, cv[j]))
        continue;

      jet32t jet;

      mesh3_copy_vert(wedge->mesh, cv[j], x);
      jet.f = dbl3_dist(x, xsrc);

      dbl3_sub(x, xsrc, jet.Df);
      dbl3_dbl_div_inplace(jet.Df, jet.f);

      dbl33 eye, Df_otimes_Df;
      dbl33_eye(eye);
      dbl3_outer(jet.Df, jet.Df, Df_otimes_Df);
      dbl33_sub(eye, Df_otimes_Df, jet.D2f);
      dbl33_dbl_div_inplace(jet.D2f, jet.f);

      eik3_add_pt_src_BCs(wedge->eik_direct, cv[j], jet);
      ++num_initialized;
    }
  }

  if (num_initialized == 0) {
    printf("No TRIAL vertices!\n");
    return JMM_ERROR_BAD_ARGUMENTS;
  }

  if (wedge->spec.verbose) {
    printf("Number of vertices in the initialization ball: %lu\n",
           num_initialized);
  }

  eik3_solve(wedge->eik_direct);

  /* Compute groundtruth data: */

  set_jet_direct_gt(wedge, sp, phip);

  return JMM_ERROR_NONE;
}

static void
jmm_3d_wedge_problem_dump_direct_jet_gt(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->jet_direct_gt, sizeof(wedge->jet_direct_gt[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

void jmm_3d_wedge_problem_dump(jmm_3d_wedge_problem_s *wedge,
                               char const *path,
                               bool dump_direct,
                               bool dump_o_face,
                               bool dump_n_face)
{
  size_t file_path_strlen = strlen(path) + 64;
  char *file_path = malloc(file_path_strlen + 1);

  /* Dump the wedge problem specification: */

  strcpy(file_path, path);
  file_path = strcat(file_path, "/spec.txt");
  jmm_3d_wedge_spec_dump(&wedge->spec, file_path);

  /* Dump the domain tetrahedron mesh's vertices and cell indices: */

  strcpy(file_path, path);
  file_path = strcat(file_path, "/verts.bin");
  mesh3_dump_verts(wedge->mesh, file_path);

  strcpy(file_path, path);
  file_path = strcat(file_path, "/cells.bin");
  mesh3_dump_cells(wedge->mesh, file_path);

  /* Create and dump the surface mesh's vertices and faces indices: */

  mesh2_s *surface_mesh = mesh3_get_surface_mesh(wedge->mesh);

  strcpy(file_path, path);
  file_path = strcat(file_path, "/surface_verts.bin");
  mesh2_dump_verts(surface_mesh, file_path);

  strcpy(file_path, path);
  file_path = strcat(file_path, "/surface_faces.bin");
  mesh2_dump_faces(surface_mesh, file_path);

  mesh2_deinit(surface_mesh);
  mesh2_dealloc(&surface_mesh);

  if (dump_direct) {

    /* Dump the direct eikonal's data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_jet.bin");
    eik3_dump_jet(wedge->eik_direct, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_state.bin");
    eik3_dump_state(wedge->eik_direct, file_path);

    /* Dump the direct eikonal's groundtruth data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_jet_gt.bin");
    jmm_3d_wedge_problem_dump_direct_jet_gt(wedge, file_path);

  }

  if (dump_o_face) {

    /* Dump the o-face reflected eikonal's data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_jet.bin");
    eik3_dump_jet(wedge->eik_o_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_state.bin");
    eik3_dump_state(wedge->eik_o_refl, file_path);

  }

  if (dump_n_face) {

    /* Dump the n-face reflected eikonal's data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_jet.bin");
    eik3_dump_jet(wedge->eik_n_refl, file_path);

  }

  /* Clean up: */

  free(file_path);
}
