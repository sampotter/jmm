#include "3d_wedge.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <array.h>
#include <eik3_transport.h>
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

  size_t nverts = mesh3_nverts(wedge->mesh);

  eik3_alloc(&wedge->eik_direct);
  eik3_init(wedge->eik_direct, wedge->mesh);

  eik3_alloc(&wedge->eik_o_refl);
  eik3_init(wedge->eik_o_refl, wedge->mesh);

  eik3_alloc(&wedge->eik_n_refl);
  eik3_init(wedge->eik_n_refl, wedge->mesh);

  wedge->jet_direct_gt = malloc(nverts*sizeof(jet32t));
  wedge->jet_o_refl_gt = malloc(nverts*sizeof(jet32t));
  wedge->jet_n_refl_gt = malloc(nverts*sizeof(jet32t));

  wedge->origin_direct = malloc(nverts*sizeof(dbl));
  wedge->origin_o_refl = malloc(nverts*sizeof(dbl));
  wedge->origin_n_refl = malloc(nverts*sizeof(dbl));

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

  free(wedge->jet_o_refl_gt);
  wedge->jet_o_refl_gt = NULL;

  free(wedge->jet_n_refl_gt);
  wedge->jet_n_refl_gt = NULL;

  free(wedge->origin_direct);
  wedge->origin_direct = NULL;

  free(wedge->origin_o_refl);
  wedge->origin_o_refl = NULL;

  free(wedge->origin_n_refl);
  wedge->origin_n_refl = NULL;
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

F_context get_context_direct(dbl sp, dbl phip, dbl n) {
  (void)n;
  return (F_context) {
    .xsrc = {sp*cos(-phip), sp*sin(-phip), 0},
    .v1 = {0, 0, 1},
    .x = {NAN, NAN, NAN}
  };
}

bool in_valid_zone_direct(dbl phi, dbl phip, dbl n) {
  (void)n;
  return !(0 < phi && phi < JMM_PI - phip);
}

F_context get_context_o_refl(dbl sp, dbl phip, dbl n) {
  (void)n;
  dbl phi_img = phip;
  return (F_context) {
    .xsrc = {sp*cos(phi_img), sp*sin(phi_img), 0},
    .v1 = {0, 0, 1},
    .x = {NAN, NAN, NAN}
  };
}

bool in_valid_zone_o_refl(dbl phi, dbl phip, dbl n) {
  dbl dphi = JMM_PI - phip;
  return -dphi <= phi && phi < n*JMM_PI/2;
}

F_context get_context_n_refl(dbl sp, dbl phip, dbl n) {
  dbl dphi = (n - 1)*JMM_PI + phip;
  dbl phi_img = phip - 2*dphi;
  return (F_context) {
    .xsrc = {sp*cos(-phi_img), sp*sin(-phi_img), 0},
    .v1 = {0, 0, 1},
    .x = {NAN, NAN, NAN}
  };
}

bool in_valid_zone_n_refl(dbl phi, dbl phip, dbl n) {
  dbl dphi = (n - 1)*JMM_PI + phip;
  return n*JMM_PI/2 < phi && 0 <= n*JMM_PI + dphi;
}

static void set_jet_gt(jmm_3d_wedge_problem_s *wedge, dbl sp, dbl phip,
                       jet32t *jet,
                       F_context (*get_context)(dbl, dbl, dbl),
                       bool (*in_valid_zone)(dbl, dbl, dbl)) {
  size_t nverts = mesh3_nverts(wedge->mesh);

  F_context context = get_context(sp, phip, wedge->spec.n);

  dbl tmin = -wedge->spec.h/2;
  dbl tmax = wedge->spec.h/2;

  for (size_t i = 0; i < nverts; ++i) {
    mesh3_copy_vert(wedge->mesh, i, context.x);

    /* Compute the cylindrical angle of x about the wedge in order to
     * determine visibility. */
    dbl phi = atan2(context.x[1], context.x[0]);

    /* Compute the radius of x in cylindrical coordinates. */
    dbl rho = hypot(context.x[0], context.x[1]);

    /* The target point is in the valid zone: */
    if (rho == 0 || in_valid_zone(phi, phip, wedge->spec.n)) {
      /* Compute the eikonal and its gradient: */
      dbl3_sub(context.x, context.xsrc, jet[i].Df);
      jet[i].f = dbl3_norm(jet[i].Df);
      dbl3_dbl_div_inplace(jet[i].Df, jet[i].f);

      /* And its Hessian: */
      dbl33 outer;
      dbl3_outer(jet[i].Df, jet[i].Df, outer);
      dbl33_eye(jet[i].D2f);
      dbl33_sub_inplace(jet[i].D2f, outer);
      dbl33_dbl_div_inplace(jet[i].D2f, jet[i].f);
    }

    /* Edge diffraction has occurred: */
    else {
      dbl t;
      if (!hybrid((hybrid_cost_func_t)dFdt, tmin, tmax, &context, &t))
        assert(false); // TODO: ?!

      /* Copy over eikonal value and gradient: */
      jet[i].f = context.F;
      dbl3_copy(context.t_out, jet[i].Df);

      /* Compute unit vector pointing from diffracting edge to x: */
      dbl3 v2 = {context.x[0]/rho, context.x[1]/rho, 0};

      /* Compute a unit vector orthogonal to v1 and v2: */
      dbl3 q1;
      dbl3_cross(context.v1, v2, q1);

      /* Compute au nit vector orthogonal to q1 and the ray
       * direction: */
      dbl3 q2;
      dbl3_cross(jet[i].Df, q1, q2);

      /* Compute the first curvature outer product: */
      dbl33 outer1;
      dbl3_outer(q1, q1, outer1);
      dbl33_dbl_div_inplace(outer1, rho);

      /* ... and the second curvature outer product: */
      dbl33 outer2;
      dbl3_outer(q2, q2, outer2);
      dbl33_dbl_div_inplace(outer2, jet[i].f);

      /* Sum them up to get the Hessian of the eikonal: */
      dbl33_add(outer1, outer2, jet[i].D2f);
    }
  }
}

jmm_error_e
jmm_3d_wedge_problem_solve(jmm_3d_wedge_problem_s *wedge, dbl sp, dbl phip,
                           dbl rfac, double omega) {
  (void)omega;

  dbl3 xsrc; // the point source
  dbl3 x; // a varying mesh vertex

  /* Compute the location of the point source */

  xsrc[0] = sp*cos(phip);
  xsrc[1] = -sp*sin(phip);
  xsrc[2] = 0;

  /* Set up and solve the direct eikonal problem */

  array_s *direct_trial_inds;
  array_alloc(&direct_trial_inds);
  array_init(direct_trial_inds, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Find all vertices which lie inside the factoring radius, and
   * initialize the vertices of all cells incident on those vertices
   * with exact data. Note that there may be some cells which
   * intersect the factoring radius, but which aren't discovered this
   * way. */
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

      jet31t jet;

      mesh3_copy_vert(wedge->mesh, cv[j], x);
      jet.f = dbl3_dist(x, xsrc);

      dbl3_sub(x, xsrc, jet.Df);
      dbl3_dbl_div_inplace(jet.Df, jet.f);

      eik3_add_trial(wedge->eik_direct, cv[j], jet);

      array_append(direct_trial_inds, &cv[j]);
    }
  }

  if (array_size(direct_trial_inds) == 0) {
    printf("No TRIAL vertices!\n");
    return JMM_ERROR_BAD_ARGUMENTS;
  }

  if (wedge->spec.verbose) {
    printf("Number of vertices in the initialization ball: %lu\n",
           array_size(direct_trial_inds));
    printf("Solving point source problem... ");
  }

  eik3_solve(wedge->eik_direct);

  if (wedge->spec.verbose)
    puts("done");

  /* Compute groundtruth data for point source problem: */
  set_jet_gt(wedge, sp, phip, wedge->jet_direct_gt,
             get_context_direct, in_valid_zone_direct);

  array_s *o_refl_trial_inds;
  array_alloc(&o_refl_trial_inds);
  array_init(o_refl_trial_inds, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* If the geometry is right for it, set up and solve the o-face
   * reflection eikonal problem: */
  if (-phip > -JMM_PI) {
    for (size_t i = 0; i < mesh3_nverts(wedge->mesh); ++i) {
      mesh3_copy_vert(wedge->mesh, i, x);

      if (x[0] < 0 || x[1] != 0)
        continue;

      jet31t jet = eik3_get_jet(wedge->eik_direct, i);

      /* Reflect gradient over the o-face */
      jet.Df[1] *= -1;

      eik3_add_trial(wedge->eik_o_refl, i, jet);

      /* Keep track of the TRIAL index for tracking origins later */
      array_append(o_refl_trial_inds, &i);
    }

    if (wedge->spec.verbose)
      printf("Computing o-face reflection... ");

    eik3_solve(wedge->eik_o_refl);

    if (wedge->spec.verbose)
      puts("done");

    /* Compute the groundtruth data for the o-face reflection: */
    set_jet_gt(wedge, sp, phip, wedge->jet_o_refl_gt,
               get_context_o_refl, in_valid_zone_o_refl);
    if (wedge->spec.verbose)
      puts("Computed groundtruth data for o-face reflection");
  }

  array_s *n_refl_trial_inds;
  array_alloc(&n_refl_trial_inds);
  array_init(n_refl_trial_inds, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Ditto for the n-face problem: */
  dbl n_radians = JMM_PI*wedge->spec.n;
  if (-phip < n_radians - JMM_PI) {
    /* Get the surface normal for the n-face */
    dbl3 n_normal = {-sin(n_radians), cos(n_radians), 0};

    /* Compute reflection matrix for the surface normal */
    dbl33 n_refl;
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        n_refl[i][j] = i == j ?
          1 - 2*n_normal[i]*n_normal[j] : -2*n_normal[i]*n_normal[j];

    for (size_t i = 0; i < mesh3_nverts(wedge->mesh); ++i) {
      /* Skip vertices that aren't on the boundary */
      if (!mesh3_bdv(wedge->mesh, i))
        continue;

      /* Check if the angle of the vertex matches the angle of the
       * n-face of the wedge */
      mesh3_copy_vert(wedge->mesh, i, x);
      if (fabs(atan2(x[1], x[0]) - n_radians) > 1e-7)
        continue;

      jet31t jet = eik3_get_jet(wedge->eik_direct, i);

      /* Reflected gradient over n-face */
      dbl33_dbl3_mul_inplace(n_refl, jet.Df);

      eik3_add_trial(wedge->eik_n_refl, i, jet);

      /* Keep track of the TRIAL index for tracking origins later */
      array_append(n_refl_trial_inds, &i);
    }

    if (wedge->spec.verbose)
      printf("Computing n-face reflection... ");

    eik3_solve(wedge->eik_n_refl);

    if (wedge->spec.verbose)
      puts("done");

    /* Compute the groundtruth data for the n-face reflection: */
    set_jet_gt(wedge, sp, phip, wedge->jet_n_refl_gt,
               get_context_n_refl, in_valid_zone_n_refl);
    if (wedge->spec.verbose)
      puts("Computed groundtruth data for n-face reflection");
  }

  puts("Finished solving eikonal equations");

  /** Compute "approximate origins": */

  /* direct arrival */

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    wedge->origin_direct[l] = NAN;

  for (size_t i = 0, l; i < array_size(direct_trial_inds); ++i) {
    array_get(direct_trial_inds, i, &l);
    wedge->origin_direct[l] = 1; /* initially `TRIAL` nodes originate
                                  * from the point source */
  }

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l) {
    if (!isnan(wedge->origin_direct[l]))
      continue; /* already set this value */

    mesh3_copy_vert(wedge->mesh, l, x);
    if (hypot(x[0], x[1]) < 1e-13)
      wedge->origin_direct[l] = 0; /* label nodes on the diff. edge */
  }

  eik3_transport_dbl(wedge->eik_direct, wedge->origin_direct, true);

  /* o-refl */

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    wedge->origin_o_refl[l] = NAN;

  if (array_size(o_refl_trial_inds) > 0) {
    for (size_t i = 0, l; i < array_size(o_refl_trial_inds); ++i) {
      array_get(o_refl_trial_inds, i, &l);
      wedge->origin_o_refl[l] = 1;
    }

    for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l) {
      mesh3_copy_vert(wedge->mesh, l, x);
      if (hypot(x[0], x[1]) < 1e-13) {
        assert(wedge->origin_o_refl[l] == 1);
        wedge->origin_o_refl[l] = 0;
      }
    }

    eik3_transport_dbl(wedge->eik_o_refl, wedge->origin_o_refl, true);
  }

  /* n-refl */

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    wedge->origin_n_refl[l] = NAN;

  if (array_size(n_refl_trial_inds) > 0) {
    for (size_t i = 0, l; i < array_size(n_refl_trial_inds); ++i) {
      array_get(n_refl_trial_inds, i, &l);
      wedge->origin_n_refl[l] = 1;
    }

    for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l) {
      mesh3_copy_vert(wedge->mesh, l, x);
      if (hypot(x[0], x[1]) < 1e-13) {
        assert(wedge->origin_n_refl[l] == 1);
        wedge->origin_n_refl[l] = 0;
      }
    }

    eik3_transport_dbl(wedge->eik_n_refl, wedge->origin_n_refl, true);
  }

  /** Clean up: */

  array_deinit(direct_trial_inds);
  array_dealloc(&direct_trial_inds);

  array_deinit(o_refl_trial_inds);
  array_dealloc(&o_refl_trial_inds);

  array_deinit(n_refl_trial_inds);
  array_dealloc(&n_refl_trial_inds);

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

static void
jmm_3d_wedge_problem_dump_direct_origin(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->origin_direct, sizeof(wedge->origin_direct[0]),
         mesh3_nverts(wedge->mesh), fp);
}

static void
jmm_3d_wedge_problem_dump_o_refl_jet_gt(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->jet_o_refl_gt, sizeof(wedge->jet_o_refl_gt[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void
jmm_3d_wedge_problem_dump_o_refl_origin(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->origin_o_refl, sizeof(wedge->origin_o_refl[0]),
         mesh3_nverts(wedge->mesh), fp);
}

static void
jmm_3d_wedge_problem_dump_n_refl_jet_gt(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->jet_n_refl_gt, sizeof(wedge->jet_n_refl_gt[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void
jmm_3d_wedge_problem_dump_n_refl_origin(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->origin_n_refl, sizeof(wedge->origin_n_refl[0]),
         mesh3_nverts(wedge->mesh), fp);
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

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_par_l.bin");
    eik3_dump_par_l(wedge->eik_direct, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_par_b.bin");
    eik3_dump_par_b(wedge->eik_direct, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_accepted.bin");
    eik3_dump_accepted(wedge->eik_direct, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_has_bc.bin");
    eik3_dump_has_bc(wedge->eik_direct, file_path);

    /* Dump the direct eikonal's groundtruth data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_jet_gt.bin");
    jmm_3d_wedge_problem_dump_direct_jet_gt(wedge, file_path);

    /* Dump "approximate origins" for direct eikonal: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_origin.bin");
    jmm_3d_wedge_problem_dump_direct_origin(wedge, file_path);

  }

  if (dump_o_face) {

    /* Dump the o-face reflected eikonal's data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_jet.bin");
    eik3_dump_jet(wedge->eik_o_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_state.bin");
    eik3_dump_state(wedge->eik_o_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_par_l.bin");
    eik3_dump_par_l(wedge->eik_o_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_par_b.bin");
    eik3_dump_par_b(wedge->eik_o_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_accepted.bin");
    eik3_dump_accepted(wedge->eik_o_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_has_bc.bin");
    eik3_dump_has_bc(wedge->eik_o_refl, file_path);

    /* Dump the o-face eikonal's groundtruth data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_jet_gt.bin");
    jmm_3d_wedge_problem_dump_o_refl_jet_gt(wedge, file_path);

    /* Dump "approximate origins" for the o-face eikonal: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_origin.bin");
    jmm_3d_wedge_problem_dump_o_refl_origin(wedge, file_path);

  }

  if (dump_n_face) {

    /* Dump the n-face reflected eikonal's data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_jet.bin");
    eik3_dump_jet(wedge->eik_n_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_state.bin");
    eik3_dump_state(wedge->eik_n_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_par_l.bin");
    eik3_dump_par_l(wedge->eik_n_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_par_b.bin");
    eik3_dump_par_b(wedge->eik_n_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_accepted.bin");
    eik3_dump_accepted(wedge->eik_n_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_has_bc.bin");
    eik3_dump_has_bc(wedge->eik_n_refl, file_path);

    /* Dump the n-face eikonal's groundtruth data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_jet_gt.bin");
    jmm_3d_wedge_problem_dump_n_refl_jet_gt(wedge, file_path);

    /* Dump "approximate origins" for the n-face eikonal: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_origin.bin");
    jmm_3d_wedge_problem_dump_n_refl_origin(wedge, file_path);

  }

  /* Clean up: */

  free(file_path);
}
