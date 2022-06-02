#include "3d_wedge.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <array.h>
#include <bmesh.h>
#include <eik3_transport.h>
#include <error.h>
#include <hybrid.h>
#include <macros.h>
#include <mat.h>
#include <mesh2.h>

#include "mesh3_extra.h"

typedef enum wedge_eik {
  WEDGE_EIK_DIRECT,
  WEDGE_EIK_O_REFL,
  WEDGE_EIK_N_REFL,
} wedge_eik_e;

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

  wedge->D2T_direct = malloc(nverts*sizeof(dbl33));
  wedge->D2T_o_refl = malloc(nverts*sizeof(dbl33));
  wedge->D2T_n_refl = malloc(nverts*sizeof(dbl33));

  wedge->A_direct = malloc(nverts*sizeof(dbl));
  wedge->A_o_refl = malloc(nverts*sizeof(dbl));
  wedge->A_n_refl = malloc(nverts*sizeof(dbl));

  wedge->jet_direct_gt = malloc(nverts*sizeof(jet32t));
  wedge->jet_o_refl_gt = malloc(nverts*sizeof(jet32t));
  wedge->jet_n_refl_gt = malloc(nverts*sizeof(jet32t));

  wedge->t_in_direct = malloc(nverts*sizeof(dbl3));
  wedge->t_in_o_refl = malloc(nverts*sizeof(dbl3));
  wedge->t_in_n_refl = malloc(nverts*sizeof(dbl3));

  wedge->t_out_direct = malloc(nverts*sizeof(dbl3));
  wedge->t_out_o_refl = malloc(nverts*sizeof(dbl3));
  wedge->t_out_n_refl = malloc(nverts*sizeof(dbl3));

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

  free(wedge->D2T_direct);
  wedge->D2T_direct = NULL;

  free(wedge->D2T_o_refl);
  wedge->D2T_o_refl = NULL;

  free(wedge->D2T_n_refl);
  wedge->D2T_n_refl = NULL;

  free(wedge->A_direct);
  wedge->A_direct = NULL;

  free(wedge->A_o_refl);
  wedge->A_o_refl = NULL;

  free(wedge->A_n_refl);
  wedge->A_n_refl = NULL;

  free(wedge->t_in_direct);
  wedge->t_in_direct = NULL;

  free(wedge->t_in_o_refl);
  wedge->t_in_o_refl = NULL;

  free(wedge->t_in_n_refl);
  wedge->t_in_n_refl = NULL;

  free(wedge->t_out_direct);
  wedge->t_out_direct = NULL;

  free(wedge->t_out_o_refl);
  wedge->t_out_o_refl = NULL;

  free(wedge->t_out_n_refl);
  wedge->t_out_n_refl = NULL;
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
    .xsrc = {sp*cos(phip), sp*sin(phip), 0},
    .v1 = {0, 0, 1},
    .x = {NAN, NAN, NAN}
  };
}

bool in_valid_zone_direct(dbl rho, dbl phi, dbl phip, dbl n) {
  (void)n;
  return fabs(rho) < 1e-13 ||
    (-(2 - n)*JMM_PI/2 < phi && phi < JMM_PI + phip);
}

F_context get_context_o_refl(dbl sp, dbl phip, dbl n) {
  (void)n;
  return (F_context) {
    .xsrc = {sp*cos(phip), -sp*sin(phip), 0},
    .v1 = {0, 0, 1},
    .x = {NAN, NAN, NAN}
  };
}

bool in_valid_zone_o_refl(dbl rho, dbl phi, dbl phip, dbl n) {
  return fabs(rho) < 1e-13 ||
    (-(2 - n)*JMM_PI/2 <= phi && phi <= JMM_PI - phip);
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

bool in_valid_zone_n_refl(dbl rho, dbl phi, dbl phip, dbl n) {
  return fabs(rho) < 1e-13 ||
    !(-(2 - n)*JMM_PI/2 < phi && phi < (2*n - 1)*JMM_PI - phip);
}

dbl get_phi(dbl3 const x) {
  dbl phi = atan2(x[1], x[0]);
  return phi < 0 ? phi + 2*JMM_PI : phi;
}

static void set_jet_gt(jmm_3d_wedge_problem_s *wedge, dbl sp, dbl phip,
                       jet32t *jet,
                       F_context (*get_context)(dbl, dbl, dbl),
                       bool (*in_valid_zone)(dbl, dbl, dbl, dbl)) {
  size_t nverts = mesh3_nverts(wedge->mesh);

  F_context context = get_context(sp, phip, wedge->spec.n);

  dbl tmin = -wedge->spec.h/2;
  dbl tmax = wedge->spec.h/2;

  for (size_t i = 0; i < nverts; ++i) {
    mesh3_copy_vert(wedge->mesh, i, context.x);

    /* Compute the cylindrical angle of x about the wedge in order to
     * determine visibility. */
    dbl phi = get_phi(context.x);

    /* Compute the radius of x in cylindrical coordinates. */
    dbl rho = hypot(context.x[0], context.x[1]);

    /* The target point is in the valid zone: */
    if (in_valid_zone(rho, phi, phip, wedge->spec.n)) {
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

static bool updated_from_diff_edge(eik3_s const *eik, size_t l) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  par3_s par = eik3_get_par(eik, l);

  size_t npar = 0;
  for (size_t i = 0; i < 3; ++i)
    npar += par.l[i] != NO_PARENT;

  if (npar == 0 || npar == 3)
    return false;
  else if (npar == 1)
    return mesh3_vert_incident_on_diff_edge(mesh, par.l[0]);
  else /* npar == 2 */
    return mesh3_is_diff_edge(mesh, par.l);
}

static bool
any_cell_vert_updated_from_diff_edge(eik3_s const *eik, size_t cv[4]) {
  for (size_t i = 0; i < 4; ++i)
    if (updated_from_diff_edge(eik, cv[i]))
      return true;
  return false;
}

static bool cell_incident_on_diff_edge(mesh3_s const *mesh, size_t cv[4]) {
  for (size_t i = 0; i < 4; ++i)
    if (mesh3_vert_incident_on_diff_edge(mesh, cv[i]))
      return true;
  return false;
}

static bool par_inc_on_diff_edge(eik3_s const *eik, size_t l) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  par3_s par = eik3_get_par(eik, l);
  for (size_t i = 0; i < 3; ++i)
    if (par.l[i] != NO_PARENT
        && par.b[i] > 1e-14
        && mesh3_vert_incident_on_diff_edge(mesh, par.l[i]))
      return true;
  return false;
}

/* Approximate the Hessian at each vertex, storing the result for
 * vertex `l` at `D2T[l]`. The user should have already allocated and
 * initialized `D2T`. Entries which are `NAN` which will be filled,
 * and those which are finite will be left alone and used to compute
 * other values. */
static void approx_D2T(eik3_s const *eik, dbl33 *D2T) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  jet31t const *jet = eik3_get_jet_ptr(eik);

  dbl33 *D2T_cell = malloc(4*mesh3_ncells(mesh)*sizeof(dbl33));

  bool *has_init = malloc(mesh3_nverts(mesh)*sizeof(bool));
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    has_init[l] = dbl33_isfinite(D2T[l]);

  /* first, compute the Hessian at each cell vertex */
  for (size_t lc = 0, lv[4]; lc < mesh3_ncells(mesh); ++lc) {
    mesh3_cv(mesh, lc, lv);

    /* copy in initial values of D2T */
    for (size_t i = 0; i < 4; ++i)
      if (has_init[lv[i]])
        dbl33_copy(D2T[lv[i]], D2T_cell[4*lc + i]);

    /* get T and DT */
    jet31t J[4];
    for (size_t i = 0; i < 4; ++i) {
      J[i] = jet[lv[i]];
    }

    /* set up A */
    dbl4 A[3];
    for (size_t i = 0; i < 3; ++i) {
      dbl4_zero(A[i]);
      A[i][i] = 1;
      A[i][3] = -1;
    }

    /* get cell verts */
    dbl43 X;
    for (size_t i = 0; i < 4; ++i)
      mesh3_copy_vert(mesh, lv[i], X[i]);

    /* set up dX */
    dbl33 dX;
    for (size_t i = 0; i < 3; ++i)
      dbl3_sub(X[i], X[3], dX[i]);

    dbl33 dXinv, dXinvT;
    dbl33_copy(dX, dXinv);
    dbl33_invert(dXinv);
    dbl33_transposed(dXinv, dXinvT);

    /* set up bb33 */
    bb33 bb;
    bb33_init_from_jets(&bb, J, X);

    /* compute Hessian at each vertex */
    for (size_t i = 0; i < 4; ++i) {
      if (has_init[lv[i]])
        continue;

      dbl4 b;
      dbl4_e(b, i);

      /* compute Hessian in affine coordinates */
      dbl33 D2T_affine;
      for (size_t p = 0; p < 3; ++p) {
        for (size_t q = 0; q < 3; ++q) {
          dbl4 a[2];
          dbl4_copy(A[p], a[0]); // blech
          dbl4_copy(A[q], a[1]); // blech
          D2T_affine[p][q] = bb33_d2f(&bb, b, a);
        }
      }

      /* transform back to Cartesian and store with cell vertex */
      dbl33 tmp;
      dbl33_mul(dXinv, D2T_affine, tmp);
      dbl33_mul(tmp, dXinvT, D2T_cell[4*lc + i]);
    }
  }

  /* zero out D2T */
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    if (!has_init[l])
      dbl33_zero(D2T[l]);

  /* number of terms in weighted average */
  size_t *N = calloc(mesh3_nverts(mesh), sizeof(size_t));

  /* accumulate each D2T_cell entry into D2T */
  for (size_t lc = 0, cv[4]; lc < mesh3_ncells(mesh); ++lc) {
    mesh3_cv(mesh, lc, cv);

    /* skip this cell if its data is invalid */
    if (dbl33_isnan(D2T_cell[4*lc]))
      continue;

    for (size_t i = 0; i < 4; ++i) {
      if (has_init[cv[i]])
        continue;

      /* If this vertex was updated from a diff edge, don't use data
       * from a cell which is incident on a diff edge... */
      if (updated_from_diff_edge(eik, cv[i]) &&
          cell_incident_on_diff_edge(mesh, cv))
        continue;

      /* If this vertex is incident on a diff edge, don't use data
       * from a cell which was updated from a diff edge */
      if (mesh3_vert_incident_on_diff_edge(mesh, cv[i]) &&
          any_cell_vert_updated_from_diff_edge(eik, cv))
        continue;

      dbl33_add_inplace(D2T[cv[i]], D2T_cell[4*lc + i]);
      ++N[cv[i]]; /* increment number of terms in average */
    }
  }

  /* normalize each entry by the number of incident cells */
  for (size_t lv = 0; lv < mesh3_nverts(mesh); ++lv) {
    if (has_init[lv])
      continue;
    size_t nvc = mesh3_nvc(mesh, lv);
    dbl33_dbl_div_inplace(D2T[lv], nvc);
  }

  free(N);
  free(has_init);
  free(D2T_cell);
}

static eik3_s *get_eik(jmm_3d_wedge_problem_s *wedge, wedge_eik_e wedge_eik) {
  if      (wedge_eik == WEDGE_EIK_DIRECT) return wedge->eik_direct;
  else if (wedge_eik == WEDGE_EIK_O_REFL) return wedge->eik_o_refl;
  else if (wedge_eik == WEDGE_EIK_N_REFL) return wedge->eik_n_refl;
  else assert(false);
}

static dbl *get_A(jmm_3d_wedge_problem_s *wedge, wedge_eik_e wedge_eik) {
  if      (wedge_eik == WEDGE_EIK_DIRECT) return wedge->A_direct;
  else if (wedge_eik == WEDGE_EIK_O_REFL) return wedge->A_o_refl;
  else if (wedge_eik == WEDGE_EIK_N_REFL) return wedge->A_n_refl;
  else assert(false);
}

static dbl33 *get_D2T(jmm_3d_wedge_problem_s *wedge, wedge_eik_e wedge_eik) {
  if      (wedge_eik == WEDGE_EIK_DIRECT) return wedge->D2T_direct;
  else if (wedge_eik == WEDGE_EIK_O_REFL) return wedge->D2T_o_refl;
  else if (wedge_eik == WEDGE_EIK_N_REFL) return wedge->D2T_n_refl;
  else assert(false);
}

static size_t const *get_accepted_ptr(jmm_3d_wedge_problem_s *wedge, wedge_eik_e wedge_eik) {
  if      (wedge_eik == WEDGE_EIK_DIRECT) return eik3_get_accepted_ptr(wedge->eik_direct);
  else if (wedge_eik == WEDGE_EIK_O_REFL) return eik3_get_accepted_ptr(wedge->eik_o_refl);
  else if (wedge_eik == WEDGE_EIK_N_REFL) return eik3_get_accepted_ptr(wedge->eik_n_refl);
  else assert(false);
}

static void prop_amp(jmm_3d_wedge_problem_s *wedge, wedge_eik_e wedge_eik,
                     dbl sp, dbl phip) {
  size_t const *accepted = get_accepted_ptr(wedge, wedge_eik);

  eik3_s const *eik = get_eik(wedge, wedge_eik);
  dbl *A = get_A(wedge, wedge_eik);
  dbl33 *D2T = get_D2T(wedge, wedge_eik);

  dbl3 x, xsrc = {sp*cos(phip), sp*sin(phip), 0};

  for (size_t i = 0, l; i < mesh3_nverts(wedge->mesh); ++i) {
    l = accepted[i];

    par3_s par = eik3_get_par(eik, l);

    /* initialize amplitude BCs */
    mesh3_copy_vert(wedge->mesh, l, x);
    dbl r = dbl3_dist(x, xsrc);
    if (dbl3_all_nan(par.b)) {
      if (wedge_eik == WEDGE_EIK_DIRECT) {
        /* point source initialization */
        A[l] = 1/r;
      } else {
        /* copy data to initialize a reflection */
        A[l] = wedge->A_direct[l];
      }
      continue;
    }

    /* transport */

    dbl3 lam, abslam;
    size_t perm[3];
    dbl33_eigvals_sym(D2T[l], lam);
    dbl3_abs(lam, abslam);
    dbl3_argsort(abslam, perm);

    dbl kappa1 = lam[perm[2]], kappa2 = lam[perm[1]];

    dbl Alam = 1;
    assert(isfinite(par.b[0]));
    for (size_t j = 0; j < 3; ++j) {
      if (isfinite(par.b[j])) {
        assert(isfinite(A[par.l[j]]));
        Alam *= pow(A[par.l[j]], par.b[j]);
      }
    }

    dbl3 xlam = {0, 0, 0};
    for (size_t j = 0; j < 3; ++j) {
      if (isfinite(par.b[j])) {
        dbl3 x_;
        mesh3_copy_vert(wedge->mesh, par.l[j], x_);
        for (size_t k = 0; k < 3; ++k)
          xlam[k] += par.b[j]*x_[k];
      }
    }
    dbl L = dbl3_dist(x, xlam);

    A[l] = Alam*exp(-L*(kappa1 + kappa2)/2);
  }
}

/* Find the index of the reflector corresponding to the o-face. */
static size_t get_o_face_index(mesh3_s const *mesh) {
  dbl mesh_eps = mesh3_get_eps(mesh);
  dbl3 x;

  size_t o_face_index, nf, (*lf)[3] = NULL;

  for (o_face_index = 0;
       o_face_index < mesh3_get_num_reflectors(mesh);
       ++o_face_index) {
    nf = mesh3_get_reflector_size(mesh, o_face_index);
    lf = malloc(nf*sizeof(size_t[3]));
    mesh3_get_reflector(mesh, o_face_index, lf);

    /* Check if the y coordinate of any of the face vertices is
     * something other than zero... this is a sufficient check for the
     * o-face with this geometry */
    bool wrong_reflector = false;
    for (size_t i = 0; i < nf; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        mesh3_copy_vert(mesh, lf[i][j], x);
        wrong_reflector |= fabs(x[1]) > mesh_eps;
      }
    }

    /* Break if we found the o-face... and don't free the faces! */
    if (!wrong_reflector)
      break;

    /* Free the faces and move to the next reflector */
    free(lf);
    lf = NULL;
  }
  assert(lf != NULL);
  assert(o_face_index < mesh3_get_num_reflectors(mesh));

  return o_face_index;
}

jmm_error_e
jmm_3d_wedge_problem_solve(jmm_3d_wedge_problem_s *wedge, dbl sp, dbl phip,
                           dbl rfac, double omega) {
  (void)omega;

  dbl3 xsrc; // the point source
  dbl3 x; // a varying mesh vertex

  /* Compute the location of the point source */

  xsrc[0] = sp*cos(phip);
  xsrc[1] = sp*sin(phip);
  xsrc[2] = 0;

  /** DIRECT */

  /* Set up and solve the direct eikonal problem */

  eik3_add_pt_src_bcs(wedge->eik_direct, xsrc, rfac);

  array_s const *direct_bc_inds = eik3_get_bc_inds(wedge->eik_direct);

  /* Solve direct eikonal */
  if (wedge->spec.verbose) {
    printf("Number of vertices in the initialization ball: %lu\n",
           array_size(direct_bc_inds));
    printf("Solving point source problem... ");
    fflush(stdout);
  }
  eik3_solve(wedge->eik_direct);
  if (wedge->spec.verbose)
    puts("done");

  /* Compute groundtruth data for point source problem: */
  set_jet_gt(wedge, sp, phip, wedge->jet_direct_gt,
             get_context_direct, in_valid_zone_direct);
  if (wedge->spec.verbose)
    puts("- computed groundtruth data for direct arrival");

  /* Cell-averaged D2T for direct eikonal */

  if (wedge->spec.verbose) {
    printf("- computing D2T using cell averaging... ");
    fflush(stdout);
  }

  /* specially initialize D2T near the point source for the direct
   * eikonal */
  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l) {
    mesh3_copy_vert(wedge->mesh, l, x);
    if (dbl3_dist(x, xsrc) <= rfac)
      dbl33_copy(wedge->jet_direct_gt[l].D2f, wedge->D2T_direct[l]);
    else
      dbl33_nan(wedge->D2T_direct[l]);
  }

  /* we also want to initialize D2T for points which are immediately
   * downwind of the diffracting edge */
  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l) {
    par3_s par = eik3_get_par(wedge->eik_direct, l);
    size_t la[3];
    size_t na = par3_get_active_inds(&par, la);
    if (na == 2 && mesh3_is_diff_edge(wedge->mesh, la)) {
      /* Get the update point */
      dbl3 xhat;
      mesh3_copy_vert(wedge->mesh, l, xhat);

      /* Find the point of diffraction */
      dbl3 xe = {0, 0, xhat[2]}; // manually project...

      /* Unit tangent vector for diffracting edge */
      dbl3 te = {0, 0, 1};

      /* Compute unit vector pointing from `xe` to `xhat` */
      dbl3 tf;
      dbl3_sub(xhat, xe, tf);
      dbl rho = dbl3_normalize(tf); // cylindrical radius for `xhat`

      /* Get eikonal jet at `xhat` */
      jet31t J = eik3_get_jet(wedge->eik_direct, l);

      /* Get the ray direction */
      /* Compute unit vector orthogonal to `te` and `tf` (this vector
       * will be orthogonal to the ray direction) */
      dbl3 q1;
      dbl3_cross(te, tf, q1);
      assert(fabs(dbl3_dot(q1, J.Df)) < 1e-13);

      /* Compute unit vector orthogonal to `q1` and the ray
       * direction */
      dbl3 q2;
      dbl3_cross(J.Df, q1, q2);

      /* Compute first curvature outer product */
      dbl33 outer1;
      dbl3_outer(q1, q1, outer1);
      dbl33_dbl_div_inplace(outer1, rho);

      /* Compute second curvature outer product */
      dbl33 outer2;
      dbl3_outer(q2, q2, outer2);
      dbl33_dbl_div_inplace(outer2, J.f);

      /* Sum them up to get the Hessian */
      dbl33_add(outer1, outer2, wedge->D2T_direct[l]);
    }
  }

  approx_D2T(wedge->eik_direct, wedge->D2T_direct);

  if (wedge->spec.verbose) {
    puts("done");
    printf("- propagating the amplitude... ");
    fflush(stdout);
  }

  prop_amp(wedge, WEDGE_EIK_DIRECT, sp, phip);

  if (wedge->spec.verbose)
    puts("done");

  /** O-REFL */

  /* Figure out which reflector is the o-face */

  size_t o_face_refl_index = get_o_face_index(wedge->mesh);
  eik3_add_refl_bcs(wedge->eik_o_refl, wedge->eik_direct, o_face_refl_index);

  if (wedge->spec.verbose) {
    printf("Computing o-face reflection... ");
    fflush(stdout);
  }
  eik3_solve(wedge->eik_o_refl);
  if (wedge->spec.verbose)
    puts("done");

  /* Compute the groundtruth data for the o-face reflection: */
  set_jet_gt(wedge, sp, phip, wedge->jet_o_refl_gt,
             get_context_o_refl, in_valid_zone_o_refl);
  if (wedge->spec.verbose)
    puts("- computed groundtruth data for o-face reflection");

  /* Cell-averaged D2T for o-refl eikonal */
  if (wedge->spec.verbose) {
    printf("- computing D2T using cell averaging... ");
    fflush(stdout);
  }

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    dbl33_nan(wedge->D2T_o_refl[l]);

  approx_D2T(wedge->eik_o_refl, wedge->D2T_o_refl);

  if (wedge->spec.verbose) {
    puts("done");
    printf("- propagating the amplitude... ");
    fflush(stdout);
  }

  prop_amp(wedge, WEDGE_EIK_O_REFL, sp, phip);

  if (wedge->spec.verbose)
    puts("done");

  /** N-REFL */

  /* Ditto for the n-face problem: */
  dbl n_radians = JMM_PI*wedge->spec.n;

  /* Get the surface normal for the n-face */
  dbl3 n_normal = {-sin(n_radians), cos(n_radians), 0};

  /* Compute reflection matrix for the surface normal */
  dbl33 n_refl;
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      n_refl[i][j] = i == j ?
        1 - 2*n_normal[i]*n_normal[j] : -2*n_normal[i]*n_normal[j];

  // TODO: add diff or refl BCs depending on problem setup...
  // TODO: move this into add_*_BCs...

  // array_s const *n_refl_trial_inds = eik3_get_trial_inds(wedge->eik_n_refl);
  array_s const *n_refl_bc_inds = eik3_get_bc_inds(wedge->eik_n_refl);

  // TODO: set refl or diff BCs depending on angle...
  eik3_add_diff_bcs(wedge->eik_n_refl, wedge->eik_direct, 0);

  if (wedge->spec.verbose) {
    printf("Computing n-face reflection... ");
    fflush(stdout);
  }
  eik3_solve(wedge->eik_n_refl);
  if (wedge->spec.verbose)
    puts("done");

  /* Compute the groundtruth data for the n-face reflection: */
  set_jet_gt(wedge, sp, phip, wedge->jet_n_refl_gt,
             get_context_n_refl, in_valid_zone_n_refl);
  if (wedge->spec.verbose)
    puts("- computed groundtruth data for n-face reflection");

  /* cell-averaged D2T for n-face reflection */
  if (wedge->spec.verbose) {
    printf("- computing D2T using cell averaging... ");
    fflush(stdout);
  }

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    dbl33_nan(wedge->D2T_n_refl[l]);

  approx_D2T(wedge->eik_n_refl, wedge->D2T_n_refl);

  if (wedge->spec.verbose) {
    puts("done");
    printf("- propagating the amplitude... ");
    fflush(stdout);
  }

  prop_amp(wedge, WEDGE_EIK_N_REFL, sp, phip);

  if (wedge->spec.verbose)
    puts("done");

  puts("Finished solving eikonal equations");

  /** Transport t_in and t_out vectors: */

  jet31t const *jet = NULL;

  /* direct arrival */

  jet = eik3_get_jet_ptr(wedge->eik_direct);

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    dbl3_nan(wedge->t_in_direct[l]);

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    dbl3_nan(wedge->t_out_direct[l]);

  for (size_t i = 0, l; i < array_size(direct_bc_inds); ++i) {
    array_get(direct_bc_inds, i, &l);
    dbl3_copy(jet[l].Df, wedge->t_in_direct[l]);
    dbl3_copy(jet[l].Df, wedge->t_out_direct[l]);
  }

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    if (par_inc_on_diff_edge(wedge->eik_direct, l))
      dbl3_copy(jet[l].Df, wedge->t_out_direct[l]);

  eik3_transport_unit_vector(wedge->eik_direct, wedge->t_in_direct, true);
  eik3_transport_unit_vector(wedge->eik_direct, wedge->t_out_direct, true);

  /* o-refl */

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    dbl3_nan(wedge->t_in_o_refl[l]);

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    dbl3_nan(wedge->t_out_o_refl[l]);

  array_s const *o_refl_bc_inds = eik3_get_bc_inds(wedge->eik_o_refl);

  if (array_size(o_refl_bc_inds) > 0) {
    for (size_t i = 0, l; i < array_size(o_refl_bc_inds); ++i) {
      array_get(o_refl_bc_inds, i, &l);
      dbl3_copy(jet[l].Df, wedge->t_in_o_refl[l]);

      if (!mesh3_vert_incident_on_diff_edge(wedge->mesh, l)) {
        dbl3_copy(jet[l].Df, wedge->t_out_o_refl[l]);
        dbl33_dbl3_mul_inplace(n_refl, wedge->t_out_o_refl[l]);
      }
    }

    for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l) {
      if (updated_from_diff_edge(wedge->eik_o_refl, l)) {
        jet31t J = eik3_get_jet(wedge->eik_o_refl, l);
        dbl3_copy(J.Df, wedge->t_out_o_refl[l]);
      }
    }

    eik3_transport_unit_vector(wedge->eik_o_refl, wedge->t_in_o_refl, true);
    eik3_transport_unit_vector(wedge->eik_o_refl, wedge->t_out_o_refl, true);
  }

  /* n-refl */

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    dbl3_nan(wedge->t_in_n_refl[l]);

  for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
    dbl3_nan(wedge->t_out_n_refl[l]);

  if (array_size(n_refl_bc_inds) > 0) {
    for (size_t i = 0, l; i < array_size(n_refl_bc_inds); ++i) {
      array_get(n_refl_bc_inds, i, &l);
      dbl3_copy(jet[l].Df, wedge->t_in_n_refl[l]);

      if (!mesh3_vert_incident_on_diff_edge(wedge->mesh, l)) {
        dbl3_copy(jet[l].Df, wedge->t_out_n_refl[l]);
        dbl33_dbl3_mul_inplace(n_refl, wedge->t_out_n_refl[l]);
      }
    }

    for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l) {
      if (updated_from_diff_edge(wedge->eik_n_refl, l)) {
        jet31t J = eik3_get_jet(wedge->eik_n_refl, l);
        dbl3_copy(J.Df, wedge->t_out_n_refl[l]);
      }
    }

    eik3_transport_unit_vector(wedge->eik_n_refl, wedge->t_in_n_refl, true);
    eik3_transport_unit_vector(wedge->eik_n_refl, wedge->t_out_n_refl, true);
  }

  /** Clean up: */

  return JMM_ERROR_NONE;
}

static void jmm_3d_wedge_problem_dump_direct_hess(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->D2T_direct, sizeof(wedge->D2T_direct[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void jmm_3d_wedge_problem_dump_o_refl_hess(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->D2T_o_refl, sizeof(wedge->D2T_o_refl[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void jmm_3d_wedge_problem_dump_n_refl_hess(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->D2T_n_refl, sizeof(wedge->D2T_n_refl[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
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
jmm_3d_wedge_problem_dump_direct_t_in(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->t_in_direct, sizeof(wedge->t_in_direct[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void
jmm_3d_wedge_problem_dump_direct_t_out(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->t_out_direct, sizeof(wedge->t_out_direct[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void
jmm_3d_wedge_problem_dump_o_refl_t_in(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->t_in_o_refl, sizeof(wedge->t_in_o_refl[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void
jmm_3d_wedge_problem_dump_o_refl_t_out(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->t_out_o_refl, sizeof(wedge->t_out_o_refl[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void
jmm_3d_wedge_problem_dump_n_refl_t_in(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->t_in_n_refl, sizeof(wedge->t_in_n_refl[0]),
         mesh3_nverts(wedge->mesh), fp);
  fclose(fp);
}

static void
jmm_3d_wedge_problem_dump_n_refl_t_out(
  jmm_3d_wedge_problem_s const *wedge,
  char const *path)
{
  FILE *fp = fopen(path, "wb");
  fwrite(wedge->t_out_n_refl, sizeof(wedge->t_out_n_refl[0]),
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
    file_path = strcat(file_path, "/direct_hess.bin");
    jmm_3d_wedge_problem_dump_direct_hess(wedge, file_path);

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

    /* Dump the t_in and t_out vectors */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_t_in.bin");
    jmm_3d_wedge_problem_dump_direct_t_in(wedge, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/direct_t_out.bin");
    jmm_3d_wedge_problem_dump_direct_t_out(wedge, file_path);

  }

  if (dump_o_face) {

    /* Dump the o-face reflected eikonal's data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_jet.bin");
    eik3_dump_jet(wedge->eik_o_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_hess.bin");
    jmm_3d_wedge_problem_dump_o_refl_hess(wedge, file_path);

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

    /* Dump the t_in and t_out vectors */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_t_in.bin");
    jmm_3d_wedge_problem_dump_o_refl_t_in(wedge, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/o_refl_t_out.bin");
    jmm_3d_wedge_problem_dump_o_refl_t_out(wedge, file_path);

  }

  if (dump_n_face) {

    /* Dump the n-face reflected eikonal's data: */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_jet.bin");
    eik3_dump_jet(wedge->eik_n_refl, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_hess.bin");
    jmm_3d_wedge_problem_dump_n_refl_hess(wedge, file_path);

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

    /* Dump the t_in and t_out vectors */

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_t_in.bin");
    jmm_3d_wedge_problem_dump_n_refl_t_in(wedge, file_path);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/n_refl_t_out.bin");
    jmm_3d_wedge_problem_dump_n_refl_t_out(wedge, file_path);

  }

  /* Clean up: */

  free(file_path);
}

static void
find_cells_for_img_grid(mesh3_s const *mesh, grid2_s const *img_grid, dbl z,
                        size_t **lc_grid, size_t (**cv_grid)[4], dbl4 **b_grid)
{
  size_t num_grid = grid2_nind(img_grid);

  /* Allocate space for cell indices, incident vertices, and
   * barycentric coordinates */
  *lc_grid = malloc(num_grid*sizeof(size_t));
  *cv_grid = malloc(num_grid*sizeof(size_t[4]));
  *b_grid = malloc(num_grid*sizeof(dbl4));

  /* Initialize everything to a bad value */
  for (size_t i = 0; i < num_grid; ++i) {
    (*lc_grid)[i] = (size_t)NO_INDEX;
    for (size_t j = 0; j < 4; ++j) {
      (*cv_grid)[i][j] = (size_t)NO_INDEX;
      (*b_grid)[i][j] = NAN;
    }
  }

  for (size_t lc = 0; lc < mesh3_ncells(mesh); ++lc) {
    tetra3 tetra = mesh3_get_tetra(mesh, lc);

    size_t offset[2];
    grid2_s subgrid = tetra3_get_covering_xy_subgrid(&tetra, img_grid, offset);

    for (size_t l = 0; l < grid2_nind(&subgrid); ++l) {
      int ind[2];
      grid2_l2ind(&subgrid, l, ind);

      int ind_orig[2] = {offset[0] + ind[0], offset[1] + ind[1]};
      size_t l_orig = grid2_ind2l(img_grid, ind_orig);
      if ((*lc_grid)[l_orig] != (size_t)NO_INDEX)
        continue;

      dbl3 x = {[2] = z};
      grid2_l2xy(img_grid, l_orig, x);
      if (!mesh3_cell_contains_point(mesh, lc, x))
        continue;

      (*lc_grid)[l_orig] = lc;
      mesh3_cv(mesh, lc, (*cv_grid)[l_orig]);
      tetra3_get_bary_coords(&tetra, x, (*b_grid)[l_orig]);
    }
  }
}

static void
dump_slice(jmm_3d_wedge_problem_s const *wedge, grid2_s const *img_grid,
           size_t const *lc_grid, size_t const (*cv_grid)[4], dbl4 const *b_grid,
           char const *file_path, field_e field, wedge_eik_e wedge_eik) {
  FILE *fp = fopen(file_path, "wb");

  if (field == FIELD_A) {
    dbl const *A = get_A((jmm_3d_wedge_problem_s *)wedge, wedge_eik);
    dbl Ab;
    for (size_t l = 0; l < grid2_nind(img_grid); ++l) {
      if (isnan(b_grid[l][0])) {
        Ab = NAN;
      } else {
        Ab = 1;
        for (size_t i = 0; i < 4; ++i)
          Ab *= pow(A[cv_grid[l][i]], b_grid[l][i]);
      }
      fwrite(&Ab, sizeof(Ab), 1, fp);
    }
  } else {
    jet31t *jet = NULL;
    jet32t *jet_gt = NULL;
    jet31t *field_jet = NULL;

    if (field == FIELD_T || field == FIELD_E_T) {
      if (wedge_eik == WEDGE_EIK_DIRECT)
        jet = eik3_get_jet_ptr(wedge->eik_direct);
      else if (wedge_eik == WEDGE_EIK_O_REFL)
        jet = eik3_get_jet_ptr(wedge->eik_o_refl);
      else if (wedge_eik == WEDGE_EIK_N_REFL)
        jet = eik3_get_jet_ptr(wedge->eik_n_refl);
      else
        assert(false);
    }

    if (field == FIELD_E_T) {
      if (wedge_eik == WEDGE_EIK_DIRECT)
        jet_gt = wedge->jet_direct_gt;
      else if (wedge_eik == WEDGE_EIK_O_REFL)
        jet_gt = wedge->jet_o_refl_gt;
      else if (wedge_eik == WEDGE_EIK_N_REFL)
        jet_gt = wedge->jet_n_refl_gt;
      else
        assert(false);
    }

    if (field == FIELD_T)
      field_jet = jet;

    if (field == FIELD_E_T) {
      field_jet = malloc(mesh3_nverts(wedge->mesh)*sizeof(jet31t));
      for (size_t l = 0; l < mesh3_nverts(wedge->mesh); ++l)
        jet31t_sub((jet31t const *)&jet_gt[l], &jet[l], &field_jet[l]);
    }

    bmesh33_s *bmesh;
    bmesh33_alloc(&bmesh);
    bmesh33_init_from_mesh3_and_jets(bmesh, wedge->mesh, field_jet);

    for (size_t l = 0; l < grid2_nind(img_grid); ++l) {
      dbl value;
      if (isnan(b_grid[l][0]))
        value = NAN;
      else
        value = bb33_f(bmesh33_get_bb_ptr(bmesh, lc_grid[l]), b_grid[l]);
      fwrite(&value, sizeof(value), 1, fp);
    }

    bmesh33_deinit(bmesh);
    bmesh33_dealloc(&bmesh);

    if (field == FIELD_E_T)
      free(field_jet);
  }

  fclose(fp);
}

void jmm_3d_wedge_problem_save_slice_plots(jmm_3d_wedge_problem_s const *wedge,
                                           char const *path,
                                           bool dump_direct,
                                           bool dump_o_face,
                                           bool dump_n_face,
                                           grid2_s const *img_grid)
{
  size_t file_path_strlen = strlen(path) + 64;
  char *file_path = malloc(file_path_strlen + 1);

  if (dump_direct || dump_o_face || dump_n_face) {
    strcpy(file_path, path);
    file_path = strcat(file_path, "/img_grid.txt");
    grid2_save(img_grid, file_path);
  }

  size_t *lc_grid = NULL;
  size_t (*cv_grid)[4] = NULL;
  dbl4 *b_grid = NULL;
  find_cells_for_img_grid(wedge->mesh, img_grid, 0, &lc_grid, &cv_grid, &b_grid);

  if (dump_direct) {
    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_direct_T.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_T,
               WEDGE_EIK_DIRECT);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_direct_A.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_A,
               WEDGE_EIK_DIRECT);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_direct_E_T.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_E_T,
               WEDGE_EIK_DIRECT);
  }

  if (dump_o_face) {
    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_o_refl_T.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_T,
               WEDGE_EIK_O_REFL);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_o_refl_A.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_A,
               WEDGE_EIK_O_REFL);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_o_refl_E_T.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_E_T,
               WEDGE_EIK_O_REFL);
  }

  if (dump_n_face) {
    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_n_refl_T.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_T,
               WEDGE_EIK_N_REFL);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_n_refl_A.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_A,
               WEDGE_EIK_N_REFL);

    strcpy(file_path, path);
    file_path = strcat(file_path, "/slice_n_refl_E_T.bin");
    dump_slice(wedge, img_grid, lc_grid, cv_grid, b_grid, file_path, FIELD_E_T,
               WEDGE_EIK_N_REFL);
  }

  free(lc_grid);
  free(b_grid);
  free(cv_grid);

  free(file_path);
}
