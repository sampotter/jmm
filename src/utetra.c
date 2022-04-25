#include "utetra.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "array.h"
#include "bb.h"
#include "eik3.h"
#include "macros.h"
#include "mat.h"
#include "mesh3.h"
#include "opt.h"
#include "util.h"

#define MAX_NITER 100

utetra_spec_s utetra_spec_empty() {
  return (utetra_spec_s) {
    .eik = NULL,
    .lhat = (size_t)NO_INDEX,
    .l = {NO_INDEX, NO_INDEX, NO_INDEX},
    .state = {UNKNOWN, UNKNOWN, UNKNOWN},
    .xhat = {NAN, NAN, NAN},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .Df = {NAN, NAN, NAN}},
      {.f = INFINITY, .Df = {NAN, NAN, NAN}},
      {.f = INFINITY, .Df = {NAN, NAN, NAN}}
    },
    .tol = NAN
  };
}

utetra_spec_s utetra_spec_from_eik_and_inds(eik3_s const *eik, size_t l,
                                            size_t l0, size_t l1, size_t l2) {
  dbl h = mesh3_get_min_edge_length(eik3_get_mesh(eik));

  state_e const *state = eik3_get_state_ptr(eik);
  return (utetra_spec_s) {
    .eik = eik,
    .lhat = l,
    .l = {l0, l1, l2},
    .state = {state[l0], state[l1], state[l2]},
    .xhat = {NAN, NAN, NAN},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .Df = {NAN, NAN, NAN}},
      {.f = INFINITY, .Df = {NAN, NAN, NAN}},
      {.f = INFINITY, .Df = {NAN, NAN, NAN}}
    },
    .tol = pow(h, 4)
  };
}

utetra_spec_s utetra_spec_from_eik_without_l(eik3_s const *eik, dbl const x[3],
                                             size_t l0, size_t l1, size_t l2) {
  dbl h = mesh3_get_min_edge_length(eik3_get_mesh(eik));

  state_e const *state = eik3_get_state_ptr(eik);
  return (utetra_spec_s) {
    .eik = eik,
    .lhat = (size_t)NO_INDEX,
    .l = {l0, l1, l2},
    .state = {state[l0], state[l1], state[l2]},
    .xhat = {x[0], x[1], x[2]},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .Df = {NAN, NAN, NAN}},
      {.f = INFINITY, .Df = {NAN, NAN, NAN}},
      {.f = INFINITY, .Df = {NAN, NAN, NAN}}
    },
    .tol = pow(h, 4)
  };
}

utetra_spec_s utetra_spec_from_ptrs(mesh3_s const *mesh, jet31t const *jet,
                                    size_t l, size_t l0, size_t l1, size_t l2) {
  dbl h = mesh3_get_min_edge_length(mesh);

  utetra_spec_s spec = {
    .eik = NULL,
    .lhat = NO_INDEX,
    .l = {NO_INDEX, NO_INDEX, NO_INDEX},
    .state = {UNKNOWN, UNKNOWN, UNKNOWN},
    .jet = {jet[l0], jet[l1], jet[l2]},
    .tol = pow(h, 4)
  };

  mesh3_copy_vert(mesh, l, spec.xhat);
  mesh3_copy_vert(mesh, l0, spec.x[0]);
  mesh3_copy_vert(mesh, l1, spec.x[1]);
  mesh3_copy_vert(mesh, l2, spec.x[2]);

  return spec;
}

struct utetra {
  dbl lam[2]; // Current iterate
  dbl f;
  dbl g[2];
  dbl H[2][2];
  dbl p[2]; // Newton step
  dbl x_minus_xb[3];
  dbl L;

  dbl tol;
  int niter;

  size_t lhat, l[3];

  /* The state of the vertices of the update base (i.e., `state[i]`
   * goes with index `l[i]`). */
  state_e state[3];

  dbl x[3]; // x[l]
  dbl X[3][3]; // X = [x[l0] x[l1] x[l2]]
  dbl Xt[3][3]; // X'
  dbl XtX[3][3]; // X'*X

  // B-coefs for 9-point triangle interpolation T on base of update
  bb32 T;
};

void utetra_alloc(utetra_s **utetra) {
  *utetra = malloc(sizeof(utetra_s));
}

void utetra_dealloc(utetra_s **utetra) {
  free(*utetra);
  *utetra = NULL;
}

void utetra_init(utetra_s *u, utetra_spec_s const *spec) {
  /* First, validate the spec */

  // TODO: previously I was returning true or false from this function
  // to indicate whether the utetra was OK to do. I'm asserting false
  // now, since I'd like to avoid doing this in this function...

  bool passed_lhat = spec->lhat != (size_t)NO_INDEX;
  bool passed_l0 = spec->l[0] != (size_t)NO_INDEX;
  bool passed_l1 = spec->l[1] != (size_t)NO_INDEX;
  bool passed_l2 = spec->l[2] != (size_t)NO_INDEX;
  bool passed_l = passed_l0 && passed_l1 && passed_l2;

  bool passed_jet0 = jet31t_is_finite(&spec->jet[0]);
  bool passed_jet1 = jet31t_is_finite(&spec->jet[1]);
  bool passed_jet2 = jet31t_is_finite(&spec->jet[2]);
  bool passed_jet = passed_jet0 && passed_jet1 && passed_jet2;
  if (passed_jet0 || passed_jet1 || passed_jet2)
    assert(passed_jet);

  bool passed_xhat = dbl3_isfinite(spec->xhat);
  assert(passed_lhat ^ passed_xhat); // pass exactly one of these

  if (passed_l0 || passed_l1 || passed_l2)
    assert(passed_l);

  bool passed_x0 = dbl3_isfinite(spec->x[0]);
  bool passed_x1 = dbl3_isfinite(spec->x[1]);
  bool passed_x2 = dbl3_isfinite(spec->x[2]);
  bool passed_x = passed_x0 && passed_x1 && passed_x2;
  if (passed_x0 || passed_x1 || passed_x2)
    assert(passed_x);

  assert(passed_l ^ passed_x); // pass exactly one of these

  bool passed_state0 = spec->state[0] != UNKNOWN;
  bool passed_state1 = spec->state[1] != UNKNOWN;
  bool passed_state2 = spec->state[2] != UNKNOWN;
  bool passed_state = passed_state0 && passed_state1 && passed_state2;
  if (passed_state0 || passed_state1 || passed_state2)
    assert(passed_state);

  assert(passed_jet ^ passed_l); // exactly one of these

  /* Initialize f with INFINITY---this needs to be done so that `u`
   * `cmp`s correctly with initialized `utetra` (i.e., if we sort an
   * array of `utetra`, uninitialized `utetra` will move to the back
   * so that they can be easily ignored). */
  u->f = INFINITY;

  /* Initialize Newton iteration variables with dummy values */
  u->lam[0] = u->lam[1] = NAN;
  u->g[0] = u->g[1] = NAN;
  u->H[0][0] = u->H[1][0] = u->H[0][1] = u->H[1][1] = NAN;
  u->p[0] = u->p[1] = NAN;
  u->x_minus_xb[0] = u->x_minus_xb[1] = u->x_minus_xb[2] = NAN;
  u->L = NAN;

  u->tol = spec->tol;

  mesh3_s const *mesh = spec->eik ? eik3_get_mesh(spec->eik) : NULL;

  u->lhat = spec->lhat;

  /* Initialize x(hat) */
  if (passed_lhat) {
    assert(mesh);
    mesh3_copy_vert(mesh, u->lhat, u->x);
  } else { // passed_xhat
    dbl3_copy(spec->xhat, u->x);
  }

  memcpy(u->l, spec->l, sizeof(size_t[3]));

  if (passed_l)
    for (size_t i = 0; i < 3; ++i)
      mesh3_copy_vert(mesh, u->l[i], u->Xt[i]);
  else  // passed_x
    for (size_t i = 0; i < 3; ++i)
      dbl3_copy(spec->x[i], u->Xt[i]);

  dbl33_transposed(u->Xt, u->X);
  dbl33_mul(u->Xt, u->X, u->XtX);

  memcpy(u->state, spec->state, sizeof(state_e[3]));

  /* init jets, depending on how they've been specified */
  jet31t jet[3];
  for (size_t i = 0; i < 3; ++i) {
    jet[i] = passed_jet ? spec->jet[i] : eik3_get_jet(spec->eik, spec->l[i]);
    /* none of them should be singular! */
    assert(jet31t_is_finite(&jet[i]));
  }

  /* Do BB interpolation to set up the coefficients of T. */
  dbl3 T, DT[3];
  for (int i = 0; i < 3; ++i) {
    T[i] = jet[i].f;
    memcpy(DT[i], jet[i].Df, sizeof(dbl3));
  }
  bb32_init_from_3d_data(&u->T, T, DT, u->Xt);
}

/* Check if the point being updated lies in the plane spanned by by
 * x0, x1, and x2. If it does, the update is degenerate. */
bool utetra_is_degenerate(utetra_s const *u) {
  dbl const *x[4] = {u->x, u->Xt[0], u->Xt[1], u->Xt[2]};
  return points_are_coplanar(x);
}

static void set_lambda(utetra_s *cf, dbl const lam[2]) {
  // TODO: question... would it make more sense to use different
  // vectors for a1 and a2? This choice seems to result in a lot of
  // numerical instability. For now I'm fixing this by replacing sums
  // and dot products involving a1 or a2 with the Neumaier equivalent.
  static dbl a1[3] = {-1, 1, 0};
  static dbl a2[3] = {-1, 0, 1};

  dbl b[3], xb[3], tmp1[3], tmp2[3][3], DL[2], D2L[2][2], DT[2], D2T[2][2];
  dbl tmp3[2];

  cf->lam[0] = lam[0];
  cf->lam[1] = lam[1];

  b[1] = lam[0];
  b[2] = lam[1];
  b[0] = 1 - b[1] - b[2];

  assert(dbl3_valid_bary_coord(b));

  dbl33_dbl3_mul(cf->X, b, xb);
  dbl3_sub(cf->x, xb, cf->x_minus_xb);
  cf->L = dbl3_norm(cf->x_minus_xb);
  assert(cf->L > 0);

  dbl33_dbl3_mul(cf->Xt, cf->x_minus_xb, tmp1);
  dbl3_dbl_div(tmp1, -cf->L, tmp1);

  DL[0] = dbl3_ndot(a1, tmp1);
  DL[1] = dbl3_ndot(a2, tmp1);
  assert(dbl2_isfinite(DL));

  dbl3_outer(tmp1, tmp1, tmp2);
  dbl33_sub(cf->XtX, tmp2, tmp2);
  dbl33_dbl_div(tmp2, cf->L, tmp2);

  dbl33_dbl3_nmul(tmp2, a1, tmp1);
  D2L[0][0] = dbl3_ndot(tmp1, a1);
  D2L[1][0] = D2L[0][1] = dbl3_ndot(tmp1, a2);
  dbl33_dbl3_nmul(tmp2, a2, tmp1);
  D2L[1][1] = dbl3_ndot(tmp1, a2);
  assert(dbl22_isfinite(D2L));

  DT[0] = bb32_df(&cf->T, b, a1);
  DT[1] = bb32_df(&cf->T, b, a2);

  D2T[0][0] = bb32_d2f(&cf->T, b, a1, a1);
  D2T[1][0] = D2T[0][1] = bb32_d2f(&cf->T, b, a1, a2);
  D2T[1][1] = bb32_d2f(&cf->T, b, a2, a2);

  cf->f = cf->L + bb32_f(&cf->T, b);
  assert(isfinite(cf->f));

  dbl2_add(DL, DT, cf->g);
  assert(dbl2_isfinite(cf->g));

  dbl22_add(D2L, D2T, cf->H);
  assert(dbl22_isfinite(cf->H));

  /**
   * Finally, compute Newton step solving the minimization problem:
   *
   *     minimize  y’*H*y/2 + [g - H*x]’*y + [x’*H*x/2 - g’*x + f(x)]
   *   subject to  x >= 0
   *               sum(x) <= 1
   *
   * perturbing the Hessian below should ensure a descent
   * direction. (It would be interesting to see if we can remove the
   * perturbation entirely.)
   */

  // Compute the trace and determinant of the Hessian
  dbl tr = cf->H[0][0] + cf->H[1][1];
  dbl det = cf->H[0][0]*cf->H[1][1] - cf->H[0][1]*cf->H[1][0];
  assert(tr != 0 && det != 0);

  // Conditionally perturb the Hessian
  dbl min_eig_doubled = tr - sqrt(tr*tr - 4*det);
  if (min_eig_doubled < 0) {
    cf->H[0][0] -= min_eig_doubled;
    cf->H[1][1] -= min_eig_doubled;
  }

  // Solve a quadratic program to find the next iterate.
  triqp2_s qp;
  dbl22_dbl2_mul(cf->H, lam, tmp3);
  dbl2_sub(cf->g, tmp3, qp.b);
  memcpy((void *)qp.A, (void *)cf->H, sizeof(dbl)*2*2);
  triqp2_solve(&qp);

  // Compute the projected Newton step from the current iterate and
  // next iterate.
  dbl2_sub(qp.x, lam, cf->p);
}

static void step(utetra_s *cf) {
  dbl const atol = 1e-15, c1 = 1e-4;

  dbl lam1[2], f, c1_times_g_dot_p, beta;

  // Get values for current iterate
  dbl lam[2] = {cf->lam[0], cf->lam[1]};
  dbl p[2] = {cf->p[0], cf->p[1]};
  f = cf->f;

  // Do backtracking line search
  beta = 1;
  c1_times_g_dot_p = c1*dbl2_dot(p, cf->g);
  dbl2_saxpy(beta, p, lam, lam1);
  set_lambda(cf, lam1);
  while (cf->f > f + beta*c1_times_g_dot_p + atol) {
    beta *= 0.9;
    dbl2_saxpy(beta, p, lam, lam1);
    set_lambda(cf, lam1);
  }

  ++cf->niter;
}

/**
 * Do a tetrahedron update starting at `lam`, writing the result to
 * `jet`. If `lam` is `NULL`, then the first iterate will be selected
 * automatically.
 */
void utetra_solve(utetra_s *cf, dbl const *lam) {
  cf->niter = 0;

  if (lam == NULL)
    set_lambda(cf, (dbl[2]) {0, 0}); // see, automatically!
  else
    set_lambda(cf, lam);

  int const max_niter = 100;
  dbl const tol = cf->tol*(1 + dbl2_norm(cf->p));
  for (int _ = 0; _ < max_niter; ++_) {
    if (dbl2_norm(cf->p) <= tol)
      break;
    step(cf);
  }
}

static void get_b(utetra_s const *cf, dbl b[3]) {
  assert(!isnan(cf->lam[0]) && !isnan(cf->lam[1]));
  b[0] = 1 - cf->lam[0] - cf->lam[1];
  b[1] = cf->lam[0];
  b[2] = cf->lam[1];
}

dbl utetra_get_value(utetra_s const *cf) {
  return cf->f;
}

void utetra_get_t(utetra_s const *u, dbl t[3]) {
  dbl3_normalized(u->x_minus_xb, t);
}

void utetra_get_jet31t(utetra_s const *cf, jet31t *jet) {
  jet->f = cf->f;

  utetra_get_t(cf, jet->Df);
}

/**
 * Compute the Lagrange multipliers for the constraint optimization
 * problem corresponding to this type of update
 */
static void get_lag_mults(utetra_s const *cf, dbl alpha[3]) {
  dbl const atol = 5e-15;
  dbl b[3] = {1 - dbl2_sum(cf->lam), cf->lam[0], cf->lam[1]};
  alpha[0] = alpha[1] = alpha[2] = 0;
  // TODO: optimize this
  if (fabs(b[0] - 1) < atol) {
    alpha[0] = 0;
    alpha[1] = -cf->g[0];
    alpha[2] = -cf->g[1];
  } else if (fabs(b[1] - 1) < atol) {
    alpha[0] = cf->g[0];
    alpha[1] = 0;
    alpha[2] = cf->g[0] - cf->g[1];
  } else if (fabs(b[2] - 1) < atol) {
    alpha[0] = cf->g[0];
    alpha[1] = cf->g[0] - cf->g[1];
    alpha[2] = 0;
  } else if (fabs(b[0]) < atol) { // b[1] != 0 && b[2] != 0
    alpha[0] = dbl2_sum(cf->g)/2;
    alpha[1] = 0;
    alpha[2] = 0;
  } else if (fabs(b[1]) < atol) { // b[0] != 0 && b[2] != 0
    alpha[0] = 0;
    alpha[1] = -cf->g[0];
    alpha[2] = 0;
  } else if (fabs(b[2]) < atol) { // b[0] != 0 && b[1] != 0
    alpha[0] = 0;
    alpha[1] = 0;
    alpha[2] = -cf->g[1];
  } else {
    assert(dbl3_valid_bary_coord(b));
  }
}

bool utetra_has_interior_point_solution(utetra_s const *cf) {
  dbl alpha[3];
  get_lag_mults(cf, alpha);
  return dbl3_maxnorm(alpha) <= 1e-15;
}

int utetra_cmp(utetra_s const **h1, utetra_s const **h2) {
  utetra_s const *u1 = *h1;
  utetra_s const *u2 = *h2;
  if (u1 == NULL && u2 == NULL) {
    return 0;
  } else if (u2 == NULL) {
    return -1;
  } else if (u1 == NULL) {
    return 1;
  } else {
    dbl f1 = utetra_get_value(u1), f2 = utetra_get_value(u2);
    if (f1 < f2) {
      return -1;
    } else if (f1 > f2) {
      return 1;
    } else {
      return 0;
    }
  }
}

bool utetra_adj_are_optimal(utetra_s const *u1, utetra_s const *u2) {
  dbl const atol = 1e-15;

  return fabs(u1->lam[0] - u2->lam[0]) <= atol
    && fabs(u1->lam[1]) <= atol
    && fabs(u2->lam[1]) <= atol
    && fabs(utetra_get_value(u1) - utetra_get_value(u2)) <= atol;
}

#if JMM_DEBUG
static bool update_inds_are_set(utetra_s const *utetra) {
  return !(
    utetra->l[0] == (size_t)NO_INDEX ||
    utetra->l[1] == (size_t)NO_INDEX ||
    utetra->l[2] == (size_t)NO_INDEX);
}
#endif

#if JMM_DEBUG
static bool all_inds_are_set(utetra_s const *utetra) {
  return utetra->lhat != (size_t)NO_INDEX && update_inds_are_set(utetra);
}
#endif

bool ray_in_tetra_cone(mesh3_s const *mesh, dbl3 p, size_t lc, size_t lv) {
  dbl3 xhat;
  mesh3_copy_vert(mesh, lv, xhat);

  size_t l[3];
  assert(mesh3_cvf(mesh, lc, lv, l));

  dbl3 x[3];
  for (size_t i = 0; i < 3; ++i)
    mesh3_copy_vert(mesh, l[i], x[i]);

  dbl33 dX;
  for (size_t i = 0; i < 3; ++i) {
    dbl3 dx;
    dbl3_sub(x[i], xhat, dx);
    dbl33_set_column(dX, i, dx);
  }

  dbl3 alpha;
  dbl33_dbl3_solve(dX, p, alpha);

  return dbl3_nonneg(alpha);
}

bool ray_in_vertex_cone(mesh3_s const *mesh, dbl3 p, size_t lv) {
  size_t nvc = mesh3_nvc(mesh, lv);
  size_t *vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(mesh, lv, vc);

  bool in_cone = false;
  for (size_t i = 0; i < nvc; ++i)
    if ((in_cone = ray_in_tetra_cone(mesh, p, vc[i], lv)))
      break;

  free(vc);

  return in_cone;
}

bool utetra_update_ray_is_physical(utetra_s const *utetra, eik3_s const *eik) {
#if JMM_DEBUG
  assert(all_inds_are_set(utetra));
#endif

  mesh3_s const *mesh = eik3_get_mesh(eik);

  dbl3 xhat;
  mesh3_copy_vert(mesh, utetra->lhat, xhat);

  dbl3 xlam;
  utetra_get_x(utetra, xlam);

  dbl3 dxlam;
  dbl3_sub(xlam, xhat, dxlam);

  dbl3 dxhat;
  dbl3_sub(xhat, xlam, dxhat);

  /* Check whether the update ray is contained in the cone spanned by
   * the cells incident on `utetra->lhat`. */

  bool ray_end_is_feasible = ray_in_vertex_cone(mesh, dxlam, utetra->lhat);

  /* Check whether the start of the update ray is contained in the set
   * of reasible ray directions at xlam */

  bool ray_start_is_feasible = false;

  size_t l_active[3];
  size_t num_active_constraints = utetra_get_active_inds(utetra, l_active);
  num_active_constraints = 3 - num_active_constraints; // TODO: fix this...

  /* interior minimizer */
  if (num_active_constraints == 0) {
    /* get cells incident on base of update */
    size_t nfc = mesh3_nfc(mesh, l_active);
    size_t *fc = malloc(nfc*sizeof(size_t));
    mesh3_fc(mesh, l_active, fc);

    /* check whether the `dxhat` points into a tetrahedron incident on
     * the base of the update */
    for (size_t i = 0, lv; i < nfc; ++i) {
      mesh3_cfv(mesh, fc[i], l_active, &lv);

      size_t lf[3];
      mesh3_cvf(mesh, fc[i], lv, lf);

      dbl3 x[3];
      for (size_t j = 0; j < 3; ++j)
        mesh3_copy_vert(mesh, lf[j], x[j]);

      dbl3 dx[3];
      dbl3_sub(mesh3_get_vert_ptr(mesh, lv), x[0], dx[0]);
      dbl3_sub(x[1], x[0], dx[1]);
      dbl3_sub(x[2], x[0], dx[2]);

      /* compute the face normal and orient it so that it points
       * outside the tetrahedron */
      dbl3 n;
      dbl3_cross(dx[1], dx[2], n);
      if (dbl3_dot(n, dx[0]) < 0)
        dbl3_negate(n);

      if (dbl3_dot(n, dxhat) > 0) {
        ray_start_is_feasible = true;
        break;
      }
    }

    free(fc);
  }

  /* edge minimizer */
  else if (num_active_constraints == 1) {
    /* gets cells incident on active edge */
    size_t nec = mesh3_nec(mesh, l_active);
    size_t *ec = malloc(nec*sizeof(size_t));
    mesh3_ec(mesh, l_active, ec);

    /* check whether `dxhat` points into a tetrahedron incident on the
     * base of the active edge */
    for (size_t i = 0, le[2]; i < nec; ++i) {
      /* get the opposite edge */
      mesh3_cee(mesh, ec[i], l_active, le);

      dbl3 x[2];
      mesh3_copy_vert(mesh, l_active[0], x[0]);
      mesh3_copy_vert(mesh, l_active[1], x[1]);

      dbl3 y[2];
      mesh3_copy_vert(mesh, le[0], y[0]);
      mesh3_copy_vert(mesh, le[1], y[1]);

      // te = (x1 - x0)/|x1 - x0|
      dbl3 te;
      dbl3_sub(x[1], x[0], te);
      dbl3_normalize(te);

      // q1 = y0 - xlam
      // q1 = q1 - (te@q1)*te
      // q1 /= |q1|
      dbl3 q1;
      dbl3_sub(y[0], xlam, q1);
      dbl te_q1 = dbl3_dot(te, q1);
      for (size_t j = 0; j < 3; ++j)
        q1[j] -= te_q1*te[j];
      dbl3_normalize(q1);

      // q2 = y1 - xlam
      // u2 = q2 - (te@q2)*te;
      // q2 = u2 - (q1@q2)*q1
      // q2 /= |q2|
      dbl3 q2, u2;
      dbl3_sub(y[1], xlam, q2);
      dbl te_q2 = dbl3_dot(te, q2);
      for (size_t j = 0; j < 3; ++j)
        u2[j] = q2[j] - te_q2*te[j];
      dbl q1_q2 = dbl3_dot(q1, q2);
      for (size_t j = 0; j < 3; ++j)
        q2[j] = u2[j] - q1_q2*q1[j];
      dbl3_normalize(q2);

      // thetamax = atan2(q2*u2, q1*u2);
      dbl thetamax = atan2(dbl3_dot(q2, u2), dbl3_dot(q1, u2));

      // theta = atan2(q2*dxhat, q1*dxhat);
      dbl theta = atan2(dbl3_dot(q2, dxhat), dbl3_dot(q1, dxhat));

      if (0 <= theta && theta <= thetamax) {
        ray_start_is_feasible = true;
        break;
      }
    }

    free(ec);
  }

  /* corner minimizer */
  else if (num_active_constraints == 2)
    ray_start_is_feasible = ray_in_vertex_cone(mesh, dxhat, l_active[0]);

  else assert(false);

  return ray_start_is_feasible && ray_end_is_feasible;
}

int utetra_get_num_interior_coefs(utetra_s const *utetra) {
  dbl const atol = 1e-14;
  dbl b[3];
  get_b(utetra, b);
  return (atol < b[0] && b[0] < 1 - atol)
       + (atol < b[1] && b[1] < 1 - atol)
       + (atol < b[2] && b[2] < 1 - atol);
}

size_t utetra_get_l(utetra_s const *utetra) {
  assert(utetra->lhat != (size_t)NO_INDEX);
  return utetra->lhat;
}

void utetra_get_update_inds(utetra_s const *utetra, size_t l[3]) {
  memcpy(l, utetra->l, sizeof(size_t[3]));
}

void utetra_set_update_inds(utetra_s *utetra, size_t const l[3]) {
  memcpy(utetra->l, l, sizeof(size_t[3]));
}

/* Get the triangle of `u`. If `u` is split, this does *not* return
 * the triangle of the split updates (since that would not be
 * well-defined!). */
static tri3 get_tri(utetra_s const *u) {
  tri3 tri;
  memcpy(tri.v, u->Xt, sizeof(dbl[3][3]));
  return tri;
}

/**
 * Check whether the optimum of `u` is incident on the base of
 * `u_other`.
 */
bool utetra_opt_inc_on_other_utetra(utetra_s const *u, utetra_s const *u_other) {
  dbl const atol = 1e-14;
  dbl xb[3];
  utetra_get_x(u, xb);
  tri3 tri = get_tri(u_other);
  return tri3_dist(&tri, xb) < atol;
}

void utetra_get_x(utetra_s const *u, dbl x[3]) {
  dbl b[3];
  get_b(u, b);
  dbl33_dbl3_mul(u->X, b, x);
}

size_t utetra_get_active_inds(utetra_s const *utetra, size_t l[3]) {
#if JMM_DEBUG
  assert(update_inds_are_set(utetra));
#endif

  dbl const atol = 1e-14;

  dbl alpha[3];
  get_lag_mults(utetra, alpha);

  size_t num_active = 0;
  l[0] = l[1] = l[2] = (size_t)NO_INDEX;
  for (int i = 0; i < 3; ++i)
    if (fabs(alpha[i]) < atol)
      l[num_active++] = utetra->l[i];
  return num_active;
}

par3_s utetra_get_parent(utetra_s const *utetra) {
  par3_s par; get_b(utetra, par.b);
  assert(dbl3_valid_bary_coord(par.b));
  memcpy(par.l, utetra->l, sizeof(size_t[3]));
  return par;
}

bool utetra_approx_hess(utetra_s const *u, dbl h, dbl33 hess) {
  dbl const atol = 1e-14;

  if (h < atol)
    return false;

  utetra_s u_ = *u;

  dbl dx[3], t[3], lam[2] = {u->lam[0], u->lam[1]};

  dbl33_zero(hess);

  /* Approximate the Hessian using central differences. For each i,
   * compute (t(x + h*e_i) - t(x - h*e_i))/(2*h) and store in the ith
   * row of hess */
  for (size_t i = 0; i < 3; ++i) {
    dbl3_zero(dx);

    /* Compute t(x + h*e_i) */
    dx[i] = h;
    dbl3_add(u->x, dx, u_.x);
    utetra_solve(&u_, lam);
    utetra_get_t(&u_, t);
    dbl3_add_inplace(hess[i], t);

    /* Compute t(x - h*e_i) */
    dx[i] = -h;
    dbl3_add(u->x, dx, u_.x);
    utetra_solve(&u_, lam);
    utetra_get_t(&u_, t);
    dbl3_sub_inplace(hess[i], t);

    /* Set H[i, :] = (t(x + h*e_i) - t(x - h*e_i))/(2*h) */
    dbl3_dbl_div_inplace(hess[i], 2*h);
  }

  /* Make sure the Hessian is symmetric */
  dbl33_symmetrize(hess);

  return true;
}

void utetra_get_other_inds(utetra_s const *utetra, size_t li, size_t l[2]) {
  assert(utetra->l[0] == li || utetra->l[1] == li || utetra->l[2]);

  for (size_t i = 0, j = 0; i < 3; ++i)
    if (utetra->l[i] != li)
      l[j++] = utetra->l[i];
}

bool utetra_index_is_active(utetra_s const *utetra, size_t l) {
  assert(l == utetra->l[0] || l == utetra->l[1] || l == utetra->l[2]);

  dbl lam0 = 1 - utetra->lam[0] - utetra->lam[1];

  if (l == utetra->l[0])
    return lam0 != 0;
  else if (l == utetra->l[1])
    return utetra->lam[0] != 0;
  else if (l == utetra->l[2])
    return utetra->lam[1] != 0;
  else
    assert(false);
}

bool utetras_have_same_minimizer(utetra_s const *u1, utetra_s const *u2) {
  size_t l1[3], l2[3];
  (void)utetra_get_active_inds(u1, l1);
  (void)utetra_get_active_inds(u2, l2);

  SORT3(l1[0], l1[1], l1[2]);
  SORT3(l2[0], l2[1], l2[2]);
  if (l1[0] != l2[0] || l1[1] != l2[1] || l1[2] != l2[2])
    return false;

  dbl x1[3], x2[3];
  utetra_get_x(u1, x1);
  utetra_get_x(u2, x2);
  return dbl3_dist(x1, x2) <= 1e-13;
}

bool utetras_have_same_inds(utetra_s const *u1, utetra_s const *u2) {
  if (u1->lhat != u2->lhat)
    return false;

  size_t l1[3];
  memcpy(l1, u1->l, sizeof(size_t[3]));
  SORT3(l1[0], l1[1], l1[2]);

  size_t l2[3];
  memcpy(l2, u2->l, sizeof(size_t[3]));
  SORT3(l2[0], l2[1], l1[2]);

  return l1[0] == l2[0] && l1[1] == l2[1] && l1[2] == l2[2];
}

#if JMM_TEST
void utetra_step(utetra_s *u) {
  step(u);
}

void utetra_get_lambda(utetra_s const *u, dbl lam[2]) {
  lam[0] = u->lam[0];
  lam[1] = u->lam[1];
}

void utetra_set_lambda(utetra_s *u, dbl const lam[2]) {
  set_lambda(u, lam);
}

size_t utetra_get_num_iter(utetra_s const *u) {
  return u->niter;
}
#endif
