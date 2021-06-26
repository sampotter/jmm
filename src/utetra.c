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
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
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
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
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
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
    },
    .tol = pow(h, 4)
  };
}

utetra_spec_s utetra_spec_from_ptrs(mesh3_s const *mesh, jet3 const *jet,
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

bool utetra_init(utetra_s *u, utetra_spec_s const *spec) {
  /* First, validate the spec */

  bool passed_lhat = spec->lhat != (size_t)NO_INDEX;
  bool passed_l0 = spec->l[0] != (size_t)NO_INDEX;
  bool passed_l1 = spec->l[1] != (size_t)NO_INDEX;
  bool passed_l2 = spec->l[2] != (size_t)NO_INDEX;
  bool passed_l = passed_l0 && passed_l1 && passed_l2;

  bool passed_jet0 = jet3_is_finite(&spec->jet[0]);
  bool passed_jet1 = jet3_is_finite(&spec->jet[1]);
  bool passed_jet2 = jet3_is_finite(&spec->jet[2]);
  bool passed_jet = passed_jet0 && passed_jet1 && passed_jet2;
  if (passed_jet0 || passed_jet1 || passed_jet2)
    assert(passed_jet);

#if JMM_DEBUG
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
#endif

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

  jet3 jet[3];
  for (size_t i = 0; i < 3; ++i)
    jet[i] = passed_jet ? spec->jet[i] : eik3_get_jet(spec->eik, spec->l[i]);

  /* Figure out which jets lack gradient information */
  bool pt_src[3];
  size_t num_pt_srcs = 0;
  for (size_t i = 0; i < 3; ++i)
    num_pt_srcs += pt_src[i] = jet3_is_point_source(&jet[i]);

  if (num_pt_srcs == 3) {
    /* If all of the jets lack point source data, just return false
     * now. */
    return false;
  } else if (num_pt_srcs > 0) {
    /* If some but not all of the jets are point sources, linearly
     * interpolate the Bezier coefficients from the available gradient
     * data. This is inaccurate but will unstick the solver in a few
     * places, especially near the boundary of the physical rays
     * emitted by diffracting edges. */
    bb32_init_from_jets(&u->T, jet, u->Xt);
  } else {
    /* If we have all the gradient data we need, do regular ol' BB
     * interpolation. */
    dbl T[3], DT[3][3];
    for (int i = 0; i < 3; ++i) {
      T[i] = jet[i].f;
      memcpy(DT[i], &jet[i].fx, sizeof(dbl[3]));
    }
    bb32_init_from_3d_data(&u->T, T, &DT[0], u->Xt);
  }

  /* If we're updating from a virtual boundary face, then we want to
   * unconditionally accept updates that we do from the wrong
   * side. This is important when we're solving the extended problem,
   * since this will allow us to propagate the eikonal from either
   * side of a face with reflection BCs. */
  if (mesh3_bdf_is_virtual(mesh, u->l))
    return true;

  // Compute the surface normal for the plane spanned by (x1 - x0, x2
  // - x0), using DT[i] to determine its orientation. Return whether x
  // is on the right side of this plane.

  dbl n[3];
  dbl dx[2][3];
  dbl3_sub(u->Xt[1], u->Xt[0], dx[0]);
  dbl3_sub(u->Xt[2], u->Xt[0], dx[1]);
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);

  int sgn[3];
  for (size_t i = 0; i < 3; ++i) {
    sgn[i] = pt_src[i] ? 0 : signum(dbl3_dot(&jet[i].fx, n));
  };

  /* Verify that the jets don't span the same plane as the base of the
   * update does... There's no way of handling this when c == 1, but
   * if c != 1, then we can try. */
  if (sgn[0] == 0 && sgn[1] == 0 && sgn[2] == 0)
    return false;

  /* If the jets point on either side of the update, this is a bit of
   * a weird an unphysical configuration, so let's just return false
   * now... */
  int sgnmax = MAX(sgn[0], MAX(sgn[1], sgn[2]));
  int sgnmin = MIN(sgn[0], MIN(sgn[1], sgn[2]));
  if (sgnmax == 1 && sgnmin == -1)
    return false;

  /* We might have some jets that lie in the update base and some that
   * don't. Use the ones that don't to determine which way to orient
   * the update base. */
  assert(!(sgnmax == 1 && sgnmin == -1));
  if (sgnmin == -1)
    dbl3_negate(n);

  dbl x0_minus_x[3];
  dbl3_sub(u->Xt[0], u->x, x0_minus_x);

  dbl dot = -dbl3_dot(x0_minus_x, n);

  return dot > 0;
}

void utetra_deinit(utetra_s *u) {
  (void)u;
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

  // I think we're actually supposed to use cf->g here instead of
  // cf->p, but maybe cf->p is more appropriate for constraint
  // Newton's method, since it takes into account the constraints
  // while cf->p doesn't.
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

void utetra_get_jet(utetra_s const *cf, jet3 *jet) {
  jet->f = cf->f;
  jet->fx = cf->x_minus_xb[0]/cf->L;
  jet->fy = cf->x_minus_xb[1]/cf->L;
  jet->fz = cf->x_minus_xb[2]/cf->L;
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

/* Check if `u` emits a "terminal ray". This is basically a ray which
 * is emitted from a node with BCs which doesn't "match the
 * singularity structure" of the problem type; equivalently, the ray
 * is emitted from the boundary of the *subset of the boundary* which
 * has BCs supplied.
 *
 * For an edge diffraction problem, this means the ray was emitted
 * from an endpoint of a diffracting edge.
 *
 * For a reflection, this means that the ray was emitted from the edge
 * of reflecting part of the boundary (that is, from the
 * silhouette).
 *
 * Note: a point source can't emit terminal rays. */
bool utetra_emits_terminal_ray(utetra_s const *u, eik3_s const *eik) {
  par3_s par = utetra_get_parent(u);

  if (!par3_is_on_BC_boundary(&par, eik))
    return false;

  ftype_e ftype = eik3_get_ftype(eik);
  assert(ftype != FTYPE_POINT_SOURCE);

  /* Pretty sure this is correct... Double check */
  if (ftype == FTYPE_EDGE_DIFFRACTION)
    return true;

  assert(ftype == FTYPE_REFLECTION);

  size_t num_active = par3_num_active(&par);
  assert(num_active < 3);

  /* If we've reached this point, we know a few things:
   *
   * 1) The Lagrange multipliers aren't too small
   * 2) The minimizer is on the "BC boundary"
   * 3) We're computing a reflection
   * 4) This isn't an interior point minimizer
   *
   * So, if there are two active vertices, then we know that we can
   * accept the solution, since it will by an interior point minimizer
   * if the domain is restricted to the active edge on the "BC
   * boundary"
   *
   * On the other hand, if there's only one active vertex, then we
   * don't know whether we should accept this update or not. We could
   * try to check whether the active Lagrange multipliers "point" to
   * another part of the "BC boundary" which is collinear with the one
   * we're working with now, but a simpler thing to do is just reject
   * the update now. This will cause it to go into the "old update"
   * list. Later, if this truly *is* a good minimizer, we'll find
   * another `utetra` that matches it---we will then use the pair of
   * them to deduce that we should accept the minimizer. */
  return num_active == 2;
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

static ray3 get_ray(utetra_s const *utetra) {
  ray3 ray;
  dbl b[3];
  get_b(utetra, b);
  dbl33_dbl3_mul(utetra->X, b, ray.org);
  dbl3_normalized(utetra->x_minus_xb, ray.dir);
  return ray;
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

static dbl get_L(utetra_s const *u) {
  return u->L;
}

bool utetra_update_ray_is_physical(utetra_s const *utetra, eik3_s const *eik) {
  assert(all_inds_are_set(utetra));

  if (utetra_updated_from_refl_BCs(utetra, eik))
    return true;

  size_t const *l = utetra->l;

  mesh3_s const *mesh = eik3_get_mesh(eik);

  // TODO: the following section where we check to see if the stuff
  // below gives "an interior ray" can be wrapped up and reused for
  // both this and the corresponding section in utri.c...

  ray3 ray = get_ray(utetra);

  /**
   * First, check if the start of the ray is "in free space". To do
   * this, we just check if we can a cell containing a point just
   * ahead of and just behind the origin of the ray.
   */

  // Get points just before and just after the start of the ray. We
  // perturb forward and backward by one half of the minimum triangle
  // altitude (taken of the entire mesh). This is to ensure that the
  // perturbations are small enough to stay inside neighboring
  // tetrahedra, but also large enough to take into consideration the
  // characteristic length scale of the mesh.
  dbl xm[3], xp[3], h = mesh3_get_min_edge_length(mesh)/4; // TODO: rename t...
  ray3_get_point(&ray, -h, xm);
  ray3_get_point(&ray, h, xp);
  assert(dbl3_isfinite(xm));
  assert(dbl3_isfinite(xp));

  array_s *cells;
  array_alloc(&cells);
  array_init(cells, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Fill `cells` with all of the cells incident on any of the
   * neighboring update indices. */
  for (int i = 0; i < 3; ++i) {
    int nvc = mesh3_nvc(mesh, l[i]);
    size_t *vc = malloc(nvc*sizeof(size_t));
    mesh3_vc(mesh, l[i], vc);
    for (int j = 0; j < nvc; ++j)
      if (!array_contains(cells, &vc[j]))
        array_append(cells, &vc[j]);
    free(vc);
  }

  bool xm_in_cell = false, xp_in_cell = false;
  for (size_t i = 0, lc; i < array_size(cells); ++i) {
    array_get(cells, i, &lc);
    xm_in_cell |= mesh3_cell_contains_point(mesh, lc, xm);
    xp_in_cell |= mesh3_cell_contains_point(mesh, lc, xp);
    if (xm_in_cell && xp_in_cell)
      break;
  }
  if (!xm_in_cell || !xp_in_cell) {
    array_deinit(cells);
    array_dealloc(&cells);
    return false;
  }

  /* Next, we'll pull out the boundary faces incident on these cells
   * and check if the ray intersects any of them. If it does, the ray
   * is unphysical. */

  array_s *bdf;
  array_alloc(&bdf);
  array_init(bdf, sizeof(size_t[3]), ARRAY_DEFAULT_CAPACITY);

  // Get the incident boundary faces
  size_t cf[4][3];
  for (size_t i = 0, lc; i < array_size(cells); ++i) {
    array_get(cells, i, &lc);
    mesh3_cf(mesh, lc, cf);
    for (size_t j = 0; j < 4; ++j) {
      if (array_contains(bdf, cf[j]))
        continue;
      /* skip faces containing the update point to make the
       * intersection test simpler */
      if (point_in_face(utetra->lhat, cf[j]))
        continue;

      /* Check if `cf[j]` indexes a real boundary face. We don't want
       * to accidentally intersect a virtual boundary face here. */
      if (mesh3_is_bdf(mesh, cf[j], false /* virtual not OK */))
        array_append(bdf, cf[j]);
    }
  }

  dbl L = get_L(utetra);

  // Check for intersections with the nearby boundary faces
  bool found_bdf_isect = false;
  size_t lf[3];
  for (size_t i = 0; i < array_size(bdf); ++i) {
    array_get(bdf, i, lf);
    tri3 tri = mesh3_get_tri(mesh, lf);
    if (ray3_and_tri3_are_parallel(&ray, &tri))
      continue;
    dbl t;
    bool isect = ray3_intersects_tri3(&ray, &tri, &t);
    if (isect && t <= L) {
      found_bdf_isect = true;
      break;
    }
  }

  array_deinit(bdf);
  array_dealloc(&bdf);

  array_deinit(cells);
  array_dealloc(&cells);

  if (found_bdf_isect)
    return false;

  /* Next, check and see if the point just before the end of the ray
   * lies in a cell. */

  dbl xhatm[3];
  ray3_get_point(&ray, L - h, xhatm);

  int nvc = mesh3_nvc(mesh, utetra->lhat);
  size_t *vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(mesh, utetra->lhat, vc);

  bool xhatm_in_cell = false;
  for (int i = 0; i < nvc; ++i) {
    xhatm_in_cell = mesh3_cell_contains_point(mesh, vc[i], xhatm);
    if (xhatm_in_cell)
      break;
  }

  free(vc);

  return xhatm_in_cell;
}

bool utetra_updated_from_refl_BCs(utetra_s const *utetra, eik3_s const *eik) {
  ftype_e ftype = eik3_get_ftype(eik);
  if (ftype != FTYPE_REFLECTION)
    return false;

  par3_s par = utetra_get_parent(utetra);
  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active]; // TODO: unused
  par3_get_active(&par, l, b);

  for (size_t i = 0; i < num_active; ++i)
    if (!eik3_has_BCs(eik, l[i]))
      return false;

  return true;
}

int utetra_get_num_interior_coefs(utetra_s const *utetra) {
  dbl const atol = 1e-14;
  dbl b[3];
  get_b(utetra, b);
  return (b[0] > atol) + (b[1] > atol) + (b[2] > atol);
}

static size_t get_shared_inds(utetra_s const *u1, utetra_s const *u2, size_t *l) {
  assert(update_inds_are_set(u1));
  assert(update_inds_are_set(u2));

  size_t const *l1 = u1->l, *l2 = u2->l;

  size_t i = 0;

  for (size_t j = 0; j < 3; ++j)
    if (l1[j] == l2[0] || l1[j] == l2[1] || l1[j] == l2[2])
      l[i++] = l1[j];

  return i; // == num_shared_inds
}

static bool get_point_for_index(utetra_s const *utetra, size_t l, dbl x[3]) {
  assert(update_inds_are_set(utetra));

  for (int i = 0; i < 3; ++i)
    if (utetra->l[i] == l) {
      memcpy(x, utetra->Xt[i], sizeof(dbl[3]));
      return true;
    }

  return false;
}

static bool get_op_ind(utetra_s const *utetra, size_t const le[2], size_t *l) {
  bool found[3] = {false, false, false};

  for (int i = 0; i < 3; ++i)
    if (utetra->l[i] == le[0] || utetra->l[i] == le[1])
      found[i] = true;

  if (found[0] + found[1] + found[2] == 2) {
    for (int i = 0; i < 3; ++i)
      if (!found[i]) {
        *l = utetra->l[i];
        return true;
      }
    die();
  }

  return false;
}

static size_t get_num_equal(utetra_s const **u, size_t n) {
  if (n < 2)
    return n;

  dbl x[3], y[3];

  /* The first thing to check is whether the jet computed by each
   * update and the start of the update ray parametrized by each
   * update is the same. This is necessary but not sufficient. */

  jet3 jet[2];

  // Prefetch the first coords and jet
  utetra_get_x(u[0], x);
  utetra_get_jet(u[0], &jet[0]);

  size_t neq = 1;

  dbl const atol = u[0]->tol;

  /* We can't handle varying tolerances at the moment. */
  for (size_t i = 1; i < n; ++i)
    assert(atol == u[i]->tol);

  // First, check that the update data for each utetra is the
  // same. This is cheaper than the topological check that follows.
  for (neq = 1; neq < n; ++neq) {
    // Get the next jet and check that it's finite
    utetra_get_jet(u[neq], &jet[1]);
    if (!jet3_is_finite(&jet[1]))
      break;

    // Get the next coords
    utetra_get_x(u[neq], y);

    // Check if the base of the update rays coincide...
    if (dbl3_dist(x, y) > atol)
      break;

    // Check if the computed jets are the same...
    if (!jet3_approx_eq(&jet[0], &jet[1], atol))
      break;

    // Swap the coords and jet that we just fetched to make way for
    // the next ones
    for (int j = 0; j < 3; ++j)
      SWAP(x[j], y[j]);
    SWAP(jet[0], jet[1]);
  }

  return neq;
}

static size_t count_unique_inds(utetra_s const **u, size_t n) {
  array_s *l;

  array_alloc(&l);
  array_init(l, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (!array_contains(l, &u[i]->l[j]))
        array_append(l, &u[i]->l[j]);

  size_t num_unique = array_size(l);

  array_deinit(l);
  array_dealloc(&l);

  return num_unique;
}

static void get_shared_and_op_inds_3(utetra_s const **u, size_t *l_shared,
                                     size_t l_op[3]) {
  assert(count_unique_inds(u, 3) == 4);

  /* Get the indices shared by two pairs of the `utetra`s. */
  size_t shared[3][2];
  get_shared_inds(u[0], u[1], shared[0]);
  get_shared_inds(u[0], u[2], shared[1]);
  get_shared_inds(u[1], u[2], shared[2]);

  /* Find the index common to `shared[0]` and `shared[1]`. This is the
   * index shared by all three `utetra`... */
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      if (shared[0][i] == shared[1][j]) {
        *l_shared = shared[0][i];
        break;
      }

  /* (double-check that `*l_shared` is in `shared[2]`) */
  assert(*l_shared == shared[2][0] || *l_shared == shared[2][1]);

  /* ... the indices *not* common to `shared1` and `shared2` should
   * therefore go to `l_op`. */
  for (size_t i = 0; i < 3; ++i)
    l_op[i] = shared[i][0] == *l_shared ? shared[i][1] : shared[i][0];

  /* (make sure that each entry of `l_op` is distinct!) */
  assert(l_op[0] != l_op[1] && l_op[1] != l_op[2]);
}

/* Checks whether `n` tetrahedron updates stored in `utetra[i]` yield
 * the same update. This is a pretty complicated check... Need to
 * explain the cases this handles a bit better. */
bool utetras_yield_same_update(utetra_s const **u, size_t n) {
  assert(n > 1);

  size_t neq = get_num_equal(u, n);

  if (neq == 1)
    return false;

  /* We want to check whether the mesh comprised of the `n` update
   * bases (triangles) are "thick" from the perspective of the update
   * ray. Basically, what we're trying to check is whether, after
   * projecting the union of the update bases into the plane normal to
   * the update ray, the ray hits an interior point or not. E.g., we
   * want to rule out rays that only *graze* the update triangles, or
   * weird sets of updates that aren't manifold at the intersection
   * point. As a reminder, this function should only be called when we
   * were unable to find a single update with an interior point
   * solution. */

  jet3 jet;
  utetra_get_jet(u[0], &jet);

  if (neq == 2) {
    // Get the indices of the shared edge, and validate that this edge
    // is incident on the base of each tetrahedron update
    size_t l_shared[2];
    size_t num_shared = get_shared_inds(u[0], u[1], l_shared);
    assert(num_shared < 3);

    if (num_shared == 1)
      return false;

    dbl x[3], dx[3], y[3], dy[3], normal[3];

    // Get the vectors defining the line spanned by the shared edge
    get_point_for_index(u[0], l_shared[0], y);
    get_point_for_index(u[0], l_shared[1], dy);
    dbl3_sub_inplace(dy, y);

    // Get the normal for the plane spanned by the edge and update ray,
    // taking advantage of the fact that all the update rays can be
    // assumed to be equal at this point
    dbl3_cross(dy, &jet.fx, normal);
    dbl3_normalize(normal); // (probably overkill)

    // For each tetrahedron update, pull out the points corresponding to
    // index that isn't incident on the shared edge, and translate them
    // so that the edge passes through the origin. Then, compute the dot
    // product between this point and the normal vector `n`, keeping
    // track of the minimum and maximum dot products.
    dbl amin = INFINITY, amax = -INFINITY;
    for (size_t i = 0, l_op; i < neq; ++i) {
      get_op_ind(u[i], l_shared, &l_op);
      get_point_for_index(u[i], l_op, x);
      dbl3_sub(x, y, dx);
      dbl a = dbl3_dot(normal, dx);
      amin = fmin(amin, a);
      amax = fmax(amax, a);
    }

    // We want to check that at least two of the "opposite update
    // indices" lie on either side of the plane spanned by the edge and
    // the update ray. The dot products `amin` and `amax` are the
    // extreme dot products values, so we just check this here.
    return amin < 0 && 0 < amax;
  }

  if (neq == 3) {
    size_t num_unique = count_unique_inds(u, neq);
    assert(num_unique >= 4);
    if (num_unique > 4)
      return false;

    size_t l_shared, l_op[3];
    get_shared_and_op_inds_3(u, &l_shared, l_op);

    /* Get the `tri3` indexed by `l_op`. This is a little circuitous
     * since we don't know the provenance of each index in `l_op` at
     * this point, so we need to do a quick search...  */
    tri3 tri;
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        if (get_point_for_index(u[j], l_op[i], tri.v[i]))
          continue;

    /* Get the update ray with which we'll try to intersect `tri`. */
    ray3 ray = get_ray(u[0]);

    dbl unused;

    /* First check the default orientation of the update ray... */
    if (ray3_intersects_tri3(&ray, &tri, &unused))
      return true;

    /* ... and if that doesn't intersect `tri`, flip it around and try
     * the other orientation, as well. We do this because what we
     * *really* want to do is intersect the line spanned by `ray` with
     * `tri`. */
    // TODO: implement a `line_intersects_tri3` function to simplify
    // this code...
    dbl3_negate(ray.dir);
    return ray3_intersects_tri3(&ray, &tri, &unused);
  }

#if JMM_DEBUG
  die();
#else
  return false;
#endif
}

size_t utetra_get_l(utetra_s const *utetra) {
  assert(utetra->lhat != (size_t)NO_INDEX);
  return utetra->lhat;
}

void utetra_set_update_inds(utetra_s *utetra, size_t l[3]) {
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
  assert(update_inds_are_set(utetra));

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

void utetra_get_t(utetra_s const *u, dbl t[3]) {
  dbl3_normalized(u->x_minus_xb, t);
}

dbl utetra_get_L(utetra_s const *u) {
  return u->L;
}

bool utetra_approx_hess(utetra_s const *u, dbl h, dbl33 hess) {
  dbl const atol = 1e-14;

  if (h < atol)
    return false;

  utetra_s u_ = *u;

  dbl dx[3], t[3], lam[2];
  utetra_get_lambda(u, lam);

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
