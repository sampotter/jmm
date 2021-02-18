#include "utetra.h"

#include <assert.h>
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

struct utetra {
  dbl lam[2]; // Current iterate

  dbl f;
  dbl g[2];
  dbl H[2][2];

  dbl p[2]; // Newton step

  dbl angles[3];

  dbl x[3]; // x[l]
  dbl X[3][3]; // X = [x[l0] x[l1] x[l2]]
  dbl Xt[3][3]; // X'
  dbl XtX[3][3]; // X'*X

  dbl L, T[3];

  // B-coefs for 9-point triangle interpolation T on base of update
  bb32 bb_T;

  dbl x_minus_xb[3];

  size_t lhat, l[3];

  int niter;
};

void utetra_alloc(utetra_s **utetra) {
  *utetra = malloc(sizeof(utetra_s));
}

void utetra_dealloc(utetra_s **utetra) {
  free(*utetra);
  *utetra = NULL;
}

void utetra_reset(utetra_s *cf) {
  cf->niter = 0;
  cf->f = INFINITY;
}

bool utetra_init_from_eik3(utetra_s *cf, eik3_s const *eik,
                           size_t l, size_t l0, size_t l1, size_t l2) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  jet3 const *jet = eik3_get_jet_ptr(eik);
  return utetra_init_from_ptrs(cf, mesh, jet, l, l0, l1, l2);
}

bool utetra_init_from_ptrs(utetra_s *cf, mesh3_s const *mesh, jet3 const *jet,
                           size_t l, size_t l0, size_t l1, size_t l2) {
  dbl x[3];
  mesh3_copy_vert(mesh, l, x);

  dbl Xt[3][3];
  mesh3_copy_vert(mesh, l0, Xt[0]);
  mesh3_copy_vert(mesh, l1, Xt[1]);
  mesh3_copy_vert(mesh, l2, Xt[2]);

  jet3 jet_[3] = {jet[l0], jet[l1], jet[l2]};

  assert(jet3_is_finite(&jet_[0]));
  assert(jet3_is_finite(&jet_[1]));
  assert(jet3_is_finite(&jet_[2]));

  cf->lhat = l;

  cf->l[0] = l0;
  cf->l[1] = l1;
  cf->l[2] = l2;

  return utetra_init(cf, x, Xt, jet_);
}

bool utetra_init(utetra_s *cf, dbl const x[3], dbl const Xt[3][3],
                 jet3 const jet[3]) {
  // Initialize lambda with a dummy value so that we can check for bad
  // accesses later
  cf->lam[0] = cf->lam[1] = NAN;

  memcpy(cf->x, x, 3*sizeof(dbl));
  memcpy(cf->Xt, Xt, 3*3*sizeof(dbl));

  dbl33_transposed(cf->Xt, cf->X);
  dbl33_mul(cf->Xt, cf->X, cf->XtX);

  dbl Xt_minus_x[3][3];
  for (int i = 0; i < 3; ++i) {
    dbl3_sub(cf->Xt[i], cf->x, Xt_minus_x[i]);
    dbl3_normalize(Xt_minus_x[i]);
  }
  for (int i = 0; i < 3; ++i) {
    cf->angles[i] = dbl3_dot(Xt_minus_x[i], Xt_minus_x[(i + 1) % 3]);
  }

  dbl DT[3][3];
  for (int i = 0; i < 3; ++i) {
    cf->T[i] = jet[i].f;
    DT[i][0] = jet[i].fx;
    DT[i][1] = jet[i].fy;
    DT[i][2] = jet[i].fz;
  }
  bb32_init_from_3d_data(&cf->bb_T, cf->T, &DT[0], cf->Xt);

  cf->niter = 0;

  // Compute the surface normal for the plane spanned by (x1 - x0, x2
  // - x0), using DT[i] to determine its orientation. Return whether x
  // is on the right side of this plane.
  dbl n[3];
  dbl dx[2][3];
  dbl3_sub(cf->Xt[1], cf->Xt[0], dx[0]);
  dbl3_sub(cf->Xt[2], cf->Xt[0], dx[1]);
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
  int sgn[3] = {
    signum(dbl3_dot(DT[0], n)),
    signum(dbl3_dot(DT[1], n)),
    signum(dbl3_dot(DT[2], n))
  };
  if (sgn[0] == -1)
    dbl3_negate(n);
  dbl dot = -dbl3_dot(Xt_minus_x[0], n) > 0;
  return dot > 0;
}

bool utetra_is_degenerate(utetra_s const *cf) {
  // Check if the point being updated lies in the plane spanned by by
  // x0, x1, and x2. If it does, the update is degenerate.
  dbl const *x[4] = {cf->x, cf->Xt[0], cf->Xt[1], cf->Xt[2]};
  return points_are_coplanar(x);
}

bool utetra_is_causal(utetra_s const *cf) {
  return cf->angles[0] >= 0 && cf->angles[1] >= 0 && cf->angles[2] >= 0;
}

/**
 * Do a tetrahedron update starting at `lam`, writing the result to
 * `jet`. This assumes that `utetra_set_lambda` has already been
 * called, so that `cf` is currently at `lam`.
 */
void utetra_solve(utetra_s *cf) {
  //dbl const rtol = 1e-15;
  //dbl const atol = 5e-15;

  dbl const c1 = 1e-2;
  dbl lam[2], lam1[2], p[2], f, c1_times_g_dot_p, beta;

  for (cf->niter = 0; cf->niter < 20; ++cf->niter) {
    // Get values for current iterate
    lam[0] = cf->lam[0]; lam[1] = cf->lam[1];
    p[0] = cf->p[0]; p[1] = cf->p[1];
    f = cf->f;

    // Do backtracking line search
    beta = 1;
    c1_times_g_dot_p = c1*dbl2_dot(p, cf->g);
    dbl2_saxpy(beta, p, lam, lam1);
    utetra_set_lambda(cf, lam1);
    while (cf->f > f + beta*c1_times_g_dot_p) {
      beta /= 2;
      dbl2_saxpy(beta, p, lam, lam1);
      utetra_set_lambda(cf, lam1);
    }
  };
}

void utetra_get_lambda(utetra_s const *cf, dbl lam[2]) {
  assert(!isnan(cf->lam[0]) && !isnan(cf->lam[1]));
  lam[0] = cf->lam[0];
  lam[1] = cf->lam[1];
}

void utetra_set_lambda(utetra_s *cf, dbl const lam[2]) {
  static dbl a1[3] = {-1, 1, 0};
  static dbl a2[3] = {-1, 0, 1};

  static dbl const atol = 1e-15;

  dbl b[3], xb[3], tmp1[3], tmp2[3][3], DL[2], D2L[2][2], DT[2], D2T[2][2];
  dbl tmp3[2];

  cf->lam[0] = lam[0];
  cf->lam[1] = lam[1];

  b[1] = lam[0];
  b[2] = lam[1];
  b[0] = 1 - b[1] - b[2];

  assert(b[0] >= -atol);
  assert(b[1] >= -atol);
  assert(b[2] >= -atol);

  dbl33_dbl3_mul(cf->X, b, xb);
  dbl3_sub(cf->x, xb, cf->x_minus_xb);
  cf->L = dbl3_norm(cf->x_minus_xb);
  assert(cf->L > 0);

  dbl33_dbl3_mul(cf->Xt, cf->x_minus_xb, tmp1);
  dbl3_dbl_div(tmp1, -cf->L, tmp1);

  DL[0] = dbl3_dot(a1, tmp1);
  DL[1] = dbl3_dot(a2, tmp1);
  assert(dbl2_isfinite(DL));

  dbl3_outer(tmp1, tmp1, tmp2);
  dbl33_sub(cf->XtX, tmp2, tmp2);
  dbl33_dbl_div(tmp2, cf->L, tmp2);

  dbl33_dbl3_mul(tmp2, a1, tmp1);
  D2L[0][0] = dbl3_dot(tmp1, a1);
  D2L[1][0] = D2L[0][1] = dbl3_dot(tmp1, a2);
  dbl33_dbl3_mul(tmp2, a2, tmp1);
  D2L[1][1] = dbl3_dot(tmp1, a2);
  assert(dbl22_isfinite(D2L));

  DT[0] = bb32_df(&cf->bb_T, b, a1);
  DT[1] = bb32_df(&cf->bb_T, b, a2);

  D2T[0][0] = bb32_d2f(&cf->bb_T, b, a1, a1);
  D2T[1][0] = D2T[0][1] = bb32_d2f(&cf->bb_T, b, a1, a2);
  D2T[1][1] = bb32_d2f(&cf->bb_T, b, a2, a2);

  cf->f = cf->L + bb32_f(&cf->bb_T, b);
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

  // Conditionally perturb the Hessian
  dbl tr = cf->H[0][0] + cf->H[1][1];
  dbl det = cf->H[0][0]*cf->H[1][1] - cf->H[0][1]*cf->H[1][0];
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

void utetra_get_bary_coords(utetra_s const *cf, dbl b[3]) {
  assert(!isnan(cf->lam[0]) && !isnan(cf->lam[1]));
  b[0] = 1 - cf->lam[0] - cf->lam[1];
  b[1] = cf->lam[0];
  b[2] = cf->lam[1];
}

dbl utetra_get_value(utetra_s const *cf) {
  return cf->f;
}

void utetra_get_gradient(utetra_s const *cf, dbl g[2]) {
  g[0] = cf->g[0];
  g[1] = cf->g[1];
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
void utetra_get_lag_mults(utetra_s const *cf, dbl alpha[3]) {
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
    assert(b[0] > -atol && b[1] > -atol && b[2] > -atol);
  }
}

int utetra_get_num_iter(utetra_s const *cf) {
  return cf->niter;
}

bool utetra_has_interior_point_solution(utetra_s const *cf) {
  dbl alpha[3];
  utetra_get_lag_mults(cf, alpha);
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
    dbl T1 = u1->f, T2 = u2->f;
    if (T1 < T2) {
      return -1;
    } else if (T1 > T2) {
      return 1;
    } else {
      return 0;
    }
  }
}

bool utetra_adj_are_optimal(utetra_s const *u1, utetra_s const *u2) {
  dbl const atol = 1e-15;
  dbl lam1[2], lam2[2];
  utetra_get_lambda(u1, lam1);
  utetra_get_lambda(u2, lam2);
  return fabs(lam1[0] - lam2[0]) <= atol
    && fabs(lam1[1]) <= atol
    && fabs(lam2[1]) <= atol
    && fabs(utetra_get_value(u1) - utetra_get_value(u2)) <= atol;
}

void utetra_get_point_on_ray(utetra_s const *utetra, dbl t, dbl xt[3]) {
  // TODO: optimize this by using utetra->x instead of computing xb
  dbl b[3], xb[3], L;
  utetra_get_bary_coords(utetra, b);
  dbl33_dbl3_mul(utetra->X, b, xb);
  L = dbl3_norm(utetra->x_minus_xb);
  dbl3_saxpy(t/L, utetra->x_minus_xb, xb, xt);
}

int utetra_get_interior_coefs(utetra_s const *utetra, bool I[3]) {
  dbl const atol = 1e-14;
  I[0] = utetra->lam[0] + utetra->lam[1] < 1 - atol;
  I[1] = utetra->lam[0] > atol;
  I[2] = utetra->lam[1] > atol;
  return I[0] + I[1] + I[2];
}

bool utetra_update_ray_is_physical(utetra_s const *utetra, eik3_s const *eik) {
  size_t const *l = utetra->l;

  mesh3_s const *mesh = eik3_get_mesh(eik);

  // TODO: it's pretty hard to say where this is what we want to do or
  // not... let's see how it goes!
  if (mesh3_bdf(mesh, l))
    return false;

  // TODO: the following section where we check to see if the stuff
  // below gives "an interior ray" can be wrapped up and reused for
  // both this and the corresponding section in utri.c...

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
  dbl xm[3], xp[3], t = mesh3_get_min_tetra_alt(mesh)/2;
  utetra_get_point_on_ray(utetra, -t, xm);
  utetra_get_point_on_ray(utetra, t, xp);

  // Find the number and location of interior coefficients.
  bool I[3];
  utetra_get_interior_coefs(utetra, I);

  bool xm_in_cell = false, xp_in_cell = false;

  array_s *cells;
  array_alloc(&cells);
  array_init(cells, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  int nvc;
  size_t *vc;

  for (int i = 0; i < 3; ++i) {
    if (!I[i]) continue;

    nvc = mesh3_nvc(mesh, l[i]);
    vc = malloc(nvc*sizeof(size_t));
    mesh3_vc(mesh, l[i], vc);

    for (int j = 0; j < nvc; ++j)
      if (!array_contains(cells, &vc[j]))
        array_append(cells, &vc[j]);

    free(vc);
  }

  dbl b[4];
  size_t lc;
  for (size_t i = 0; i < array_size(cells); ++i) {
    array_get(cells, i, &lc);
    xm_in_cell |= mesh3_dbl3_in_cell(mesh, lc, xm, b);
    xp_in_cell |= mesh3_dbl3_in_cell(mesh, lc, xp, b);
    if (xm_in_cell && xp_in_cell)
      break;
  }

  // If we didn't find a containing cell, we can conclude that the ray
  // is unphysical!
  if (!xm_in_cell || !xp_in_cell)
    return false;

  /**
   * Next, check and see if the point just before the end of the ray
   * lies in a cell.
   */

  size_t lhat = utetra->lhat;

  dbl xhatm[3];
  dbl3_saxpy(-t/utetra->L, utetra->x_minus_xb, utetra->x, xhatm);

  nvc = mesh3_nvc(mesh, lhat);
  vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(mesh, lhat, vc);

  bool xhatm_in_cell = false;
  for (int i = 0; i < nvc; ++i) {
    xhatm_in_cell = mesh3_dbl3_in_cell(mesh, vc[i], xhatm, b);
    if (xhatm_in_cell)
      break;
  }

  free(vc);

  return xhatm_in_cell;
}

int utetra_get_num_interior_coefs(utetra_s const *utetra) {
  dbl const atol = 1e-14;
  return (utetra->lam[0] > atol) + (utetra->lam[1] > atol)
    + (utetra->lam[0] + utetra->lam[1] < 1 - atol);
}

void utetra_get_update_inds(utetra_s const *utetra, size_t l[3]) {
  memcpy(l, utetra->l, sizeof(size_t[3]));
}

bool utetra_has_shadow_solution(utetra_s const *utetra, eik3_s const *eik) {
  dbl const atol = 1e-14;

  dbl b[3];
  utetra_get_bary_coords(utetra, b);

  bool has_shadow_solution = true;

  // Check if active => shadow
  for (int i = 0; i < 3; ++i)
    if (fabs(b[i]) > atol)
      has_shadow_solution &= eik3_is_shadow(eik, utetra->l[i]);

  return has_shadow_solution;
}

bool utetras_yield_same_update(utetra_s const **utetra, int n) {
  dbl const atol = 1e-14;

  dbl b[2][3];
  jet3 jet[2];

  // Prefetch the first coords and jet
  utetra_get_bary_coords(utetra[0], b[0]);
  utetra_get_jet(utetra[0], &jet[0]);

  for (int i = 1; i < n; ++i) {
    // Get the next jet and check that it's finite
    utetra_get_jet(utetra[i], &jet[1]);
    if (!jet3_is_finite(&jet[1]))
      return false;

    // Get the next coords
    utetra_get_bary_coords(utetra[i], b[1]);

    // Check if the coords and jet agree up to `atol` and return early
    // if they don't
    if (fabs(b[0][0] - b[1][0]) > atol ||
        fabs(b[0][1] - b[1][1]) > atol ||
        fabs(b[0][2] - b[1][2]) > atol ||
        !jet3_approx_eq(&jet[0], &jet[1], atol))
      return false;

    // Swap the coords and jet that we just fetched to make way for
    // the next ones
    for (int j = 0; j < 3; ++j) SWAP(b[0][j], b[1][j]);
    SWAP(jet[0], jet[1]);
  }

  return true;
}
