#include "update.h"

#include <assert.h>
#include <string.h>

#include "bb.h"
#include "mat.h"
#include "mesh3.h"
#include "opt.h"

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

  dbl T[3];

  // B-coefs for 9-point triangle interpolation T on base of update
  dbl Tc[10];

  dbl x_minus_xb[3];

  int niter;
};

void utetra_alloc(utetra_s **utetra) {
  *utetra = malloc(sizeof(utetra_s));
}

void utetra_dealloc(utetra_s **utetra) {
  free(*utetra);
  *utetra = NULL;
}

void utetra_init(utetra_s *cf, mesh3_s const *mesh, jet3 const *jet,
                   size_t l, size_t l0, size_t l1, size_t l2) {
  assert(jet3_is_finite(&jet[l0]));
  assert(jet3_is_finite(&jet[l1]));
  assert(jet3_is_finite(&jet[l2]));

  mesh3_get_vert(mesh, l, cf->x);

  mesh3_get_vert(mesh, l0, cf->Xt[0]);
  mesh3_get_vert(mesh, l1, cf->Xt[1]);
  mesh3_get_vert(mesh, l2, cf->Xt[2]);

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

  /**
   * Compute Bernstein-Bezier coefficients before transposing Xt and
   * computing XtX
   */

  dbl DT[3][3];

  cf->T[0] = jet[l0].f;
  DT[0][0] = jet[l0].fx;
  DT[0][1] = jet[l0].fy;
  DT[0][2] = jet[l0].fz;

  cf->T[1] = jet[l1].f;
  DT[1][0] = jet[l1].fx;
  DT[1][1] = jet[l1].fy;
  DT[1][2] = jet[l1].fz;

  cf->T[2] = jet[l2].f;
  DT[2][0] = jet[l2].fx;
  DT[2][1] = jet[l2].fy;
  DT[2][2] = jet[l2].fz;

  bb3tri_interp3(cf->T, &DT[0], cf->Xt, cf->Tc);

  cf->niter = 0;
}

bool utetra_is_causal(utetra_s const *cf) {
  return cf->angles[0] >= 0 && cf->angles[1] >= 0 && cf->angles[2] >= 0;
}

void utetra_reset(utetra_s *cf) {
  cf->niter = 0;
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

  assert(cf->f > cf->T[0]);
  assert(cf->f > cf->T[1]);
  assert(cf->f > cf->T[2]);
}

void utetra_get_lambda(utetra_s *cf, dbl lam[2]) {
  lam[0] = cf->lam[0];
  lam[1] = cf->lam[1];
}

void utetra_set_lambda(utetra_s *cf, dbl const lam[2]) {
  static dbl a1[3] = {-1, 1, 0};
  static dbl a2[3] = {-1, 0, 1};

  static dbl const atol = 1e-15;

  dbl b[3], xb[3], tmp1[3], tmp2[3][3], L, DL[2], D2L[2][2], DT[2], D2T[2][2];
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
  L = dbl3_norm(cf->x_minus_xb);

  dbl33_dbl3_mul(cf->Xt, cf->x_minus_xb, tmp1);
  dbl3_dbl_div(tmp1, -L, tmp1);

  DL[0] = dbl3_dot(a1, tmp1);
  DL[1] = dbl3_dot(a2, tmp1);

  dbl3_outer(tmp1, tmp1, tmp2);
  dbl33_sub(cf->XtX, tmp2, tmp2);
  dbl33_dbl_div(tmp2, L, tmp2);

  dbl33_dbl3_mul(tmp2, a1, tmp1);
  D2L[0][0] = dbl3_dot(tmp1, a1);
  D2L[1][0] = D2L[0][1] = dbl3_dot(tmp1, a2);
  dbl33_dbl3_mul(tmp2, a2, tmp1);
  D2L[1][1] = dbl3_dot(tmp1, a2);

  DT[0] = dbb3tri(cf->Tc, b, a1);
  DT[1] = dbb3tri(cf->Tc, b, a2);

  D2T[0][0] = d2bb3tri(cf->Tc, b, a1, a1);
  D2T[1][0] = D2T[0][1] = d2bb3tri(cf->Tc, b, a1, a2);
  D2T[1][1] = d2bb3tri(cf->Tc, b, a2, a2);

  cf->f = L + bb3tri(cf->Tc, b);
  dbl2_add(DL, DT, cf->g);
  dbl22_add(D2L, D2T, cf->H);

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

dbl utetra_get_value(utetra_s const *cf) {
  return cf->f;
}

void utetra_get_gradient(utetra_s const *cf, dbl g[2]) {
  g[0] = cf->g[0];
  g[1] = cf->g[1];
}

void utetra_get_jet(utetra_s const *cf, jet3 *jet) {
  jet->f = cf->f;
  dbl L = dbl3_norm(cf->x_minus_xb);
  jet->fx = cf->x_minus_xb[0]/L;
  jet->fy = cf->x_minus_xb[1]/L;
  jet->fz = cf->x_minus_xb[2]/L;
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
