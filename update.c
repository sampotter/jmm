#include "update.h"

#include <assert.h>
#include <string.h>

#include "bb.h"
#include "mat.h"
#include "mesh3.h"

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

  /**
   * Compute Bernstein-Bezier coefficients before transposing Xt and
   * computing XtX
   */

  dbl T[3];
  dbl DT[3][3];

  T[0] = jet[l0].f;
  DT[0][0] = jet[l0].fx;
  DT[0][1] = jet[l0].fy;
  DT[0][2] = jet[l0].fz;

  T[1] = jet[l1].f;
  DT[1][0] = jet[l1].fx;
  DT[1][1] = jet[l1].fy;
  DT[1][2] = jet[l1].fz;

  T[2] = jet[l2].f;
  DT[2][0] = jet[l2].fx;
  DT[2][1] = jet[l2].fy;
  DT[2][2] = jet[l2].fz;

  bb3tri_interp3(T, &DT[0], cf->Xt, cf->Tc);
}

void utetra_set_lambda(utetra_s *cf, dbl const *lambda) {
  static dbl a1[3] = {-1, 1, 0};
  static dbl a2[3] = {-1, 0, 1};

  static dbl const atol = 5e-16;

  dbl b[3], xb[3], tmp1[3], tmp2[3][3], L, DL[2], D2L[2][2], DT[2], D2T[2][2];

  b[1] = lambda[0];
  b[2] = lambda[1];
  b[0] = 1 - b[1] - b[2];

  assert(b[0] >= -atol);
  assert(b[1] >= -atol);
  assert(b[2] >= -atol);

  b[0] = fmax(0.0, b[0]);
  b[1] = fmax(0.0, b[1]);
  b[2] = fmax(0.0, b[2]);

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
   * Finally, compute Newton step, making sure to perturb the Hessian
   * if it's indefinite.
   */

  // Conditionally perturb the Hessian
  dbl tr = cf->H[0][0] + cf->H[1][1];
  dbl det = cf->H[0][0]*cf->H[1][1] - cf->H[0][1]*cf->H[1][0];
  dbl min_eig_doubled = tr - sqrt(tr*tr - 4*det);
  if (min_eig_doubled < 0) {
    cf->H[0][0] -= min_eig_doubled;
    cf->H[1][1] -= min_eig_doubled;
  }

  // Compute the Newton step
  dbl22_dbl2_solve(cf->H, cf->g, cf->p);
  dbl2_negate(cf->p);

  // This may get overwritten below.
  cf->g_dot_p_active = dbl2_dot(cf->g, cf->p);

  /**
   * Project Newton step based on active constraints
   */

  if (b[0] <= atol && cf->p[0] + cf->p[1] > 0) {
    dbl newp[2] = {
      (cf->p[0] - cf->p[1])/2,
      (cf->p[1] - cf->p[0])/2
    };
    cf->p[0] = newp[0];
    cf->p[1] = newp[1];
    // cf->g_dot_p_active = -dbl2_norm_sq(cf->g);
    // cf->g_dot_p_active -= dbl22_trace(cf->H) - cf->H[0][1] - cf->H[1][0];
  }

  if (b[1] <= atol) {
    cf->p[0] = fmax(0.0, cf->p[0]);
    // cf->g_dot_p_active = -cf->g[1]*cf->g[1]/cf->H[1][1];
  }

  if (b[2] <= atol) {
    cf->p[1] = fmax(0.0, cf->p[1]);
    // cf->g_dot_p_active = -cf->g[0]*cf->g[0]/cf->H[0][0];
  }
}

/**
 * Do a tetrahedron update starting at `lam`, writing the result to
 * `jet`. This assumes that `utetra_set_lambda` has already been
 * called, so that `cf` is currently at `lam`.
 */
void utetra_solve(utetra_s *cf, dbl lam[2], jet3 *jet) {
  dbl const tscale = 0.5; // Step size scaling parameter
  dbl const c1 = 1e-4; // Constant for backtracking line search
  dbl const rtol = 1e-15, atol = 1e-15;

  dbl lam1[2], dlam[2];
  dbl f = cf->f; // Current value of cost function
  dbl t = 1; // Initial step size
  dbl tc; // Breakpoint used to find Cauchy point
  dbl c1_times_g_dot_p;
  dbl Df; // Directional derivative used in Cauchy point calculation
  dbl denom; // Denominator used in Cauchy point calculations

  // See page 16 of Kelley
  dbl tol = rtol*dbl2_norm(cf->g) + atol;

  cf->niter = 0;

  /**
   * Newton iteration
   *
   * TODO: in most of the places we call set_lambda below, we only
   * need some of the stuff computed by set_lambda (usually just
   * cf->f)... definitely don't want to waste time computing extra
   * quantities!
   */
  while (dbl2_maxnorm(cf->p) >= tol) {
    // c1_times_g_dot_p = c1*dbl2_dot(cf->g, cf->p);
    c1_times_g_dot_p = c1*cf->g_dot_p_active;
    assert(c1_times_g_dot_p < 0);

    /**
     * Find the Cauchy point
     */

    denom = dbl2_sum(cf->p);
    tc = (1 - lam[0] - lam[1])/denom;
    if (fabs(denom) > rtol && 0 < tc && tc < 1) {
      dbl2_saxpy(tc, cf->p, lam, lam1);
      utetra_set_lambda(cf, lam1);
      dbl2_sub(lam1, lam, dlam);
      Df = dbl2_dot(cf->g, dlam)/tc;
      if (Df >= 0) {
        t = tc;
        goto backtrack;
      } else {
        goto reset;
      }
    }

    denom = cf->p[0];
    tc = -lam[0]/denom;
    if (fabs(denom) > rtol && 0 < tc && tc < 1) {
      dbl2_saxpy(tc, cf->p, lam, lam1);
      utetra_set_lambda(cf, lam1);
      dbl2_sub(lam1, lam, dlam);
      Df = dbl2_dot(cf->g, dlam)/tc;
      if (Df >= 0) {
        t = tc;
        goto backtrack;
      } else {
        goto reset;
      }
    }

    denom = cf->p[1];
    tc = -lam[1]/denom;
    if (fabs(denom) > rtol && 0 < tc && tc < 1) {
      dbl2_saxpy(tc, cf->p, lam, lam1);
      utetra_set_lambda(cf, lam1);
      dbl2_sub(lam1, lam, dlam);
      Df = dbl2_dot(cf->g, dlam)/tc;
      if (Df >= 0) {
        t = tc;
        goto backtrack;
      } else {
        goto reset;
      }
    }

    // We didn't trip any breakpoints, compute lam1
    dbl2_saxpy(t, cf->p, lam, lam1);
    utetra_set_lambda(cf, lam1);

  backtrack:
    while (cf->f > f + t*c1_times_g_dot_p + atol) {
      t *= tscale;
      dbl2_saxpy(t, cf->p, lam, lam1);
      utetra_set_lambda(cf, lam1);
    }

  reset: // Reset for next iteration
    t = 1;
    lam[0] = lam1[0];
    lam[1] = lam1[1];
    f = cf->f;
    ++cf->niter;
  }

  dbl DT[3];
  memcpy((void *)DT, (void *)cf->x_minus_xb, sizeof(dbl)*3);
  dbl3_normalize(DT);

  jet->f = cf->f;
  jet->fx = DT[0];
  jet->fy = DT[1];
  jet->fz = DT[2];
}
