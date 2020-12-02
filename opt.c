#include "opt.h"

#include <assert.h>
#include <math.h>

#include "hybrid.h"
#include "vec.h"

void eqp_bary_2_0(dbl const *G, dbl const *c, dbl *x) {
  x[0] = 0.0;
  x[1] = -c[1]/G[2];
}

void eqp_bary_2_1(dbl const *G, dbl const *c, dbl *x) {
  x[0] = -c[0]/G[0];
  x[1] = 0.0;
}

void eqp_bary_2_2(dbl const *G, dbl const *c, dbl *x) {
  dbl pZ = ((G[2] - G[0])/2 + c[1] - c[0])/(G[0] - 2*G[1] + G[2]);
  x[0] = 0.5 + pZ;
  x[1] = 0.5 - pZ;
}

static void compute_p(dbl const *G, dbl const *c, dbl const *x, dbl *p) {
  dbl det = G[0]*G[2] - G[1]*G[1];
  p[0] = -(x[0] + (G[2]*c[0] - G[1]*c[1])/det);
  p[1] = -(x[1] + (G[0]*c[1] - G[1]*c[0])/det);
}

static void compute_y(dbl const *G, dbl const *c, dbl const *x, dbl *y) {
  y[0] = G[0]*x[0] + G[1]*x[1] + c[0];
  y[1] = G[1]*x[0] + G[2]*x[1] + c[1];
}

static int get_num_active(bool const *active) {
  return (int)((active[0] ? 1 : 0) +
               (active[1] ? 1 : 0) +
               (active[2] ? 1 : 0));
}

void iqp_bary_2(dbl const *G, dbl const *c, dbl const *x0, dbl *x, bool *error,
                dbl tol, int niters) {
  assert(x != NULL);

  if (error) {
    *error = false;
  }

  dbl xprev[2], p[2] = {0, 0}, y[2], alpha, alpha_new;
  if (x0) {
    x[0] = x0[0];
    x[1] = x0[1];
  } else {
    x[0] = x[1] = 0.0;
  }

  bool active[3], ind[3];
  active[0] = fabs(x[0]) <= tol;
  active[1] = fabs(x[1]) <= tol;
  active[2] = fabs(1 - x[0] - x[1]) <= tol;

  int k = 0, num_active, argmin;
  while (true) {
    xprev[0] = x[0];
    xprev[1] = x[1];
    num_active = get_num_active(active);

    compute_y(G, c, x, y);

    if (num_active == 0) {
      compute_p(G, c, x, p);
    } else if (num_active == 1) {
      if (active[0]) eqp_bary_2_0(G, c, p);
      else if (active[1]) eqp_bary_2_1(G, c, p);
      else if (active[2]) eqp_bary_2_2(G, c, p);
      else assert(false);
      p[0] -= xprev[0];
      p[1] -= xprev[1];
    } else if (num_active == 2) {
      p[0] = p[1] = 0.0;
    }

    if (fmax(fabs(p[0]), fabs(p[1])) <= tol) {
      if (num_active == 0) {
        break;
      } else if (num_active == 1) {
        if (active[0]) {
          if (y[0] >= 0) break;
          else active[0] = false;
        } else if (active[1]) {
          if (y[1] >= 0) break;
          else active[1] = false;
        } else if (active[2]) {
          if (y[0] + y[1] <= 0) break;
          else active[2] = false;
        } else {
          assert(false);
        }
      } else if (num_active == 2) {
        if (active[0] && active[1]) {
          if (y[0] >= 0 && y[1] >= 0) break;
          else active[y[0] < y[1] ? 0 : 1] = false;
        } else if (active[0] && active[2]) {
          if (y[1] <= 0 && y[1] <= y[0]) break;
          else active[y[0] < 0 ? 0 : 2] = false;
        } else if (active[1] && active[2]) {
          if (y[0] <= 0 && y[0] <= y[1]) break;
          else active[y[1] < 0 ? 1 : 2] = false;
        } else {
          assert(false);
        }
      } else {
        assert(false);
      }
    } else {
      alpha = 1;
      argmin = -1;
      ind[0] = !active[0] && p[0] < 0;
      ind[1] = !active[1] && p[1] < 0;
      ind[2] = !active[2] && p[0] + p[1] > 0;
      if (ind[0] || ind[1] || ind[2]) {
        if (ind[0]) {
          alpha_new = -x[0]/p[0];
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 0;
          }
        }
        if (ind[1]) {
          alpha_new = -x[1]/p[1];
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 1;
          }
        }
        if (ind[2]) {
          alpha_new = (1 - x[0] - x[1])/(p[0] + p[1]);
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 2;
          }
        }
        alpha = fmax(0.0, fmin(alpha, 1.0));
        if (alpha < 1) {
          active[argmin] = true;
        }
      }
      x[0] += alpha*p[0];
      x[1] += alpha*p[1];
    }

    if (++k == niters) {
      if (error) *error = true;
      break;
    }
  }
}

dbl get_min_eigval(dbl const *A) {
  dbl half_tr = (A[0] + A[2])/2;
  dbl det = A[0]*A[2] - A[1]*A[1];
  return half_tr - sqrt(half_tr*half_tr - det);
}

void sqp_bary_3_2(costfunc_s const *cf, dbl const *xinit, dbl *x,
                  dbl *f, bool *error, dbl tol, int niters)
{
  if (error) *error = false;

  int k = 0;
  dbl x0[2], x1[2], f0, f1, lambda_min;
  dbl d2f[3], df[2], h[2], c[2], alpha;
  bool iqp_error;

  x1[0] = xinit ? xinit[0] : 1./3;
  x1[1] = xinit ? xinit[1] : 1./3;

  cf->func(x1, &f1, cf->wkspc);

  while (true) {
    // Compute Hessian
    cf->hess(x1, d2f, cf->wkspc);

    // ... perturb it if it isn't positive definite
    lambda_min = get_min_eigval(d2f);
    if (lambda_min < 0) {
      d2f[0] -= 1.1*lambda_min;
      d2f[2] -= 1.1*lambda_min;
    }

    cf->grad(x1, df, cf->wkspc);

    c[0] = df[0] - d2f[0]*x1[0] - d2f[1]*x1[1];
    c[1] = df[1] - d2f[1]*x1[0] - d2f[2]*x1[1];

    x0[0] = x1[0];
    x0[1] = x1[1];

    iqp_bary_2(d2f, c, x0, x1, &iqp_error, tol, 10);
    assert(!iqp_error);

    // Compute descent step: h = x1 - x0 = x_{n+1} - x_{n}
    dbl2_sub(x1, x0, h);

    alpha = 1;
    dbl lhs, rhs = f1 + 1e-4*alpha*dbl2_dot(df, h);
  repeat:
    dbl2_saxpy(alpha, h, x0, x1);
    cf->func(x1, &lhs, cf->wkspc);
    if (alpha > tol && lhs > rhs) {
      alpha /= 2;
      goto repeat;
    }

    dbl2_saxpy(alpha, h, x0, x1);
    f0 = f1;
    cf->func(x1, &f1, cf->wkspc);

    if (fabs(f1 - f0) <= tol*fmax(f0, f1) + tol) {
      break;
    }

    if (fmax(fabs(x1[0] - x0[0]), fabs(x1[1] - x0[1])) <=
        tol*fmax(fmax(x0[0], x0[1]), fmax(x1[0], x1[1])) + tol) {
      break;
    }

    // Check if we've reached our max number of iterations
    if (++k == niters) {
      if (error) *error = true;
      break;
    }
  }

  if (x != NULL) {
    x[0] = x1[0];
    x[1] = x1[1];
  }

  if (f != NULL) {
    *f = f1;
  }
}
