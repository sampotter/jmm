#include <jmm/uline.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <jmm/eik3.h>
#include <jmm/log.h>
#include <jmm/mesh3.h>

static size_t MAX_NUM_ITER = 100;

struct uline {
  /* Line update parameters: */
  eik3_s const *eik;
  size_t lhat;
  size_t l0;
  dbl tol;

  stype_e stype;
  sfunc_s const *sfunc;

  /* Internal parameters used to compute update: */
  dbl3 xm;
  dbl3 x0, xhat;
  dbl3 phipm;
  dbl L;
  dbl T0;
  dbl s0, shat;

  /* Computed values: */
  dbl f;
  dbl3 gradf;
};

void uline_alloc(uline_s **u) {
  *u = malloc(sizeof(uline_s));
}

void uline_dealloc(uline_s **u) {
  free(*u);
  *u = NULL;
}

static void set_internal_params(uline_s *u) {
  dbl3_nan(u->xm);

  dbl3_sub(u->xhat, u->x0, u->phipm);

  u->L = dbl3_norm(u->phipm);

  dbl3_dbl_div_inplace(u->phipm, u->L);

  if (u->stype == STYPE_CONSTANT) {
    u->s0 = 1;
    u->shat = 1;
  } else if (u->stype == STYPE_FUNC_PTR) {
    u->s0 = u->sfunc->funcs.s(u->x0);
    u->shat = u->sfunc->funcs.s(u->xhat);
  } else {
    assert(false); // TODO: implement
  }
}

void uline_init(uline_s *u, eik3_s const *eik, size_t lhat, size_t l0) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  u->eik = eik;
  u->stype = eik3_get_stype(u->eik);
  u->sfunc = eik3_get_sfunc(u->eik);

  u->lhat = lhat;
  u->l0 = l0;

  mesh3_copy_vert(mesh, u->lhat, u->xhat);
  mesh3_copy_vert(mesh, u->l0, u->x0);

  u->tol = mesh3_get_vertex_tol(mesh, l0);

  u->T0 = eik3_get_T(u->eik, u->l0);

  set_internal_params(u);
}

void uline_init_from_points(uline_s *u, eik3_s const *eik, dbl3 const xhat, dbl3 const x0, dbl tol, dbl T0) {
  u->eik = eik;
  u->stype = eik3_get_stype(u->eik);
  u->sfunc = eik3_get_sfunc(u->eik);

  u->lhat = (size_t)NO_INDEX;
  u->l0 = (size_t)NO_INDEX;

  dbl3_copy(xhat, u->xhat);
  dbl3_copy(x0, u->x0);

  u->tol = tol;

  u->T0 = T0;

  set_internal_params(u);
}

void set_xm(uline_s *u, dbl3 const xm) {
  assert(u->stype == STYPE_FUNC_PTR);

  dbl3_copy(xm, u->xm);

  dbl3 phip0;
  for (size_t i = 0; i < 3; ++i)
    phip0[i] = (-3*u->x0[i] + 4*u->xm[i] - u->xhat[i])/u->L;

  dbl3 phipL;
  for (size_t i = 0; i < 3; ++i)
    phipL[i] = (u->x0[i] - 4*u->xm[i] + 3*u->xhat[i])/u->L;

  dbl phip0_norm = dbl3_norm(phip0);
  dbl phipL_norm = dbl3_norm(phipL);

  dbl sm = u->sfunc->funcs.s(u->xm);

  dbl3 gradsm;
  u->sfunc->funcs.Ds(u->xm, gradsm);

  /* NOTE: norm of u->phipm is 1 */
  u->f = u->T0 + (u->L/6)*(u->s0*phip0_norm + 4*sm + u->shat*phipL_norm);

  for (size_t i = 0; i < 3; ++i)
    u->gradf[i] = (2./3)*(
      u->L*gradsm[i]
      + u->s0*phip0[i]/phip0_norm - u->shat*phipL[i]/phipL_norm);
}

static void solve_stype_constant(uline_s *u) {
  assert(u->stype == STYPE_CONSTANT);

  u->f = u->T0 + dbl3_dist(u->xhat, u->x0);
  dbl3_nan(u->gradf);
}

static void solve_stype_func_ptr(uline_s *u) {
  assert(u->stype == STYPE_FUNC_PTR);

  dbl3 xm0, xm, gradf0;
  dbl alpha, f0;
  size_t num_iter_outer, num_iter_inner;

  /* Initialize xm to be the midpoint between x0 and xhat */
  dbl3_avg(u->x0, u->xhat, xm0);
  set_xm(u, xm0);

  /* Gradient descent on xm with backtracking line search */
  num_iter_outer = 0;
  do {
    f0 = u->f;
    dbl3_copy(u->gradf, gradf0);
    num_iter_inner = 0;
    alpha = 1;
  line_search:
    dbl3_saxpy(-alpha, gradf0, xm0, xm);
    set_xm(u, xm);
    if (num_iter_inner < MAX_NUM_ITER && u->f >= f0) {
      alpha /= 2;
      ++num_iter_inner;
      goto line_search;
    }
    dbl3_copy(xm, xm0);
    ++num_iter_outer;
  } while (dbl3_norm(u->gradf) > u->tol);
}

void uline_solve(uline_s *u) {
  if (u->stype == STYPE_CONSTANT) {
    solve_stype_constant(u);
  } else if (u->stype == STYPE_FUNC_PTR) {
    solve_stype_func_ptr(u);
  } else {
    assert(false);
  }
}

dbl uline_get_value(uline_s const *u) {
  return u->f;
}

static void get_topt_stype_constant(uline_s const *u, dbl3 topt) {
  assert(u->stype == STYPE_CONSTANT);

  dbl3_sub(u->xhat, u->x0, topt);
  dbl3_normalize(topt);
}

static void get_topt_stype_func_ptr(uline_s const *u, dbl3 topt) {
  assert(u->stype == STYPE_FUNC_PTR);

  dbl3 phipL;
  for (size_t i = 0; i < 3; ++i)
    phipL[i] = (u->x0[i] - 4*u->xm[i] + 3*u->xhat[i])/u->L;
  dbl3_normalized(phipL, topt);
}

void uline_get_topt(uline_s const *u, dbl3 topt) {
  if (u->stype == STYPE_CONSTANT) {
    get_topt_stype_constant(u, topt);
  } else if (u->stype == STYPE_FUNC_PTR) {
    get_topt_stype_func_ptr(u, topt);
  } else {
    assert(false);
  }
}

jet31t uline_get_jet(uline_s const *u) {
  jet31t jet;
  jet.f = u->f;
  uline_get_topt(u, jet.Df);
  dbl3_dbl_mul_inplace(jet.Df, u->shat);
  return jet;
}
