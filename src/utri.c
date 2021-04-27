#include "utri.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "array.h"
#include "bb.h"
#include "eik3.h"
#include "hybrid.h"
#include "mesh3.h"
#include "vec.h"

utri_spec_s utri_spec_empty() {
  return (utri_spec_s) {
    .eik = NULL,
    .lhat = (size_t)NO_INDEX,
    .l = {(size_t)NO_INDEX, (size_t)NO_INDEX},
    .state = {UNKNOWN, UNKNOWN},
    .xhat = {NAN, NAN, NAN},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
    },
    .orig_index = (size_t)NO_INDEX
  };
}

utri_spec_s utri_spec_from_eik(eik3_s const *eik, size_t l, size_t l0, size_t l1) {
  state_e const *state = eik3_get_state_ptr(eik);
  return (utri_spec_s) {
    .eik = eik,
    .lhat = l,
    .l = {l0, l1},
    .state = {state[l0], state[l1]},
    .xhat = {NAN, NAN, NAN},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
    }
  };
}

utri_spec_s utri_spec_from_eik_without_l(eik3_s const *eik, dbl const x[3],
                                         size_t l0, size_t l1) {
  state_e const *state = eik3_get_state_ptr(eik);
  return (utri_spec_s) {
    .eik = eik,
    .lhat = (size_t)NO_INDEX,
    .l = {l0, l1},
    .state = {state[l0], state[l1]},
    .xhat = {x[0], x[1], x[2]},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
    }
  };
}

utri_spec_s utri_spec_from_raw_data(dbl const x[3], dbl const Xt[2][3], jet3 jet[2]) {
  return (utri_spec_s) {
    .eik = NULL,
    .lhat = NO_INDEX,
    .l = {NO_INDEX, NO_INDEX},
    .state = {UNKNOWN, UNKNOWN},
    .xhat = {x[0], x[1], x[2]},
    .x = {
      {Xt[0][0], Xt[0][1], Xt[0][2]},
      {Xt[1][0], Xt[1][1], Xt[1][2]}
    },
    .jet = {jet[0], jet[1]}
  };
}

struct utri {
  dbl lam;
  dbl f;
  dbl Df;
  dbl x_minus_xb[3];
  dbl L;

  size_t l, l0, l1;
  dbl x[3];
  dbl x0[3];
  dbl x1[3];
  dbl x1_minus_x0[3];
  state_e state[2];
  bb31 T;

  size_t num_shadow;

  size_t orig_index; // original index
};

void utri_alloc(utri_s **utri) {
  *utri = malloc(sizeof(utri_s));
}

void utri_dealloc(utri_s **utri) {
  free(*utri);
  *utri = NULL;
}

void utri_set_lambda(utri_s *utri, dbl lam) {
  utri->lam = lam;

  dbl xb[3];
  dbl3_saxpy(lam, utri->x1_minus_x0, utri->x0, xb);
  dbl3_sub(utri->x, xb, utri->x_minus_xb);
  utri->L = dbl3_norm(utri->x_minus_xb);

  dbl dL_dlam = -dbl3_dot(utri->x1_minus_x0, utri->x_minus_xb)/utri->L;

  dbl b[2] = {1 - lam, lam};
  dbl T = bb31_f(&utri->T, b);

  utri->f = T + utri->L;

  dbl a[2] = {-1, 1};
  dbl dT_dlam = bb31_df(&utri->T, b, a);

  utri->Df = dT_dlam + dL_dlam;
}

static dbl utri_hybrid_f(dbl lam, utri_s *utri) {
  utri_set_lambda(utri, lam);
  return utri->Df;
}

void utri_init(utri_s *u, utri_spec_s const *spec) {
  bool passed_lhat = spec->lhat != (size_t)NO_INDEX;
  bool passed_l0 = spec->l[0] != (size_t)NO_INDEX;
  bool passed_l1 = spec->l[1] != (size_t)NO_INDEX;
  bool passed_l = passed_l0 && passed_l1;

  bool passed_jet0 = spec->state[0] == VALID && jet3_is_finite(&spec->jet[0]);
  bool passed_jet1 = spec->state[1] == VALID && jet3_is_finite(&spec->jet[1]);
  bool passed_jet = passed_jet0 && passed_jet1;

#if JMM_DEBUG
  /* Validate spec before doing anything else */

  bool passed_xhat = dbl3_isfinite(spec->xhat);
  assert(passed_lhat ^ passed_xhat); // exactly one of these

  if (passed_l0 || passed_l1)
    assert(passed_l);

  bool passed_x0 = dbl3_isfinite(spec->x[0]);
  bool passed_x1 = dbl3_isfinite(spec->x[1]);
  bool passed_x = passed_x0 && passed_x1;
  if (passed_x0 || passed_x1)
    assert(passed_x);

  assert(passed_l ^ passed_x); // pass exactly one of these

  bool passed_state0 = spec->state[0] != UNKNOWN;
  bool passed_state1 = spec->state[1] != UNKNOWN;
  bool passed_state = passed_state0 && passed_state1;
  if (passed_state0 || passed_state1)
    assert(passed_state);

  if (passed_jet0 || passed_jet1)
    assert(passed_jet);

  assert(passed_jet ^ passed_l); // exactly one of these
#endif

  /* Initialize `u` */

  u->f = INFINITY;

  u->lam = u->f = u->Df = u->L = NAN;
  u->x_minus_xb[0] = u->x_minus_xb[1] = u->x_minus_xb[2] = NAN;

  mesh3_s const *mesh = spec->eik ? eik3_get_mesh(spec->eik) : NULL;

  u->l = spec->lhat;

  if (passed_lhat) {
    assert(mesh);
    mesh3_copy_vert(mesh, u->l, u->x);
  } else {
    dbl3_copy(spec->xhat, u->x);
  }

  u->l0 = spec->l[0];
  u->l1 = spec->l[1];

  if (passed_l) {
    mesh3_copy_vert(mesh, u->l0, u->x0);
    mesh3_copy_vert(mesh, u->l1, u->x1);
  } else { // passed_x
    dbl3_copy(spec->x[0], u->x0);
    dbl3_copy(spec->x[1], u->x1);
  }

  dbl3_sub(u->x1, u->x0, u->x1_minus_x0);

  memcpy(u->state, spec->state, sizeof(state_e[2]));

  jet3 jet[2];
  for (size_t i = 0; i < 2; ++i)
    jet[i] = passed_jet ? spec->jet[i] : eik3_get_jet(spec->eik, spec->l[i]);

  dbl T[2], DT[2][3];
  for (int i = 0; i < 2; ++i) {
    T[i] = jet[i].f;
    memcpy(DT[i], &jet[i].fx, sizeof(dbl[3]));
  }

  dbl Xt[2][3];
  dbl3_copy(u->x0, Xt[0]);
  dbl3_copy(u->x1, Xt[1]);
  bb31_init_from_3d_data(&u->T, T, DT, Xt);

  u->orig_index = spec->orig_index;
}

void utri_deinit(utri_s *u) {
  (void)u;
}

void utri_solve(utri_s *utri) {
  dbl lam, f[2];

  if (hybrid((hybrid_cost_func_t)utri_hybrid_f, 0, 1, utri, &lam))
    return;

  utri_set_lambda(utri, 0);
  f[0] = utri->f;

  utri_set_lambda(utri, 1);
  f[1] = utri->f;

  assert(f[0] != f[1]);

  if (f[0] < f[1])
    utri_set_lambda(utri, 0);
  else
    utri_set_lambda(utri, 1);
}

static void get_update_inds(utri_s const *utri, size_t l[2]) {
  l[0] = utri->l0;
  l[1] = utri->l1;
}

static void get_bary_coords(utri_s const *utri, dbl b[2]) {
  b[0] = 1 - utri->lam;
  b[1] = utri->lam;
}

par3_s utri_get_par(utri_s const *u) {
  par3_s par = {.l = {[2] = NO_PARENT}, .b = {[2] = NAN}};
  get_update_inds(u, par.l);
  get_bary_coords(u, par.b);
  return par;
}

dbl utri_get_value(utri_s const *u) {
  return u->f;
}

void utri_get_jet(utri_s const *utri, jet3 *jet) {
  jet->f = utri->f;
  jet->fx = utri->x_minus_xb[0]/utri->L;
  jet->fy = utri->x_minus_xb[1]/utri->L;
  jet->fz = utri->x_minus_xb[2]/utri->L;
}

static dbl get_lag_mult(utri_s const *utri) {
  dbl const atol = 1e-15;
  if (utri->lam < atol) {
    return utri->Df;
  } else if (utri->lam > 1 - atol) {
    return -utri->Df;
  } else {
    return 0;
  }
}

static ray3 get_ray(utri_s const *utri) {
  ray3 ray;
  dbl3_saxpy(utri->lam, utri->x1_minus_x0, utri->x0, ray.org);
  dbl3_normalized(utri->x_minus_xb, ray.dir);
  return ray;
}

static dbl get_L(utri_s const *u) {
  return u->L;
}

bool utri_update_ray_is_physical(utri_s const *utri, eik3_s const *eik) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  size_t l[2] = {utri->l0, utri->l1};

  /**
   * First, we check if the ray is spuriously emanating from the
   * boundary. We do this by finding the ray direction, and checking
   * to see if the point just before the start of the ray is inside
   * the mesh.
   */

  // TODO: the following section where we check to see if the stuff
  // below gives "an interior ray" can be wrapped up and reused for
  // both this and the corresponding section in utetra.c...

  // TODO: we can accelerate this a bit by verifying that both l[0]
  // and l[1] are boundary indices...

  ray3 ray = get_ray(utri);

  // Get points just before and just after the start of the ray. We
  // perturb forward and backward by one half of the minimum triangle
  // altitude (taken of the entire mesh). This is to ensure that the
  // perturbations are small enough to stay inside neighboring
  // tetrahedra, but also large enough to take into consideration the
  // characteristic length scale of the mesh.
  dbl xm[3], xp[3], h = mesh3_get_min_edge_length(mesh)/4;
  ray3_get_point(&ray, -h, xm);
  ray3_get_point(&ray, h, xp);

  array_s *cells;
  array_alloc(&cells);
  array_init(cells, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  int nvc;
  size_t *vc;

  for (int i = 0; i < 2; ++i) {
    nvc = mesh3_nvc(mesh, l[i]);
    vc = malloc(nvc*sizeof(size_t));
    mesh3_vc(mesh, l[i], vc);

    for (int j = 0; j < nvc; ++j)
      if (!array_contains(cells, &vc[j]))
        array_append(cells, &vc[j]);

    free(vc);
  }

  size_t lc;
  bool xm_in_cell = false, xp_in_cell = false;
  for (size_t i = 0; i < array_size(cells); ++i) {
    array_get(cells, i, &lc);
    xm_in_cell |= mesh3_cell_contains_point(mesh, lc, xm);
    xp_in_cell |= mesh3_cell_contains_point(mesh, lc, xp);
    if (xm_in_cell && xp_in_cell)
      break;
  }

  // Free the space used for the adjacent cells since we're done with
  // it now
  array_deinit(cells);
  array_dealloc(&cells);

  // If we didn't find a containing cell, we can conclude that the ray
  // is unphysical!
  if (!xm_in_cell || !xp_in_cell)
    return false;

  // TODO: need to check for boundary faces here! See implementation
  // of `utetra_update_ray_is_physical`.

  /**
   * Next, check and see if the point just before the end of the ray
   * lies in a cell.
   */

  dbl L = get_L(utri);

  dbl xhatm[3];
  ray3_get_point(&ray, L - h, xhatm);

  /* If `utri->l` is unset, we're probably computing a jet for a cut
   * point, in which case we can fairly safely assume we're OK
   * here. Return. */
  if (utri->l == (size_t)NO_INDEX)
    return true;

  nvc = mesh3_nvc(mesh, utri->l);
  vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(mesh, utri->l, vc);

  bool xhatm_in_cell = false;
  for (int i = 0; i < nvc; ++i) {
    xhatm_in_cell = mesh3_cell_contains_point(mesh, vc[i], xhatm);
    if (xhatm_in_cell)
      break;
  }

  free(vc);

  return xhatm_in_cell;
}

int utri_cmp(utri_s const **h1, utri_s const **h2) {
  utri_s const *u1 = *h1;
  utri_s const *u2 = *h2;

  if (u1 == NULL && u2 == NULL) {
    return 0;
  } else if (u2 == NULL) {
    return -1;
  } else if (u1 == NULL) {
    return 1;
  } else {
    dbl T1 = utri_get_value(u1), T2 = utri_get_value(u2);
    if (T1 < T2) {
      return -1;
    } else if (T1 > T2) {
      return 1;
    } else {
      return 0;
    }
  }
}

bool utri_has_interior_point_solution(utri_s const *utri) {
  dbl const atol = 1e-14;
  return (atol < utri->lam && utri->lam < 1 - atol)
    || fabs(get_lag_mult(utri)) <= atol;
}

bool utri_has_orig_index(utri_s const *utri) {
  return utri->orig_index != (size_t)NO_INDEX;
}

size_t utri_get_orig_index(utri_s const *utri) {
  return utri->orig_index;
}

bool utri_is_finite(utri_s const *u) {
  return isfinite(u->f);
}

size_t utri_get_active_ind(utri_s const *utri) {
  dbl const atol = 1e-15;
  if (utri->lam < atol)
    return utri->l0;
  if (utri->lam > 1 - atol)
    return utri->l1;
  return (size_t)NO_INDEX;
}

size_t utri_get_inactive_ind(utri_s const *utri) {
  dbl const atol = 1e-15;
  if (utri->lam < atol)
    return utri->l1;
  if (utri->lam > 1 - atol)
    return utri->l0;
  return (size_t)NO_INDEX;
}

bool utri_contains_update_ind(utri_s const *u, size_t l) {
  return l == u->l0 || l == u->l1;
}

size_t utri_get_l(utri_s const *utri) {
  return utri->l;
}

bool utri_opt_inc_on_other_utri(utri_s const *u, utri_s const *other_u) {
  dbl const atol = 1e-14;
  dbl xlam[3];
  dbl3_saxpy(u->lam, u->x1_minus_x0, u->x0, xlam);
  dbl dist0 = dbl3_dist(other_u->x0, xlam);
  dbl dist1 = dbl3_dist(other_u->x1, xlam);
  return dist0 < atol || dist1 < atol;
}

static dbl get_lambda(utri_s const *u) {
  return u->lam;
}

bool utris_yield_same_update(utri_s const *utri1, utri_s const *utri2) {
  dbl const atol = 1e-14;

  if (utri1 == NULL || utri2 == NULL)
    return false;

  dbl lam1 = get_lambda(utri1);
  dbl lam2 = get_lambda(utri2);

  if (fabs(lam1 - lam2) > atol)
    return false;

  jet3 jet1, jet2;
  utri_get_jet(utri1, &jet1);
  utri_get_jet(utri2, &jet2);

  return jet3_approx_eq(&jet1, &jet2, atol);
}

#if JMM_TEST
bool utri_is_causal(utri_s const *utri) {
  dbl dx0[3], dx1[3];
  dbl3_sub(utri->x0, utri->x, dx0);
  dbl3_sub(utri->x1, utri->x, dx1);

  dbl cos01 = dbl3_dot(dx0, dx1)/(dbl3_norm(dx0)*dbl3_norm(dx1));

  return cos01 >= 0;
}

dbl utri_get_lambda(utri_s const *utri) {
  return get_lambda(utri);
}
#endif
