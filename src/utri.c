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

struct utri {
  dbl x[3];
  dbl x0[3];
  dbl x1[3];
  dbl x1_minus_x0[3];

  bb31 T;

  dbl cos01;

  dbl lam;

  dbl f;
  dbl Df;

  dbl x_minus_xb[3];
  dbl L;

  size_t l, l0, l1;

  utri_s *split;

  state_e state[2];
  size_t num_shadow;

  int i; // original index
};

static bool is_split(utri_s const *u) {
  return u->split != NULL;
}

void utri_alloc(utri_s **utri) {
  *utri = malloc(sizeof(utri_s));
}

void utri_dealloc(utri_s **utri) {
  free(*utri);
  *utri = NULL;
}

void utri_set_lambda(utri_s *utri, dbl lam) {
  assert(!is_split(utri));

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

static void init_split_utri(utri_s *u, eik3_s const *eik) {
  u->split = NULL;

  if (u->num_shadow != 1)
    return;

  dbl t;
  assert(eik3_get_cutedge_t(eik, u->l0, u->l1, &t));

  dbl const atol = 1e-15;

  if (t < atol) {
    if (u->state[0] == VALID) {
      u->state[0] = SHADOW;
      u->num_shadow = 2;
    } else {
      u->state[0] = VALID;
      u->num_shadow = 0;
    }
    return;
  }

  if (t > 1 - atol) {
    if (u->state[1] == VALID) {
      u->state[1] = SHADOW;
      u->num_shadow = 2;
    } else {
      u->state[1] = VALID;
      u->num_shadow = 0;
    }
    return;
  }

  u->split = malloc(sizeof(utri_s *));
  utri_alloc(&u->split);

  dbl Xt[2][3];
  dbl3_copy(u->x0, Xt[0]);
  dbl3_saxpy(t, u->x1_minus_x0, u->x0, Xt[1]);

  jet3 jet[2];
  jet[0] = eik3_get_jet(eik, u->l0);
  assert(eik3_get_cutedge_jet(eik, u->l0, u->l1, &jet[1]));

  utri_init(u->split, u->x, Xt, jet);

  u->split->l = NO_INDEX;
  u->split->l0 = NO_INDEX;
  u->split->l1 = NO_INDEX;

  state_e *state = u->split->state;
  state[0] = state[1] = VALID;
  u->split->num_shadow = 0;
}

void utri_init_from_eik3(utri_s *utri, eik3_s const *eik, size_t l,
                         size_t l0, size_t l1) {
  dbl x[3];
  mesh3_copy_vert(eik3_get_mesh(eik), l, x);

  utri->l = l;

  utri_init_from_eik3_without_l(utri, eik, x, l0, l1);
}

void utri_init_from_eik3_without_l(utri_s *utri, eik3_s const *eik,
                                   dbl const x[3], size_t l0, size_t l1) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  dbl Xt[2][3];
  mesh3_copy_vert(mesh, l0, Xt[0]);
  mesh3_copy_vert(mesh, l1, Xt[1]);

  jet3 jet[2] = {
    eik3_get_jet(eik, l0),
    eik3_get_jet(eik, l1)
  };

  assert(jet3_is_finite(&jet[0]));
  assert(jet3_is_finite(&jet[1]));

  utri_init(utri, x, Xt, jet);

  utri->l0 = l0;
  utri->l1 = l1;

  state_e const *state = eik3_get_state_ptr(eik);

  assert(state[l0] == VALID || state[l0] == SHADOW);
  assert(state[l1] == VALID || state[l1] == SHADOW);

  utri->state[0] = state[l0];
  utri->state[1] = state[l1];

  utri->num_shadow = (utri->state[0] == SHADOW) + (utri->state[1] == SHADOW);

  init_split_utri(utri, eik);

  assert(0 <= utri->num_shadow && utri->num_shadow <= 2);
}

void utri_init(utri_s *utri, dbl const x[3], dbl const Xt[2][3],
               jet3 const jet[2]) {
  utri->f = INFINITY;
  utri->i = NO_INDEX;

  memcpy(utri->x, x, 3*sizeof(dbl));
  memcpy(utri->x0, Xt[0], 3*sizeof(dbl));
  memcpy(utri->x1, Xt[1], 3*sizeof(dbl));

  dbl3_sub(utri->x1, utri->x0, utri->x1_minus_x0);

  dbl dx0[3], dx1[3];
  dbl3_sub(utri->x0, utri->x, dx0);
  dbl3_sub(utri->x1, utri->x, dx1);
  utri->cos01 = dbl3_dot(dx0, dx1)/(dbl3_norm(dx0)*dbl3_norm(dx1));

  dbl f[2] = {jet[0].f, jet[1].f};
  dbl Df[2][3] = {
    {jet[0].fx, jet[0].fy, jet[0].fz},
    {jet[1].fx, jet[1].fy, jet[1].fz}
  };
  bb31_init_from_3d_data(&utri->T, f, Df, Xt);
}

void utri_solve(utri_s *utri) {
  if (is_split(utri))
    utri_solve(utri->split);
  else
    (void)hybrid((hybrid_cost_func_t)utri_hybrid_f, 0, 1, utri);
}

static void get_update_inds(utri_s const *utri, size_t l[2]) {
  l[0] = utri->l0;
  l[1] = utri->l1;
}

static void get_bary_coords(utri_s const *utri, dbl b[2]) {
  assert(!is_split(utri));
  b[0] = 1 - utri->lam;
  b[1] = utri->lam;
}

static void get_x(utri_s const *u, dbl x[3]) {
  if (is_split(u))
    return get_x(u->split, x);
  dbl3_saxpy(u->lam, u->x1_minus_x0, u->x0, x);
}

par3_s utri_get_par(utri_s const *u) {
  par3_s par = {.l = {[2] = NO_PARENT}, .b = {[2] = NAN}};

  get_update_inds(u, par.l);

  if (is_split(u)) {
    dbl xt[3]; get_x(u, xt);
    par.b[1] = dbl3_dot(u->x1_minus_x0, xt);
    par.b[0] = 1 - par.b[1];
  } else {
    get_bary_coords(u, par.b);
  }

  return par;
}

dbl utri_get_value(utri_s const *u) {
  return is_split(u) ? utri_get_value(u->split) : u->f;
}

void utri_get_jet(utri_s const *utri, jet3 *jet) {
  if (is_split(utri))
    return utri_get_jet(utri->split, jet);

  jet->f = utri->f;
  jet->fx = utri->x_minus_xb[0]/utri->L;
  jet->fy = utri->x_minus_xb[1]/utri->L;
  jet->fz = utri->x_minus_xb[2]/utri->L;
}

static dbl get_lag_mult(utri_s const *utri) {
  assert(!is_split(utri));
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
  if (is_split(utri))
    return get_ray(utri->split);

  ray3 ray;
  dbl3_saxpy(utri->lam, utri->x1_minus_x0, utri->x0, ray.org);
  dbl3_normalized(utri->x_minus_xb, ray.dir);
  return ray;
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
  dbl xm[3], xp[3], t = mesh3_get_min_tetra_alt(mesh)/2;
  ray3_get_point(&ray, -t/2, xm);
  ray3_get_point(&ray, t/2, xp);

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
    xm_in_cell |= mesh3_dbl3_in_cell(mesh, lc, xm, NULL);
    xp_in_cell |= mesh3_dbl3_in_cell(mesh, lc, xp, NULL);
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

  dbl xhatm[3];
  ray3_get_point(&ray, utri->L - t, xhatm);

  nvc = mesh3_nvc(mesh, utri->l);
  vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(mesh, utri->l, vc);

  bool xhatm_in_cell = false;
  for (int i = 0; i < nvc; ++i) {
    xhatm_in_cell = mesh3_dbl3_in_cell(mesh, vc[i], xhatm, NULL);
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
  if (is_split(utri))
    return utri_has_interior_point_solution(utri->split);
  dbl const atol = 1e-14;
  return (atol < utri->lam && utri->lam < 1 - atol)
    || fabs(get_lag_mult(utri)) <= atol;
}

void utri_set_orig_index(utri_s *utri, int i) {
  assert(!is_split(utri));

  utri->i = i;
}

int utri_get_orig_index(utri_s const *utri) {
  assert(!is_split(utri));

  return utri->i;
}

bool utri_is_finite(utri_s const *u) {
  return isfinite(is_split(u) ? u->split->f : u->f);
}

static dbl get_lambda(utri_s const *u) {
  return is_split(u) ? u->split->lam : u->lam;
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
  assert(!is_split(utri));

  return utri->cos01 >= 0;
}

dbl utri_get_lambda(utri_s const *utri) {
  return get_lambda(utri);
}
#endif
