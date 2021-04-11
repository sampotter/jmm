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

  int i; // original index
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

void utri_init_from_eik3(utri_s *utri, eik3_s const *eik, size_t l,
                         size_t l0, size_t l1) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  dbl x[3];
  mesh3_copy_vert(mesh, l, x);

  dbl Xt[2][3];
  mesh3_copy_vert(mesh, l0, Xt[0]);
  mesh3_copy_vert(mesh, l1, Xt[1]);

  jet3 jet[2] = {
    eik3_get_jet(eik, l0),
    eik3_get_jet(eik, l1)
  };

  assert(jet3_is_finite(&jet[0]));
  assert(jet3_is_finite(&jet[1]));

  utri->l = l;
  utri->l0 = l0;
  utri->l1 = l1;

  utri_init(utri, x, Xt, jet);
}

void utri_init(utri_s *utri, dbl const x[3], dbl const Xt[2][3],
               jet3 const jet[2]) {
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

bool utri_is_causal(utri_s const *utri) {
  return utri->cos01 >= 0;
}

void utri_solve(utri_s *utri) {
  (void)hybrid((hybrid_cost_func_t)utri_hybrid_f, 0, 1, utri);
}

dbl utri_get_lambda(utri_s const *utri) {
  return utri->lam;
}

void utri_get_bary_coords(utri_s const *utri, dbl b[2]) {
  b[0] = 1 - utri->lam;
  b[1] = utri->lam;
}

dbl utri_get_value(utri_s const *utri) {
  return utri->f;
}

void utri_get_jet(utri_s const *utri, jet3 *jet) {
  jet->f = utri->f;
  jet->fx = utri->x_minus_xb[0]/utri->L;
  jet->fy = utri->x_minus_xb[1]/utri->L;
  jet->fz = utri->x_minus_xb[2]/utri->L;
}

dbl utri_get_lag_mult(utri_s const *utri) {
  dbl const atol = 1e-15;
  if (utri->lam < atol) {
    return utri->Df;
  } else if (utri->lam > 1 - atol) {
    return -utri->Df;
  } else {
    return 0;
  }
}

void utri_get_point_on_ray(utri_s const *utri, dbl t, dbl xt[3]) {
  dbl xlam[3];
  dbl3_saxpy(utri->lam, utri->x1_minus_x0, utri->x0, xlam);
  dbl3_saxpy(t/utri->L, utri->x_minus_xb, xlam, xt);
}

bool utri_update_ray_is_physical(utri_s const *utri, eik3_s const *eik) {
  dbl const atol = 1e-14;

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

  // Get points just before and just after the start of the ray. We
  // perturb forward and backward by one half of the minimum triangle
  // altitude (taken of the entire mesh). This is to ensure that the
  // perturbations are small enough to stay inside neighboring
  // tetrahedra, but also large enough to take into consideration the
  // characteristic length scale of the mesh.
  dbl xm[3], xp[3], t = mesh3_get_min_tetra_alt(mesh)/2;
  utri_get_point_on_ray(utri, -t/2, xm);
  utri_get_point_on_ray(utri, t/2, xp);

  array_s *cells;
  array_alloc(&cells);
  array_init(cells, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  bool I[2] = {fabs(1 - utri->lam) > atol, fabs(utri->lam) > atol};

  int nvc;
  size_t *vc;

  for (int i = 0; i < 2; ++i) {
    if (!I[i]) continue;

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

  /**
   * Next, check and see if the point just before the end of the ray
   * lies in a cell.
   */

  dbl xhatm[3];
  dbl3_saxpy(-t/utri->L, utri->x_minus_xb, utri->x, xhatm);

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

void utri_reset(utri_s *utri) {
  utri->f = INFINITY;
  utri->i = NO_INDEX;
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

bool utri_has_interior_point_solution(utri_s const *utri) {
  dbl const atol = 1e-14;
  return (atol < utri->lam && utri->lam < 1 - atol)
    || fabs(utri_get_lag_mult(utri)) <= atol;
}

void utri_get_update_inds(utri_s const *utri, size_t l[2]) {
  l[0] = utri->l0;
  l[1] = utri->l1;
}

void utri_set_orig_index(utri_s *utri, int i) {
  utri->i = i;
}

int utri_get_orig_index(utri_s const *utri) {
  return utri->i;
}

bool utri_is_finite(utri_s const *utri) {
  return isfinite(utri->f);
}

bool utris_yield_same_update(utri_s const *utri1, utri_s const *utri2) {
  dbl const atol = 1e-14;

  jet3 jet1, jet2;
  utri_get_jet(utri1, &jet1);
  utri_get_jet(utri2, &jet2);

  return fabs(utri1->lam - utri2->lam) <= atol
    && jet3_approx_eq(&jet1, &jet2, atol);
}
