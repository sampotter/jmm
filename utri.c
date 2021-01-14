#include "utri.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

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

  dbl Tc[4];

  dbl cos01;

  dbl lam;

  dbl f;
  dbl Df;

  dbl x_minus_xb[3];
  dbl L;
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
  dbl T = bb3(utri->Tc, b);

  utri->f = T + utri->L;

  dbl a[2] = {-1, 1};
  dbl dT_dlam = dbb3(utri->Tc, b, a);

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
  bb3_interp3(f, Df, Xt, utri->Tc);
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

bool utri_update_ray_is_physical(utri_s const *utri, eik3_s const *eik,
                                 size_t const l[2]) {
  dbl const atol = 1e-14;

  jet3 *jet = eik3_get_jet_ptr(eik);

  dbl *t0 = &jet[l[0]].fx, *t1 = &jet[l[1]].fx;

  // TODO: decide how to handle this case...
  assert(fabs(dbl3_dot(t0, t1)) > atol);

  // First, get plane spanned by jets at base of update.
  dbl n[3];
  dbl3_cross(t0, t1, n); // already normalized

  // If the ray is in the plane spanned by the two jets, then it
  // hasn't diffracted
  //
  // TODO: could we upgrade this requirement to "ray is in the cone"?
  if (fabs(dbl3_dot(n, utri->x_minus_xb)) <= utri->L*atol)
    return true;

  // If the ray *isn't* in that plane, and the base of the update
  // isn't a diffracting edge... we have a problem!
  mesh3_s const *mesh = eik3_get_mesh(eik);
  return mesh3_is_diff_edge(mesh, l);
}
