#include "eik3_transport.h"

#include <assert.h>
#include <math.h>

#include "mesh3.h"

static void transport_dbl(eik3_s const *eik, size_t l0, dbl *values) {
  par3_s par = eik3_get_par(eik, l0);

  if (par3_is_empty(&par))
    return;

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  values[l0] = 0;
  for (size_t i = 0; i < num_active; ++i) {
    assert(isfinite(values[l[i]]));
    values[l0] += b[i]*values[l[i]];
  }
}

void eik3_transport_dbl(eik3_s const *eik, dbl *values, bool skip_filled) {
  assert(eik3_is_solved(eik));

  mesh3_s const *mesh = eik3_get_mesh(eik);

  size_t nverts = mesh3_nverts(mesh);

  size_t const *accepted = eik3_get_accepted_ptr(eik);

  for (size_t i = 0; i < nverts; ++i) {
    size_t l0 = accepted[i];
    if (skip_filled && !isnan(values[l0]))
      continue;
    transport_dbl(eik, l0, values);
  }
}

static void transport_dblz(eik3_s const *eik, size_t l0, dblz *values) {
  par3_s par = eik3_get_par(eik, l0);

  if (par3_is_empty(&par))
    return;

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  values[l0] = 0;
  for (size_t i = 0; i < num_active; ++i) {
    assert(isfinite(creal(values[l[i]])));
    assert(isfinite(cimag(values[l[i]])));
    values[l0] += b[i]*values[l[i]];
  }
}

void eik3_transport_dblz(eik3_s const *eik, dblz *values, bool skip_filled) {
  assert(eik3_is_solved(eik));

  mesh3_s const *mesh = eik3_get_mesh(eik);

  size_t const *accepted = eik3_get_accepted_ptr(eik);

  size_t nverts = mesh3_nverts(mesh);
  for (size_t i = 0; i < nverts; ++i) {
    size_t l0 = accepted[i];
    dblz z = values[l0];
    if (skip_filled && !isnan(creal(z)) && !isnan(cimag(z)))
      continue;
    transport_dblz(eik, l0, values);
  }
}

static void transport_curvature(eik3_s const *eik, size_t l0, dbl *kappa) {
  if (isfinite(kappa[l0]))
    return;

  par3_s par = eik3_get_par(eik, l0);

  if (par3_is_empty(&par))
    return;

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  bool finite[num_active];

  size_t num_finite = 0;
  for (size_t i = 0; i < num_active; ++i) {
    num_finite += finite[i] = isfinite(kappa[l[i]]);
    if (!finite[i])
      assert(isinf(kappa[l[i]]));
  }

  if (num_finite < num_active) {
    assert(num_finite == 1);
    assert(num_active == 2);

    dbl rho = 0;
    for (size_t i = 0; i < num_active; ++i)
      rho += b[i]/kappa[l[i]];
    kappa[l0] = 1/rho;

    return;
  }

  kappa[l0] = 0;
  for (size_t i = 0; i < num_active; ++i) {
    assert(isfinite(kappa[l[i]]));
    kappa[l0] += b[i]*kappa[l[i]];
  }
}

void eik3_transport_curvature(eik3_s const *eik, dbl *kappa, bool skip_filled) {
  assert(eik3_is_solved(eik));

  mesh3_s const *mesh = eik3_get_mesh(eik);
  size_t nverts = mesh3_nverts(mesh);

  size_t const *accepted = eik3_get_accepted_ptr(eik);

  for (size_t i = 0; i < nverts; ++i) {
    size_t l0 = accepted[i];
    if (skip_filled && !isnan(kappa[l0]))
      continue;
    transport_curvature(eik, l0, kappa);
  }
}

static void slerp2(dbl const b[2], dbl3 const p[2], dbl3 q) {
  dbl theta = acos(dbl3_dot(p[0], p[1]));
  dbl sin_theta = sin(theta);
  dbl c[2] = {sin(b[0]*theta)/sin_theta, sin(b[1]*theta)/sin_theta};
  for (size_t i = 0; i < 3; ++i)
    q[i] = c[0]*p[0][i] + c[1]*p[1][i];
}

static void slerp3(dbl const b[3], dbl3 const p[3], dbl3 q) {
  /* initialize q to nlerp */
  for (size_t i = 0; i < 3; ++i)
    q[i] = b[0]*p[0][i] + b[1]*p[1][i] + b[2]*p[2][i];
  dbl3_normalize(q);

  dbl3 q0;
  dbl3_copy(q, q0);

  /* projected gradient descent for slerp */
  while (true) {
    /* take projected gradient step */
    dbl c[3];
    for (size_t i = 0; i < 3; ++i) {
      dbl dot = dbl3_dot(p[i], q);
      c[i] = acos(dot)*b[i]/sqrt(1 - dot*dot);
    }
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        q[j] += c[i]*p[i][j];
    dbl3_normalize(q);

    /* check for convergence */
    if (acos(dbl3_dot(q, q0)) < 1e-10)
      break;

    /* prepare for the next iteration */
    dbl3_copy(q, q0);
  }
}

static void transport_unit_vector(eik3_s const *eik, size_t l0, dbl3 *t) {
  par3_s par = eik3_get_par(eik, l0);
  if (par3_is_empty(&par))
    return;

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  dbl3 t_[3];
  for (size_t i = 0; i < num_active; ++i)
    dbl3_copy(t[l[i]], t_[i]);

  if (num_active == 1)
    dbl3_copy(t_[0], t[l0]);
  else if (num_active == 2)
    slerp2(b, t_, t[l0]);
  else if (num_active == 3)
    slerp3(b, t_, t[l0]);
  else
    assert(false);
}

void eik3_transport_unit_vector(eik3_s const *eik, dbl3 *t, bool skip_filled) {
  assert(eik3_is_solved(eik));

  mesh3_s const *mesh = eik3_get_mesh(eik);
  size_t nverts = mesh3_nverts(mesh);

  size_t const *accepted = eik3_get_accepted_ptr(eik);

  for (size_t i = 0; i < nverts; ++i) {
    size_t l0 = accepted[i];
    if (skip_filled && dbl3_isfinite(t[l0]))
      continue;
    transport_unit_vector(eik, l0, t);
  }
}
