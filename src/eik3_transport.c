#include <jmm/eik3_transport.h>

#include <assert.h>
#include <math.h>

#include <jmm/mesh3.h>

#include "log.h"
#include "util.h"

static void transport_dbl(eik3_s const *eik, size_t l0, dbl *values) {
  par3_s par = eik3_get_par(eik, l0);

  if (par3_is_empty(&par))
    return;

  size_t num_active = par3_num_active(&par);
  uint3 l = {(size_t)NO_INDEX, (size_t)NO_INDEX, (size_t)NO_INDEX};
  dbl3 b = {NAN, NAN, NAN};
  par3_get_active(&par, l, b);

  dbl3 par_values;
  dbl3_nan(par_values);
  dbl3_gather(values, l, par_values);

  values[l0] = 0;
  for (size_t i = 0; i < num_active; ++i) {
    assert(isfinite(par_values[i]));
    values[l0] += b[i]*par_values[i];
  }

  /* We want to be able to assert that values[l0] is contained in the
   * closed interval determined by par_values afterwards. Floating
   * point roundoff can mess this up. So we fix this here. */
  dbl nanmin = dbl3_nanmin(par_values);
  dbl nanmax = dbl3_nanmax(par_values);
  values[l0] = clamp(values[l0], nanmin, nanmax);
}

void eik3_transport_dbl(eik3_s const *eik, dbl *values, bool skip_filled) {
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

static void nlerp3(dbl const b[3], dbl3 const p[3], dbl3 q) {
  for (size_t i = 0; i < 3; ++i)
    q[i] = b[0]*p[0][i] + b[1]*p[1][i] + b[2]*p[2][i];
  dbl3_normalize(q);
}

static void slerp3(dbl const b[3], dbl3 const p[3], dbl3 q) {
  nlerp3(b, p, q);

  dbl3 q0;
  dbl3_copy(q, q0);

  /* projected gradient descent for slerp */
  size_t it = 0, max_it = 100;
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
    dbl dot = clamp(dbl3_dot(q, q0), -1, 1);
    if (fabs(1 - dot) < 1e-10)
      break;

    /* prepare for the next iteration */
    dbl3_copy(q, q0);

    ++it;

    /* Failed to converge: reset to nlerp and bail */
    if (it == max_it) {
      log_warn("slerp3 failed to converge");
      nlerp3(b, p, q);
      break;
    }
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
