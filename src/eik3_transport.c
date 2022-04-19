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
