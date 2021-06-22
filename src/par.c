#include "par.h"

#include <assert.h>
#include <math.h>

#include "array.h"
#include "eik3.h"
#include "macros.h"
#include "mesh3.h"

par3_s make_par3(size_t l[3], dbl b[3]) {
  par3_s par = {.l = {l[0], l[1], l[2]}, .b = {b[0], b[1], b[2]}};
  for (int i = 0; i < 2; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      if (par.b[i] < par.b[j]) {
        SWAP(par.l[i], par.l[j]);
        SWAP(par.b[i], par.b[j]);
      }
    }
  }
  dbl const atol = 1e-15;
  assert(par.b[0] > atol);
  if (par.b[1] <= atol) {
    par.l[1] = NO_PARENT;
    par.b[1] = NAN;
  }
  if (par.b[2] <= atol) {
    par.l[2] = NO_PARENT;
    par.b[2] = NAN;
  }
  return par;
}

void par3_init_empty(par3_s *par) {
  par->l[0] = par->l[1] = par->l[2] = NO_PARENT;
  par->b[0] = par->b[1] = par->b[2] = NAN;
}

void par3_set(par3_s *par, size_t const *l, dbl const *b, int n) {
  for (int i = 0; i < n; ++i) {
    par->l[i] = l[i];
    par->b[i] = b[i];
  }
}

size_t par3_size(par3_s const *par) {
  size_t size = (int)(par->l[0] != NO_PARENT) + (int)(par->l[1] != NO_PARENT)
    + (int)(par->l[2] != NO_PARENT);
  assert(0 <= size && size <= 3);
  return size;
}

void par3_get_xb(par3_s const *par, eik3_s const *eik, dbl xb[3]) {
  dbl const *x[3];
  mesh3_get_vert_ptrs(eik3_get_mesh(eik), par->l, 3, x);

  for (int i = 0; i < 3; ++i) {
    xb[i] = 0;
    for (int j = 0; j < 3; ++j) {
      xb[i] += par->b[j]*x[j][i];
    }
  }
}

bool par3_is_empty(par3_s const *par) {
  return par->l[0] == NO_PARENT;
}

size_t par3_num_active(par3_s const *par) {
  dbl const atol = 1e-14;
  size_t num_active = 0;
  for (size_t i = 0; i < 3; ++i)
    num_active += par->l[i] != NO_PARENT && par->b[i] > atol;
  return num_active;
}

void par3_get_active(par3_s const *par, size_t *l, dbl *b) {
  dbl const atol = 1e-14;
  size_t j = 0;
  for (size_t i = 0; i < 3; ++i)
    if (par->l[i] != NO_PARENT && par->b[i] > atol) {
      l[j] = par->l[i];
      b[j++] = par->b[i];
    }

  /* Normalize the active coefficients */
  dbl b_norm1 = 0;
  for (size_t i = 0; i < j; ++i)
    b_norm1 += fabs(b[i]);
  for (size_t i = 0; i < j; ++i)
    b[i] /= b_norm1;
}

bool par3_is_on_BC_boundary(par3_s const *par, eik3_s const *eik) {
  ftype_e ftype = eik3_get_ftype(eik);
  if (ftype == FTYPE_POINT_SOURCE)
    return false;

  size_t num_active = par3_num_active(par);
  size_t l[2];
  dbl b[num_active];
  par3_get_active(par, l, b);

  mesh3_s const *mesh = eik3_get_mesh(eik);

  for (size_t i = 0; i < num_active; ++i)
    if (!mesh3_bdv(mesh, l[i]))
      return false;

  if (ftype == FTYPE_REFLECTION && num_active == 1) {
    size_t num_inc_bdf = mesh3_get_num_inc_bdf(mesh, l[0], false);
    size_t (*lf)[3] = malloc(num_inc_bdf*sizeof(size_t[3]));
    mesh3_get_inc_bdf(mesh, l[0], lf, false);

    size_t num_inc_bdf_w_BCs = 0;
    for (size_t i = 0; i < num_inc_bdf; ++i)
      num_inc_bdf_w_BCs += eik3_has_BCs(eik, lf[i][0])
        && eik3_has_BCs(eik, lf[i][1]) && eik3_has_BCs(eik, lf[i][2]);

    free(lf);

    assert(num_inc_bdf > 0);
    assert(num_inc_bdf_w_BCs <= num_inc_bdf);

    // TODO: this may not be totally correct... easy to imagine some
    // corner cases where this may fail
    return num_inc_bdf_w_BCs < num_inc_bdf;
  }

  if (ftype == FTYPE_REFLECTION && num_active == 2) {
    if (!mesh3_bde(mesh, l))
      return false;

    size_t nlf = mesh3_get_num_inc_bdf(mesh, l[0], false);
    size_t (*lf)[3] = malloc(nlf*sizeof(size_t[3]));
    mesh3_get_inc_bdf(mesh, l[0], lf, false);

    size_t num_inc_bdf = 0;
    size_t num_inc_bdf_w_BCs = 0;

    for (size_t i = 0; i < nlf; ++i) {
      if (!edge_in_face(l, lf[i]))
        continue;

      ++num_inc_bdf;

      num_inc_bdf_w_BCs += eik3_has_BCs(eik, lf[i][0])
        && eik3_has_BCs(eik, lf[i][1]) && eik3_has_BCs(eik, lf[i][2]);
    }

    free(lf);

    assert(num_inc_bdf > 0);
    assert(num_inc_bdf <= 2);
    assert(num_inc_bdf_w_BCs <= num_inc_bdf);

    // TODO: this may not be totally correct... easy to imagine some
    // corner cases where this may fail
    return num_inc_bdf_w_BCs < num_inc_bdf;
  }

  if (ftype == FTYPE_REFLECTION && num_active == 3)
    return false;

  if (ftype == FTYPE_EDGE_DIFFRACTION) {
    if (num_active == 2)
      return false;

    assert(num_active == 1);

    size_t num_inc_diff_edges = mesh3_get_num_inc_diff_edges(mesh, l[0]);
    if (num_inc_diff_edges <= 1)
      return false;

    size_t nvv = mesh3_nvv(mesh, l[0]);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(mesh, l[0], vv);

    size_t num_inc_diff_edges_w_BCs = 0;
    for (size_t i = 0; i < nvv; ++i) {
      l[1] = vv[i];
      num_inc_diff_edges_w_BCs += eik3_has_bde_bc(eik, l);
    }

    free(vv);

    return num_inc_diff_edges_w_BCs == 1;
  }

  assert(false);
}
