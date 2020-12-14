#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "jet.h"

typedef struct mesh3 mesh3_s;

typedef struct costfunc {
  dbl f;
  dbl g[2];
  dbl H[2][2];

  dbl p[2]; // Newton step

  dbl x[3]; // x[l]
  dbl X[3][3]; // X = [x[l0]'; x[l1]'; x[l2]']
  dbl Xt[3][3];
  dbl XXt[3][3];

  // B-coefs for 9-point triangle interpolation T on base of update
  dbl Tc[10];

  dbl x_minus_xb[3];

  int niter;
} costfunc_s;

void costfunc_init(costfunc_s *cf, mesh3_s const *mesh, jet3 const *jet,
                   size_t l, size_t l0, size_t l1, size_t l2);
void costfunc_set_lambda(costfunc_s *cf, dbl const *lambda);
void tetra(costfunc_s *cf, dbl const *lambda, jet3 *jet);

/**
 * A structure managing a jet marching method solving the eikonal
 * equation in 3D on an unstructured tetrahedron mesh.
 *
 * NOTE: this is just for s = 1 at the moment. Will extend this to
 * handle s != later.
 */
typedef struct eik3 eik3_s;

void eik3_alloc(eik3_s **eik);
void eik3_dealloc(eik3_s **eik);
void eik3_init(eik3_s *eik, mesh3_s const *mesh);
void eik3_deinit(eik3_s *eik);
void eik3_step(eik3_s *eik);
void eik3_solve(eik3_s *eik);
void eik3_add_trial(eik3_s *eik, size_t ind, jet3 jet);
void eik3_add_valid(eik3_s *eik, size_t ind, jet3 jet);
bool eik3_is_valid(eik3_s const *eik, size_t ind);
jet3 *eik3_get_jet_ptr(eik3_s const *eik);

#ifdef __cplusplus
}
#endif
