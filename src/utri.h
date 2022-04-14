#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "jet.h"
#include "par.h"

typedef struct utri_spec {
  eik3_s const *eik;
  size_t lhat, l[2];
  state_e state[2];
  dbl xhat[3], x[2][3];
  jet31t jet[2];
  size_t orig_index;
} utri_spec_s;

utri_spec_s utri_spec_empty();
utri_spec_s utri_spec_from_eik(eik3_s const *eik, size_t l, size_t l0, size_t l1);
utri_spec_s utri_spec_from_eik_without_l(eik3_s const *eik, dbl const x[3],
                                         size_t l0, size_t l1);
utri_spec_s utri_spec_from_raw_data(dbl3 const x, dbl3 const Xt[2], jet31t const jet[2]);

typedef struct utri utri_s;

void utri_alloc(utri_s **utri);
void utri_dealloc(utri_s **utri);
void utri_init(utri_s *u, utri_spec_s const *spec);
void utri_solve(utri_s *utri);
par3_s utri_get_par(utri_s const *u);
dbl utri_get_value(utri_s const *utri);
void utri_get_jet31t(utri_s const *utri, jet31t *jet);
void utri_get_jet32t(utri_s const *utri, jet32t *jet);
bool utri_emits_terminal_ray(utri_s const *utri, eik3_s const *eik);
bool utri_update_ray_is_physical(utri_s const *utri, eik3_s const *eik);
int utri_cmp(utri_s const **h1, utri_s const **h2);
bool utri_has_interior_point_solution(utri_s const *utri);
bool utri_has_orig_index(utri_s const *utri);
size_t utri_get_orig_index(utri_s const *utri);
bool utri_is_finite(utri_s const *utri);
size_t utri_get_active_ind(utri_s const *utri);
size_t utri_get_inactive_ind(utri_s const *utri);
bool utri_contains_update_ind(utri_s const *utri, size_t l);
size_t utri_get_l(utri_s const *utri);
bool utri_opt_inc_on_other_utri(utri_s const *u, utri_s const *other_u);
void utri_get_update_inds(utri_s const *u, size_t l[2]);
void utri_get_t(utri_s const *u, dbl t[3]);
dbl utri_get_L(utri_s const *u);
dbl utri_get_b(utri_s const *u);
bool utri_inc_on_diff_edge(utri_s const *u, mesh3_s const *mesh);
void utri_get_xb(utri_s const *u, dbl xb[3]);
bool utri_is_degenerate(utri_s const *u);
bool utri_approx_hess(utri_s const *u, dbl h, dbl33 hess);
bool utri_inc_on_refl_BCs(utri_s const *u, eik3_s const *eik);
bool utri_accept_refl_BCs_update(utri_s const *u, eik3_s const *eik);

bool utris_yield_same_update(utri_s const *utri1, utri_s const *utri2);

#if JMM_TEST
bool utri_is_causal(utri_s const *utri);
dbl utri_get_lambda(utri_s const *utri);
#endif

#ifdef __cplusplus
}
#endif
