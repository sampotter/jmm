#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "geom.h"
#include "jet.h"
#include "par.h"

typedef struct utetra utetra_s;

void utetra_alloc(utetra_s **cf);
void utetra_dealloc(utetra_s **cf);
void utetra_reset(utetra_s *cf);
bool utetra_init_from_eik3(utetra_s *cf, eik3_s const *eik,
                           size_t l, size_t l0, size_t l1, size_t l2);
bool utetra_init_from_ptrs(utetra_s *cf, mesh3_s const *mesh, jet3 const *jet,
                           size_t l, size_t l0, size_t l1, size_t l2);
bool utetra_init_no_inds(utetra_s *cf, dbl const x[3], dbl const Xt[3][3],
                         jet3 const jet[3]);
bool utetra_init(utetra_s *cf, dbl const x[3], dbl const Xt[3][3],
                 jet3 const jet[3]);
void utetra_deinit(utetra_s *u);
bool utetra_is_degenerate(utetra_s const *cf);
void utetra_step(utetra_s *cf);
void utetra_solve(utetra_s *cf);
void utetra_get_lambda(utetra_s const *cf, dbl lam[2]);
void utetra_set_lambda(utetra_s *cf, dbl const lam[2]);
void utetra_get_bary_coords(utetra_s const *cf, dbl b[3]);
dbl utetra_get_value(utetra_s const *cf);
void utetra_get_gradient(utetra_s const *cf, dbl g[2]);
void utetra_get_jet(utetra_s const *cf, jet3 *jet);
void utetra_get_lag_mults(utetra_s const *cf, dbl alpha[3]);
int utetra_get_num_iter(utetra_s const *cf);
bool utetra_has_interior_point_solution(utetra_s const *cf);
int utetra_cmp(utetra_s const **h1, utetra_s const **h2);
bool utetra_adj_are_optimal(utetra_s const *u1, utetra_s const *u2);
bool utetra_diff(utetra_s const *utetra, mesh3_s const *mesh, size_t const l[3]);
ray3 utetra_get_ray(utetra_s const *utetra);
void utetra_get_point_on_ray(utetra_s const *utetra, dbl t, dbl xt[3]);
int utetra_get_interior_coefs_mask(utetra_s const *utetra, bool I[3]);
bool utetra_inds_are_set(utetra_s const *utetra);
bool utetra_update_ray_is_physical(utetra_s const *utetra, eik3_s const *eik);
int utetra_get_num_interior_coefs(utetra_s const *utetra);
void utetra_get_update_inds(utetra_s const *utetra, size_t l[3]);
bool utetra_has_shadow_solution(utetra_s const *utetra, eik3_s const *eik);
size_t utetra_get_num_shared_inds(utetra_s const *u1, utetra_s const *u2);
size_t utetra_get_shared_inds(utetra_s const *u1, utetra_s const *u2, size_t *l);
bool utetra_contains_inds(utetra_s const *u, size_t const *l, size_t n);
bool utetras_yield_same_update(utetra_s const **u, size_t n);
size_t utetra_get_l(utetra_s const *utetra);
void utetra_set_update_inds(utetra_s *utetra, size_t l[3]);
void utetra_set_l0(utetra_s *utetra, size_t l0);
void utetra_set_l1(utetra_s *utetra, size_t l1);
void utetra_set_l2(utetra_s *utetra, size_t l2);
bool utetra_opt_inc_on_other_utetra(utetra_s const *u1, utetra_s const *u2);
void utetra_get_x(utetra_s const *u, dbl x[3]);
void utetra_get_interior_coefs(utetra_s const *utetra, size_t *l);
size_t utetra_get_active_inds(utetra_s const *utetra, size_t l[2]);
par3_s utetra_get_parent(utetra_s const *utetra);
bool utetra_get_point_for_index(utetra_s const *utetra, size_t l, dbl x[3]);
bool utetra_get_op_ind(utetra_s const *utetra, size_t const le[2], size_t *l);

#ifdef __cplusplus
}
#endif
