#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "jet.h"
#include "par.h"

typedef struct utri utri_s;

void utri_alloc(utri_s **utri);
void utri_dealloc(utri_s **utri);
void utri_init_from_eik3(utri_s *utri, eik3_s const *eik, size_t l,
                         size_t l0, size_t l1);
void utri_init_from_eik3_without_l(utri_s *utri, eik3_s const *eik,
                                   dbl const x[3], size_t l0, size_t l1);
void utri_init(utri_s *utri, dbl const x[3], dbl const Xt[2][3],
               jet3 const jet[2]);
void utri_solve(utri_s *utri);
par3_s utri_get_par(utri_s const *u);
dbl utri_get_value(utri_s const *utri);
void utri_get_jet(utri_s const *utri, jet3 *jet);
bool utri_update_ray_is_physical(utri_s const *utri, eik3_s const *eik);
int utri_cmp(utri_s const **h1, utri_s const **h2);
bool utri_has_interior_point_solution(utri_s const *utri);
void utri_set_orig_index(utri_s *utri, int i);
int utri_get_orig_index(utri_s const *utri);
bool utri_is_finite(utri_s const *utri);

bool utris_yield_same_update(utri_s const *utri1, utri_s const *utri2);

#if JMM_TEST
bool utri_is_causal(utri_s const *utri);
dbl utri_get_lambda(utri_s const *utri);
#endif

#ifdef __cplusplus
}
#endif
