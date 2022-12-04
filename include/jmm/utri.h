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
void utri_init(utri_s *u, eik3_s const *eik, size_t lhat, size_t const l[2]);
void utri_solve(utri_s *utri);
par3_s utri_get_par(utri_s const *u);
dbl utri_get_value(utri_s const *utri);
void utri_get_jet31t(utri_s const *utri, jet31t *jet);
bool utri_ray_is_occluded(utri_s const *utri, eik3_s const *eik);
bool utri_has_interior_point_solution(utri_s const *utri);
bool utri_active_vert_is_terminal_diff_vert(utri_s const *utri,eik3_s const *eik);
bool utri_is_backwards(utri_s const *utri, eik3_s const *eik);
size_t utri_get_active_ind(utri_s const *utri);
size_t utri_get_inactive_ind(utri_s const *utri);
size_t utri_get_l(utri_s const *utri);
bool utri_is_degenerate(utri_s const *u);
bool utri_has_inds(utri_s const *u, size_t lhat, uint2 const l);

bool utris_have_same_inds(utri_s const *u1, utri_s const *u2);

#if JMM_TEST
bool utri_is_causal(utri_s const *utri);
dbl utri_get_lambda(utri_s const *utri);
#endif

#ifdef __cplusplus
}
#endif
