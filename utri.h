#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "jet.h"

typedef struct eik3 eik3_s;

typedef struct utri utri_s;

void utri_alloc(utri_s **utri);
void utri_dealloc(utri_s **utri);
void utri_init_from_eik3(utri_s *utri, eik3_s const *eik, size_t l,
                         size_t l0, size_t l1);
void utri_init(utri_s *utri, dbl const x[3], dbl const Xt[2][3],
               jet3 const jet[2]);
bool utri_is_causal(utri_s const *utri);
void utri_solve(utri_s *utri);
dbl utri_get_lambda(utri_s const *utri);
void utri_get_bary_coords(utri_s const *utri, dbl b[2]);
dbl utri_get_value(utri_s const *utri);
void utri_get_jet(utri_s const *utri, jet3 *jet);

#ifdef __cplusplus
}
#endif
