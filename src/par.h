#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "def.h"

typedef struct par2 {
  size_t l[2];
  dbl b[2];
} par2_s;

void par2_init_empty(par2_s *par);
bool par2_is_empty(par2_s const *par);

typedef struct par3 {
  size_t l[3];
  dbl b[3];
} par3_s;

par3_s make_par3(size_t l[3], dbl b[3]);
void par3_init_empty(par3_s *par);
void par3_set(par3_s *par, size_t const *l, dbl const *b, int n);
size_t par3_size(par3_s const *par);
void par3_get_xb(par3_s const *par, mesh3_s const *mesh, dbl xb[3]);
bool par3_is_empty(par3_s const *par);
size_t par3_num_active(par3_s const *par);
size_t par3_get_active_inds(par3_s const *par, size_t l[3]);
size_t par3_get_active_and_inactive_inds(par3_s const *par, uint3 la, uint3 li);
size_t par3_get_active(par3_s const *par, size_t *l, dbl *b);

#ifdef __cplusplus
}
#endif
