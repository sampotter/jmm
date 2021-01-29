#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct par3 {
  size_t l[3];
  dbl b[3];
} par3_s;

par3_s make_par3(size_t l[3], dbl b[3]);
void par3_init_empty(par3_s *par);
void par3_set(par3_s *par, size_t const *l, dbl const *b, int n);
int par3_size(par3_s const *par);

#ifdef __cplusplus
}
#endif
