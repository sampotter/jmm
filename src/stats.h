#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct runstd {
  uint64_t k;
  dbl M, Mprev;
  dbl S, Sprev;
} runstd_s;

void runstd_init(runstd_s *runstd);
void runstd_update(runstd_s *runstd, dbl x);
dbl runstd_get_mean(runstd_s const *runstd);
dbl runstd_get_std(runstd_s const *runstd);

#ifdef __cplusplus
}
#endif
