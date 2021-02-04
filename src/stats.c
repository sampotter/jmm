#include "stats.h"

#include <math.h>

void runstd_init(runstd_s *runstd) {
  runstd->k = 0;
}

void runstd_update(runstd_s *runstd, dbl x) {
  if (++runstd->k == 1) {
    runstd->M = runstd->Mprev = x;
    runstd->Sprev = 0;
  } else {
    runstd->M = runstd->Mprev + (x - runstd->Mprev)/runstd->k;
    runstd->S = runstd->Sprev + (x - runstd->Mprev)*(x - runstd->M);
    runstd->Mprev = runstd->M;
    runstd->Sprev = runstd->S;
  }
}

dbl runstd_get_mean(runstd_s const *runstd) {
  return runstd->k == 0 ? 0 : runstd->M;
}

dbl runstd_get_var(runstd_s const *runstd) {
  return runstd->k <= 1 ? 0 : runstd->S/(runstd->k - 1);
}

dbl runstd_get_std(runstd_s const *runstd) {
  return sqrt(runstd_get_var(runstd));
}
