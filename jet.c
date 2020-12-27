#include "jet.h"

#include <math.h>
#include <string.h>

bool jet3_eq(jet3 const *jet1, jet3 const *jet2) {
  return !memcmp((void *)jet1, (void *)jet2, sizeof(jet3));
}

bool jet3_is_finite(jet3 const *jet) {
  return isfinite(jet->f)
    && isfinite(jet->fx) && isfinite(jet->fy) && isfinite(jet->fz);
}
