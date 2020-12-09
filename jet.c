#include "jet.h"

#include <math.h>

bool jet3_is_finite(jet3 const *jet) {
  return isfinite(jet->f)
    && isfinite(jet->fx) && isfinite(jet->fy) && isfinite(jet->fz);
}
