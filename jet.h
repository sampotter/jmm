#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct jet {
  dbl f, fx, fy, fxy;
} jet_s;

typedef struct {
  dbl f;
  dbl fx, fy, fz;
} jet3;

#ifdef __cplusplus
}
#endif
