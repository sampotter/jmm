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

bool jet3_eq(jet3 const *jet1, jet3 const *jet2);
bool jet3_is_finite(jet3 const *jet);

#ifdef __cplusplus
}
#endif
