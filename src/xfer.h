#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "def.h"
#include "jet.h"

void xfer(mesh3_s const *mesh, jet3 const *jet, grid3_s const *grid, dbl *y);

#ifdef __cplusplus
}
#endif
