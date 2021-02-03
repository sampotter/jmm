#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "jet.h"

typedef struct mesh3 mesh3_s;
typedef struct grid3 grid3_s;

void xfer(mesh3_s const *mesh, jet3 const *jet, grid3_s const *grid, dbl *y);

#ifdef __cplusplus
}
#endif
