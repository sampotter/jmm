#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

// TODO: make it so that the verts/cells can be copied or not (in
// which case mesh3 is just a "view")

typedef struct mesh3 mesh3_s;

void mesh3_alloc(mesh3_s **mesh);
void mesh3_dealloc(mesh3_s **mesh);
void mesh3_init(mesh3_s *mesh,
                dbl const *verts, size_t nverts,
                int const *cells, size_t ncells);
void mesh3_deinit(mesh3_s *mesh);

#ifdef __cplusplus
}
#endif
