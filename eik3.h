#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "jet.h"

typedef struct field3 field3_s;
typedef struct mesh3 mesh3_s;

/**
 * A structure managing a jet marching method solving the eikonal
 * equation in 3D on an unstructured tetrahedron mesh.
 */
typedef struct eik3 eik3_s;

void eik3_alloc(eik3_s **eik);
void eik3_dealloc(eik3_s **eik);
void eik3_init(eik3_s *eik, field3_s const *s, mesh3_s const *mesh);
void eik3_deinit(eik3_s *eik);
void eik3_step(eik3_s *eik);
void eik3_solve(eik3_s *eik);
void eik3_add_trial(eik3_s *eik, size_t ind, jet3_s jet);
void eik3_add_valid(eik3_s *eik, size_t ind, jet3_s jet);

#ifdef __cplusplus
}
#endif
