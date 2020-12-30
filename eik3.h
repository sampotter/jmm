#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "jet.h"

typedef struct mesh3 mesh3_s;

/**
 * A structure managing a jet marching method solving the eikonal
 * equation in 3D on an unstructured tetrahedron mesh.
 *
 * NOTE: this is just for s = 1 at the moment. Will extend this to
 * handle s != later.
 */
typedef struct eik3 eik3_s;

void eik3_alloc(eik3_s **eik);
void eik3_dealloc(eik3_s **eik);
void eik3_init(eik3_s *eik, mesh3_s const *mesh);
void eik3_deinit(eik3_s *eik);
size_t eik3_peek(eik3_s const *eik);
size_t eik3_step(eik3_s *eik);
void eik3_solve(eik3_s *eik);
void eik3_add_trial(eik3_s *eik, size_t ind, jet3 jet);
void eik3_add_valid(eik3_s *eik, size_t ind, jet3 jet);
bool eik3_is_point_source(eik3_s const *eik, size_t ind);
bool eik3_is_far(eik3_s const *eik, size_t ind);
bool eik3_is_trial(eik3_s const *eik, size_t ind);
bool eik3_is_valid(eik3_s const *eik, size_t ind);
mesh3_s const *eik3_get_mesh(eik3_s const *eik);
jet3 eik3_get_jet(eik3_s const *eik, size_t l);
jet3 *eik3_get_jet_ptr(eik3_s const *eik);
state_e *eik3_get_state_ptr(eik3_s const *eik);
int eik3_get_num_full_updates(eik3_s const *eik);
int *eik3_get_full_update_ptr(eik3_s const *eik);

#ifdef __cplusplus
}
#endif
