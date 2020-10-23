#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "bicubic.h"
#include "heap.h"
#include "jet.h"
#include "vec.h"

typedef struct field2 field2_s;

typedef struct eik eik_s;

void eik_alloc(eik_s **eik);
void eik_dealloc(eik_s **eik);
void eik_init(eik_s *eik, field2_s const *slow, ivec2 shape, dvec2 xymin, dbl h);
void eik_deinit(eik_s *eik);
void eik_step(eik_s *eik);
void eik_solve(eik_s *eik);
void eik_add_trial(eik_s *eik, ivec2 ind, jet_s jet);
void eik_add_valid(eik_s *eik, ivec2 ind, jet_s jet);
void eik_make_bd(eik_s *eik, ivec2 ind);
ivec2 eik_get_shape(eik_s const *eik);
jet_s eik_get_jet(eik_s *eik, ivec2 ind);
jet_s *eik_get_jets_ptr(eik_s const *eik);
state_e eik_get_state(eik_s const *eik, ivec2 ind);
state_e *eik_get_states_ptr(eik_s const *eik);
dbl eik_T(eik_s *eik, dvec2 xy);
dbl eik_Tx(eik_s *eik, dvec2 xy);
dbl eik_Ty(eik_s *eik, dvec2 xy);
dbl eik_Txy(eik_s *eik, dvec2 xy);
bool eik_can_build_cell(eik_s const *eik, ivec2 indc);
void eik_build_cells(eik_s *eik);
bicubic_s eik_get_bicubic(eik_s const *eik, ivec2 indc);
bicubic_s *eik_get_bicubics_ptr(eik_s const *eik);
heap_s *eik_get_heap(eik_s const *eik);

#ifdef __cplusplus
}
#endif
