#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "bicubic.h"
#include "heap.h"
#include "jet.h"
#include "vec.h"

typedef struct sjs sjs_s;

void sjs_alloc(sjs_s **sjs);
void sjs_dealloc(sjs_s **sjs);
void sjs_init(sjs_s *sjs, ivec2 shape, dvec2 xymin, dbl h);
void sjs_deinit(sjs_s *sjs);
void sjs_step(sjs_s *sjs);
void sjs_solve(sjs_s *sjs);
void sjs_add_trial(sjs_s *sjs, ivec2 ind, jet_s jet);
void sjs_add_valid(sjs_s *sjs, ivec2 ind, jet_s jet);
void sjs_make_bd(sjs_s *sjs, ivec2 ind);
ivec2 sjs_get_shape(sjs_s const *sjs);
jet_s sjs_get_jet(sjs_s *sjs, ivec2 ind);
state_e sjs_get_state(sjs_s const *sjs, ivec2 ind);
state_e *sjs_get_states_ptr(sjs_s const *sjs);
dbl sjs_T(sjs_s *sjs, dvec2 xy);
dbl sjs_Tx(sjs_s *sjs, dvec2 xy);
dbl sjs_Ty(sjs_s *sjs, dvec2 xy);
dbl sjs_Txy(sjs_s *sjs, dvec2 xy);
bool sjs_can_build_cell(sjs_s const *sjs, ivec2 indc);
void sjs_build_cells(sjs_s *sjs);
bicubic_s sjs_get_bicubic(sjs_s const *sjs, ivec2 indc);
bicubic_s *sjs_get_bicubics_ptr(sjs_s const *sjs);
heap_s *sjs_get_heap(sjs_s const *sjs);

#ifdef __cplusplus
}
#endif
