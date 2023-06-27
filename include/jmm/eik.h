#pragma once

#include <stdbool.h>

#include "bicubic.h"
#include "common.h"
#include "grid2.h"
#include "heap.h"
#include "jet.h"
#include "par.h"
#include "vec.h"

typedef struct eik eik_s;

void eik_alloc(eik_s **eik);
void eik_dealloc(eik_s **eik);
void eik_init(eik_s *eik, field2_s const *slow, grid2_s const *grid);
void eik_deinit(eik_s *eik);
size_t eik_peek(eik_s const *eik);
void eik_step(eik_s *eik);
void eik_solve(eik_s *eik);
void eik_add_trial(eik_s *eik, int2 ind, jet21p jet);
void eik_add_valid(eik_s *eik, int2 ind, jet21p jet);
void eik_make_bd(eik_s *eik, int2 ind);
void eik_get_shape(eik_s const *eik, int2 shape);
jet21p eik_get_jet(eik_s *eik, int2 ind);
jet21p *eik_get_jets_ptr(eik_s const *eik);
state_e eik_get_state(eik_s const *eik, int2 ind);
state_e *eik_get_states_ptr(eik_s const *eik);
dbl eik_T(eik_s const *eik, dbl2 xy);
dbl eik_Tx(eik_s const *eik, dbl2 xy);
dbl eik_Ty(eik_s const *eik, dbl2 xy);
dbl eik_Txx(eik_s const *eik, dbl2 xy);
dbl eik_Txy(eik_s const *eik, dbl2 xy);
dbl eik_Tyy(eik_s const *eik, dbl2 xy);
bool eik_can_build_cell(eik_s const *eik, int2 indc);
void eik_build_cells(eik_s *eik);
bicubic_s eik_get_bicubic(eik_s const *eik, int2 indc);
bicubic_s *eik_get_bicubics_ptr(eik_s const *eik);
heap_s *eik_get_heap(eik_s const *eik);
par2_s eik_get_par(eik_s const *eik, int2 ind);
bool eik_has_par(eik_s const *eik, int2 ind);
size_t const *eik_get_accepted_ptr(eik_s const *eik);
