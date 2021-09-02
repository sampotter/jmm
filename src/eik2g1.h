#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "grid2.h"
#include "jet.h"

void eik2g1_alloc(eik2g1_s **eik);
void eik2g1_dealloc(eik2g1_s **eik);
void eik2g1_init(eik2g1_s *eik, grid2_s const *grid);
void eik2g1_deinit(eik2g1_s *eik);
size_t eik2g1_peek(eik2g1_s const *eik);
size_t eik2g1_step(eik2g1_s *eik);
void eik2g1_solve(eik2g1_s *eik);
bool eik2g1_is_solved(eik2g1_s const *eik);
void eik2g1_add_trial(eik2g1_s *eik, int2 const ind, jet22t jet);
void eik2g1_add_valid(eik2g1_s *eik, int2 const ind, jet22t jet);
state_e const *eik2g1_get_state_ptr(eik2g1_s const *eik);
jet22t const *eik2g1_get_jet_ptr(eik2g1_s const *eik);
