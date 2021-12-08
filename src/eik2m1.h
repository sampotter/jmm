#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "mesh22.h"
#include "jet.h"
#include "par.h"

void eik2m1_alloc(eik2m1_s **eik);
void eik2m1_dealloc(eik2m1_s **eik);
void eik2m1_init(eik2m1_s *eik, mesh22_s const *mesh);
void eik2m1_deinit(eik2m1_s *eik);
size_t eik2m1_peek(eik2m1_s const *eik);
size_t eik2m1_step(eik2m1_s *eik);
void eik2m1_solve(eik2m1_s *eik);
bool eik2m1_is_solved(eik2m1_s const *eik);
void eik2m1_add_trial(eik2m1_s *eik, size_t l, jet22t jet);
void eik2m1_add_valid(eik2m1_s *eik, size_t l, jet22t jet);
bool eik2m1_is_valid(eik2m1_s const *eik, size_t l);
bool eik2m1_is_trial(eik2m1_s const *eik, size_t l);
bool eik2m1_is_far(eik2m1_s const *eik, size_t l);
state_e const *eik2m1_get_state_ptr(eik2m1_s const *eik);
jet22t const *eik2m1_get_jet_ptr(eik2m1_s const *eik);
par2_s eik2m1_get_par(eik2m1_s const *eik, size_t l);
par2_s const *eik2m1_get_par_ptr(eik2m1_s const *eik);
