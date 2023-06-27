#pragma once

#include "common.h"
#include "mesh22.h"
#include "jet.h"

void eik2mp_alloc(eik2mp_s **eik);
void eik2mp_dealloc(eik2mp_s **eik);
void eik2mp_init(eik2mp_s *eik, mesh22_s const *mesh);
void eik2mp_deinit(eik2mp_s *eik);
size_t eik2mp_peek(eik2mp_s const *eik);
size_t eik2mp_step(eik2mp_s *eik);
void eik2mp_solve(eik2mp_s *eik);
bool eik2mp_is_solved(eik2mp_s const *eik);
void eik2mp_add_trial(eik2mp_s *eik, size_t l, jet21t jet);
void eik2mp_add_valid(eik2mp_s *eik, size_t l, jet21t jet);
bool eik2mp_is_valid(eik2mp_s const *eik, size_t l);
bool eik2mp_is_trial(eik2mp_s const *eik, size_t l);
bool eik2mp_is_far(eik2mp_s const *eik, size_t l);
state_e const *eik2mp_get_state_ptr(eik2mp_s const *eik);
jet21t const *eik2mp_get_jet_ptr(eik2mp_s const *eik);
