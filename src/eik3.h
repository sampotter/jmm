#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "array.h"
#include "bb.h"
#include "common.h"
#include "error.h"
#include "jet.h"
#include "par.h"

// TODO: put functions we want to go in the "public API" here

void eik3_alloc(eik3_s **eik);
void eik3_dealloc(eik3_s **eik);
void eik3_init(eik3_s *eik, mesh3_s const *mesh);
void eik3_deinit(eik3_s *eik);
bool eik3_is_initialized(eik3_s const *eik);

void eik3_dump_jet(eik3_s const *eik, char const *path);
void eik3_dump_state(eik3_s const *eik, char const *path);
void eik3_dump_par_l(eik3_s const *eik, char const *path);
void eik3_dump_par_b(eik3_s const *eik, char const *path);
void eik3_dump_accepted(eik3_s const *eik, char const *path);
void eik3_dump_has_bc(eik3_s const *eik, char const *path);

void eik3_do_diff_utri(eik3_s *eik, size_t l, size_t l0, size_t l1);
void eik3_do_utetra(eik3_s *eik, size_t l, size_t l0, size_t l1, size_t l2);

size_t eik3_peek(eik3_s const *eik);
jmm_error_e eik3_step(eik3_s *eik, size_t *l0);
jmm_error_e eik3_solve(eik3_s *eik);
bool eik3_is_solved(eik3_s const *eik);

bool eik3_is_far(eik3_s const *eik, size_t ind);
bool eik3_is_trial(eik3_s const *eik, size_t ind);
bool eik3_is_valid(eik3_s const *eik, size_t ind);

mesh3_s const *eik3_get_mesh(eik3_s const *eik);
array_s const *eik3_get_trial_inds(eik3_s const *eik);
array_s const *eik3_get_bc_inds(eik3_s const *eik);

dbl eik3_get_T(eik3_s const *eik, size_t l);
jet31t eik3_get_jet(eik3_s const *eik, size_t l);
jet31t *eik3_get_jet_ptr(eik3_s const *eik);
state_e *eik3_get_state_ptr(eik3_s const *eik);
par3_s eik3_get_par(eik3_s const *eik, size_t l);
bool eik3_has_par(eik3_s const *eik, size_t l);
bool eik3_has_BCs(eik3_s const *eik, size_t l);
size_t const *eik3_get_accepted_ptr(eik3_s const *eik);
size_t eik3_num_bc(eik3_s const *eik);

void eik3_add_trial(eik3_s *eik, size_t l, jet31t jet);
void eik3_add_bc(eik3_s *eik, size_t l, jet31t jet);
void eik3_add_pt_src_bcs(eik3_s *eik, dbl3 const xsrc, dbl tau0);
void eik3_add_diff_bcs(eik3_s *eik, eik3_s const *eik_in, size_t diff_index);
void eik3_add_refl_bcs(eik3_s *eik, eik3_s const *eik_in, size_t refl_index);
bool eik3_has_diff_bc(eik3_s const *eik, size_t const le[2]);
void eik3_get_diff_bc(eik3_s const *eik, size_t const le[2], bb31 *T);
void eik3_set_jet(eik3_s *eik, size_t l, jet31t jet);
void eik3_set_par(eik3_s *eik, size_t l, par3_s par);
void eik3_get_edge_T(eik3_s const *eik, size_t const le[2], bb31 *T);
bool eik3_updated_from_diff_edge(eik3_s const *eik, size_t l);

dbl eik3_get_max_T(eik3_s const *eik);

void eik3_get_origins(eik3_s const *eik, dbl *origin);
void eik3_get_D2T(eik3_s const *eik, dbl33 *D2T);
void eik3_get_spreading_factor(eik3_s const *eik, dbl33 const *D2T, dbl *spread);
void eik3_get_t_in(eik3_s const *eik, dbl3 *t_in);
void eik3_get_t_out(eik3_s const *eik, dbl3 *t_out);

#ifdef __cplusplus
}
#endif
