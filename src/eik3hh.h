#pragma once

#include "jet.h"
#include "mesh3.h"

void eik3hh_alloc(eik3hh_s **hh);
void eik3hh_init_pt_src(eik3hh_s *hh, mesh3_s const *mesh,
                                dbl3 const xsrc, dbl rfac);
void eik3hh_deinit(eik3hh_s *hh);
void eik3hh_dealloc(eik3hh_s **hh);
void eik3hh_solve(eik3hh_s *hh);
size_t eik3hh_num_bc(eik3hh_s const *hh);
void eik3hh_dump_xy_slice(eik3hh_s const *hh,
                          grid2_to_mesh3_mapping_s const *mapping,
                          field_e field, char const *path);
eik3_s const *eik3hh_get_eik_ptr(eik3hh_s const *hh);
jet31t const *eik3hh_get_jet_ptr(eik3hh_s const *hh);
