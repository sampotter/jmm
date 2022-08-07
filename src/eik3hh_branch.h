#pragma once

#include "grid2.h"

typedef enum eik3hh_branch_type {
  EIK3HH_BRANCH_TYPE_UNINITIALIZED,
  EIK3HH_BRANCH_TYPE_PT_SRC,
  EIK3HH_BRANCH_TYPE_REFL,
} eik3hh_branch_type_e;

void eik3hh_branch_alloc(eik3hh_branch_s **branch);
void eik3hh_branch_init_pt_src(eik3hh_branch_s *branch, eik3hh_s const *hh,
                               dbl3 const xsrc);
void eik3hh_branch_deinit(eik3hh_branch_s *branch, bool free_children);
void eik3hh_branch_dealloc(eik3hh_branch_s **hh);
void eik3hh_branch_solve(eik3hh_branch_s *branch);
eik3_s *eik3hh_branch_get_eik(eik3hh_branch_s *branch);
void eik3hh_branch_dump_xy_slice(eik3hh_branch_s const *branch,
                                 grid2_to_mesh3_mapping_s const *mapping,
                                 field_e field, char const *path);
