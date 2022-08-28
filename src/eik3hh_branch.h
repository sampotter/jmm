#pragma once

#include "array.h"
#include "camera.h"
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
void eik3hh_branch_solve(eik3hh_branch_s *branch, bool verbose);
bool eik3hh_branch_is_solved(eik3hh_branch_s const *branch);
eik3_s *eik3hh_branch_get_eik(eik3hh_branch_s *branch);
array_s *eik3hh_branch_get_children(eik3hh_branch_s *branch);
dbl const *eik3hh_branch_get_spread(eik3hh_branch_s const *branch);
dbl const *eik3hh_branch_get_org(eik3hh_branch_s const *branch);
size_t eik3hh_branch_get_earliest_refl(eik3hh_branch_s const *branch);
array_s *eik3hh_branch_get_visible_refls(eik3hh_branch_s const *branch);
eik3hh_branch_s *eik3hh_branch_add_refl(eik3hh_branch_s const *branch,
                                        size_t refl_index);
void eik3hh_branch_dump_jet(eik3hh_branch_s const *branch, char const *path);
void eik3hh_branch_dump_org(eik3hh_branch_s const *branch, char const *path);
void eik3hh_branch_dump_spread(eik3hh_branch_s const *branch, char const *path);
void eik3hh_branch_dump_xy_slice(eik3hh_branch_s const *branch,
                                 grid2_to_mesh3_mapping_s const *mapping,
                                 field_e field, char const *path);
void eik3hh_branch_render_frames(eik3hh_branch_s const *branch,
                                 camera_s const *camera,
                                 dbl t0, dbl t1, dbl frame_rate,
                                 bool verbose);
