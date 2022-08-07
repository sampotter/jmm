#pragma once

#include "camera.h"
#include "jet.h"
#include "mesh3.h"

void eik3hh_alloc(eik3hh_s **hh);
void eik3hh_init_with_pt_src(eik3hh_s *hh, mesh3_s const *mesh, dbl c,
                             dbl rfac, dbl3 const xsrc);
void eik3hh_deinit(eik3hh_s *hh);
void eik3hh_dealloc(eik3hh_s **hh);
void eik3hh_add_pt_src(eik3hh_s *hh, dbl3 const xsrc);
mesh3_s const *eik3hh_get_mesh(eik3hh_s const *hh);
dbl eik3hh_get_rfac(eik3hh_s const *hh);
eik3hh_branch_s *eik3hh_get_root_branch(eik3hh_s *hh);
void eik3hh_render_frames(eik3hh_s const *hh, camera_s const *camera,
                          dbl t0, dbl t1, dbl frame_rate, bool verbose);
