#pragma once

#include "3d_wedge.h"

#ifdef __cplusplus
extern "C"
#endif
jmm_error_e mesh3_init_from_3d_wedge_spec(mesh3_s *mesh,
                                          jmm_3d_wedge_spec_s const *spec);
