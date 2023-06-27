#pragma once

#include "common.h"
#include "jet.h"

typedef struct uline uline_s;

void uline_alloc(uline_s **u);
void uline_dealloc(uline_s **u);
void uline_init(uline_s *u, eik3_s const *eik, size_t lhat, size_t l0);
void uline_init_from_points(uline_s *u, eik3_s const *eik, dbl3 const xhat, dbl3 const x0, dbl tol, dbl T0);
void uline_solve(uline_s *u);
dbl uline_get_value(uline_s const *u);
void uline_get_topt(uline_s const *u, dbl3 topt);
jet31t uline_get_jet(uline_s const *u);
