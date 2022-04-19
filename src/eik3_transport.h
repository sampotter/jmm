#pragma once

#include "eik3.h"

void eik3_transport_dbl(eik3_s const *eik, dbl *values, bool skip_filled);
void eik3_transport_dblz(eik3_s const *eik, dblz *values, bool skip_filled);
void eik3_transport_curvature(eik3_s const *eik, dbl *kappa, bool skip_filled);
