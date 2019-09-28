#pragma once

#include "def.h"

#include <math.h>

typedef struct {
  dbl x;
  dbl y;
} dvec2;

dvec2 dvec2_ccomb(dvec2 v0, dvec2 v1, dbl t);
dbl dvec2_dist(dvec2 v0, dvec2 v1);

typedef struct {
  int i;
  int j;
} ivec2;
