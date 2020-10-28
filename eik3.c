#include "eik3.h"

#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

struct eik3 {
  field3_s const *s;
  mesh3_s const *mesh;
  state_e *state;
};

void eik3_alloc(eik3_s **eik) {
  *eik = malloc(sizeof(eik3_s));
}

void eik3_dealloc(eik3_s **eik) {
  free(*eik);
  *eik = NULL;
}

// void eik3_init(eik3_s *eik, field3_s const *s, ivec3 shape, dvec3 pmin, dbl h) {

// }
