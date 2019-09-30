#include "heap.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

void heap_basic() {
  heap_s *heap = NULL;
  heap_alloc(&heap);

  dbl values[4] = {4, 3, 2, 1}, *values_ptr;
  values_ptr = values;

  int positions[4] = {-1, -1, -1, -1}, *positions_ptr;
  positions_ptr = positions;

  value_b value = ^(int l) {
    return values_ptr[l];
  };

  setpos_b setpos = ^(int l, int pos) {
    positions_ptr[l] = pos;
  };

  heap_init(heap, 4, value, setpos);
  assert(heap_size(heap) == 0);

  heap_insert(heap, 0);
  assert(heap_size(heap) == 1);
  assert(heap_front(heap) == 0);

  heap_insert(heap, 1);
  assert(heap_size(heap) == 2);
  assert(heap_front(heap) == 1);

  heap_insert(heap, 2);
  assert(heap_size(heap) == 3);
  assert(heap_front(heap) == 2);

  heap_insert(heap, 3);
  assert(heap_size(heap) == 4);
  assert(heap_front(heap) == 3);

  heap_pop(heap);
  assert(heap_front(heap) == 2);
  assert(heap_size(heap) == 3);

  heap_pop(heap);
  assert(heap_front(heap) == 1);
  assert(heap_size(heap) == 2);

  heap_pop(heap);
  assert(heap_front(heap) == 0);
  assert(heap_size(heap) == 1);

  heap_pop(heap);
  assert(heap_size(heap) == 0);

  heap_deinit(heap);
  heap_dealloc(&heap);

  assert(heap == NULL);
}

dbl one(dvec2 xy) {
  return 1.0;
}

dvec2 grad_one(dvec2 p) {
  dvec2 zero = {0, 0};
  return zero;
}

bool approx_eq(dbl x, dbl y, dbl tol) {
  return fabs(x - y) < tol*fmax(fabs(x), fabs(y)) + tol;
}

void sjs_basic() {
  ivec2 shape = {3, 3}, ind = {1, 1};
  dvec2 xymin = {-1, -1};
  dbl h = 1, r = SQRT2;

  sjs_s *sjs;
  sjs_alloc(&sjs);
  sjs_init(sjs, shape, xymin, h, one, grad_one);

  int nf, nfc;
  sjs_add_fac_pt_src(sjs, ind, r, &nf, &nfc);
  assert(nf == 9);
  assert(nfc == 4);

  sjs_solve(sjs);

  void (^check)(dbl, dbl) = ^(dbl x, dbl y) {
    dvec2 xy = {x, y};
    dbl T = hypot(x, y);
    dbl Tx = x/T;
    dbl Ty = y/T;
    dbl Txy = -x*y/pow(T, 3);
    assert(approx_eq(sjs_T(sjs, xy), T, EPS));
    assert(approx_eq(sjs_Tx(sjs, xy), Tx, EPS));
    assert(approx_eq(sjs_Ty(sjs, xy), Ty, EPS));
    assert(approx_eq(sjs_Txy(sjs, xy), Txy, EPS));
  };

  check(-1, -1);
  check(-1,  0);
  check(-1,  1);
  check( 0, -1);

  // Handle (0, 0) as special case
  {
    dvec2 xy = {0, 0};
    assert(approx_eq(sjs_T(sjs, xy), 0, EPS));
    assert(isnan(sjs_Tx(sjs, xy)));
    assert(isnan(sjs_Ty(sjs, xy)));
    assert(isnan(sjs_Txy(sjs, xy)));
  }

  check( 0,  1);
  check( 1, -1);
  check( 1,  0);
  check( 1,  1);

  sjs_deinit(sjs);
  sjs_dealloc(&sjs);
}

int main() {
  heap_basic();
  sjs_basic();
}
