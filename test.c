#include "heap.h"

#include <assert.h>
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
  return hypot(xy.x, xy.y);
}

dvec2 grad_one(dvec2 p) {
  dvec2 zero = {0, 0};
  return zero;
}

void sjs_basic() {
  ivec2 shape = {2, 2};
  ivec2 ind = {1, 1};
  dbl h = 1, r = 1;

  sjs_s *sjs;
  sjs_alloc(&sjs);
  sjs_init(sjs, shape, h, one, grad_one);
  int nf = sjs_add_fac_pt_src(sjs, ind, r);
  assert(nf == 4);
  sjs_solve(sjs);
  sjs_deinit(sjs);
  sjs_dealloc(&sjs);
}

int main() {
  heap_basic();
}
