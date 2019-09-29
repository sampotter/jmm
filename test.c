#include <cgreen/cgreen.h>

#include "heap.h"

Describe(heap);
BeforeEach(heap) {}
AfterEach(heap) {}

Ensure (heap, basic) {
  heap_s *heap = NULL;
  heap_alloc(&heap);
  assert_that(heap, is_not_null);
  assert_that(heap, is_null);
}
