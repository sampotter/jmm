#include <cgreen/cgreen.h>

#include "edgemap.h"

Describe(edgemap);
BeforeEach(edgemap) {}
AfterEach(edgemap) {}

Ensure (edgemap, basic_test) {
  edgemap_s *edgemap;
  edgemap_alloc(&edgemap);
  edgemap_init(edgemap, sizeof(int));

  int value = 1;
  edgemap_set(edgemap, make_edge(1, 0), &value);

  value = 2;
  edgemap_set(edgemap, make_edge(2, 0), &value);

  value = 3;
  edgemap_set(edgemap, make_edge(1, 2), &value);

  edgemap_get(edgemap, make_edge(0, 1), &value);
  assert_that(value, is_equal_to(1));

  edgemap_deinit(edgemap);
  edgemap_dealloc(&edgemap);
}
