#include <cgreen/cgreen.h>

#include "alist.h"

Describe(alist);
BeforeEach(alist) {}
AfterEach(alist) {}

Ensure (alist, basic_test) {
  alist_s *lst;
  alist_alloc(&lst);
  alist_init(lst, sizeof(int), sizeof(float), 8);

  assert_that(alist_is_empty(lst));

  int k;
  float v;

  k = 4;
  v = 1.0;
  alist_append(lst, &k, &v);
  assert_that(alist_is_empty(lst), is_false);
  assert_that(alist_size(lst), is_equal_to(1));

  k = 3;
  v = 2.0;
  alist_append(lst, &k, &v);
  assert_that(alist_size(lst), is_equal_to(2));

  k = 2;
  v = 3.0;
  alist_append(lst, &k, &v);
  assert_that(alist_size(lst), is_equal_to(3));

  k = 1;
  v = 4.0;
  alist_append(lst, &k, &v);
  assert_that(alist_size(lst), is_equal_to(4));

  k = 4;
  assert_that(alist_contains(lst, &k));
  assert_that(alist_find(lst, &k), is_equal_to(0));

  k = 3;
  assert_that(alist_contains(lst, &k));
  assert_that(alist_find(lst, &k), is_equal_to(1));

  k = 2;
  assert_that(alist_contains(lst, &k));
  assert_that(alist_find(lst, &k), is_equal_to(2));

  k = 1;
  assert_that(alist_contains(lst, &k));
  assert_that(alist_find(lst, &k), is_equal_to(3));

  k = -10;
  assert_that(alist_contains(lst, &k), is_false);

  alist_get_by_index(lst, 0, &v);
  assert_that(v, is_equal_to(1.0));

  alist_get_by_index(lst, 1, &v);
  assert_that(v, is_equal_to(2.0));

  alist_get_by_index(lst, 2, &v);
  assert_that(v, is_equal_to(3.0));

  alist_get_by_index(lst, 3, &v);
  assert_that(v, is_equal_to(4.0));

  k = 4;
  alist_get_by_key(lst, &k, &v);
  assert_that(v, is_equal_to(1.0));

  k = 3;
  alist_get_by_key(lst, &k, &v);
  assert_that(v, is_equal_to(2.0));

  k = 2;
  alist_get_by_key(lst, &k, &v);
  assert_that(v, is_equal_to(3.0));

  k = 1;
  alist_get_by_key(lst, &k, &v);
  assert_that(v, is_equal_to(4.0));

  v = -10.0;
  alist_set_by_index(lst, 0, &v);
  alist_get_by_index(lst, 0, &v);
  assert_that(v, is_equal_to(-10.0));

  alist_deinit(lst);
  alist_dealloc(&lst);
}
