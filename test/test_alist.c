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
  assert_that(alist_get_by_key(lst, &k, &v), is_false);

  assert_that(alist_get_by_index(lst, 0, &v));
  assert_that(v, is_equal_to(1.0));

  assert_that(alist_get_by_index(lst, 1, &v));
  assert_that(v, is_equal_to(2.0));

  assert_that(alist_get_by_index(lst, 2, &v));
  assert_that(v, is_equal_to(3.0));

  assert_that(alist_get_by_index(lst, 3, &v));
  assert_that(v, is_equal_to(4.0));

  k = 4;
  assert_that(alist_get_by_key(lst, &k, &v));
  assert_that(v, is_equal_to(1.0));

  k = 3;
  assert_that(alist_get_by_key(lst, &k, &v));
  assert_that(v, is_equal_to(2.0));

  k = 2;
  assert_that(alist_get_by_key(lst, &k, &v));
  assert_that(v, is_equal_to(3.0));

  k = 1;
  assert_that(alist_get_by_key(lst, &k, &v));
  assert_that(v, is_equal_to(4.0));

  v = -10.0;
  alist_set_by_index(lst, 0, &v);
  assert_that(alist_get_by_index(lst, 0, &v));
  assert_that(v, is_equal_to(-10.0));

  assert_that(alist_get_key(lst, 0, &k));
  assert_that(k, is_equal_to(4));

  assert_that(alist_get_key(lst, 1, &k));
  assert_that(k, is_equal_to(3));

  assert_that(alist_get_key(lst, 2, &k));
  assert_that(k, is_equal_to(2));

  assert_that(alist_get_key(lst, 3, &k));
  assert_that(k, is_equal_to(1));

  assert_that(alist_get_key(lst, 4, &k), is_false);
  assert_that(alist_get_key(lst, 5, &k), is_false);

  assert_that(alist_get_pair(lst, 0, &k, &v));
  assert_that(k, is_equal_to(4));
  assert_that(v, is_equal_to(-10.0));

  assert_that(alist_get_pair(lst, 1, &k, &v));
  assert_that(k, is_equal_to(3));
  assert_that(v, is_equal_to(2.0));

  assert_that(alist_get_pair(lst, 2, &k, &v));
  assert_that(k, is_equal_to(2));
  assert_that(v, is_equal_to(3.0));

  assert_that(alist_get_pair(lst, 3, &k, &v));
  assert_that(k, is_equal_to(1));
  assert_that(v, is_equal_to(4.0));

  assert_that(alist_get_pair(lst, 4, &k, &v), is_false);
  assert_that(alist_get_pair(lst, 5, &k, &v), is_false);

  alist_deinit(lst);
  alist_dealloc(&lst);
}
