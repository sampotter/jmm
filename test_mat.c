#include <cgreen/cgreen.h>

Describe(mat);
BeforeEach(mat) {}
AfterEach(mat) {}

Ensure(mat, dummy) {
  assert_that(1 == 1);
}
