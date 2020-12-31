#include <cgreen/cgreen.h>

#include "mesh3.h"

Describe(mesh3);
BeforeEach(mesh3) {}
AfterEach(mesh3) {}

#define SET_UP_CUBE_MESH()                      \
  dbl verts[24] = {                             \
    0, 0, 0,                                    \
    0, 0, 1,                                    \
    0, 1, 0,                                    \
    0, 1, 1,                                    \
    1, 0, 0,                                    \
    1, 0, 1,                                    \
    1, 1, 0,                                    \
    1, 1, 1                                     \
  };                                            \
  size_t cells[20] = {                          \
    0, 2, 3, 6,                                 \
    0, 4, 5, 6,                                 \
    0, 3, 5, 6,                                 \
    0, 1, 3, 5,                                 \
    3, 5, 6, 7                                  \
  };                                            \
  mesh3_s *mesh;                                \
  mesh3_alloc(&mesh);                           \
  mesh3_init(mesh, verts, 8, cells, 5);

#define TEAR_DOWN_MESH()                        \
  mesh3_deinit(mesh);                           \
  mesh3_dealloc(&mesh);

Ensure (mesh3, nvc_works) {
  SET_UP_CUBE_MESH();

  int nvc_gt[8] = {4, 1, 1, 4, 1, 4, 4, 1};

  int nvc;
  for (int i = 0; i < 8; ++i) {
    nvc = mesh3_nvc(mesh, i);
    assert_that(nvc, is_equal_to(nvc_gt[i]));
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, nvv_works) {
  SET_UP_CUBE_MESH();

  int nvv_gt[8] = {6, 3, 3, 6, 3, 6, 6, 3};

  int nvv;
  for (int i = 0; i < 8; ++i) {
    nvv = mesh3_nvv(mesh, i);
    assert_that(nvv, is_equal_to(nvv_gt[i]));
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, nve_works) {
  SET_UP_CUBE_MESH();



  TEAR_DOWN_MESH();
}
