#include <cgreen/cgreen.h>

#include "mesh3.h"
#include "util.h"

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
  mesh3_init(mesh, verts, 8, cells, 5, true);

#define TEAR_DOWN_MESH()                        \
  mesh3_deinit(mesh);                           \
  mesh3_dealloc(&mesh);

Ensure (mesh3, get_nverts_works_for_cube) {
  SET_UP_CUBE_MESH();

  assert_that(mesh3_nverts(mesh), is_equal_to(8));

  TEAR_DOWN_MESH();
}

Ensure (mesh3, nvc_works_for_cube) {
  SET_UP_CUBE_MESH();

  int nvc_gt[8] = {4, 1, 1, 4, 1, 4, 4, 1};

  int nvc;
  for (int i = 0; i < 8; ++i) {
    nvc = mesh3_nvc(mesh, i);
    assert_that(nvc, is_equal_to(nvc_gt[i]));
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, nvf_works_for_cube) {
  SET_UP_CUBE_MESH();

  int nvf_gt[8] = {4, 1, 1, 4, 1, 4, 4, 1};

  int nvf;
  for (int i = 0; i < 8; ++i) {
    nvf = mesh3_nvf(mesh, i);
    assert_that(nvf, is_equal_to(nvf_gt[i]));
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, nvv_works_for_cube) {
  SET_UP_CUBE_MESH();

  int nvv_gt[8] = {6, 3, 3, 6, 3, 6, 6, 3};

  int nvv;
  for (int i = 0; i < 8; ++i) {
    nvv = mesh3_nvv(mesh, i);
    assert_that(nvv, is_equal_to(nvv_gt[i]));
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, ncc_works_for_cube) {
  SET_UP_CUBE_MESH();

  int ncc_gt[5] = {1, 1, 4, 1, 1};

  int ncc;
  for (int i = 0; i < 5; ++i) {
    ncc = mesh3_ncc(mesh, i);
    assert_that(ncc, is_equal_to(ncc_gt[i]));
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, cc_works_for_cube) {
  SET_UP_CUBE_MESH();

  int ncc[5] = {1, 1, 4, 1, 1};

  size_t cc_gt[5][4] = {
    {2, NO_INDEX, NO_INDEX, NO_INDEX},
    {2, NO_INDEX, NO_INDEX, NO_INDEX},
    {0, 1, 3, 4},
    {2, NO_INDEX, NO_INDEX, NO_INDEX},
    {2, NO_INDEX, NO_INDEX, NO_INDEX}
  };

  size_t cc[4];
  for (int i = 0; i < 5; ++i) {
    mesh3_cc(mesh, i, cc);
    qsort(cc, ncc[i], sizeof(size_t), (compar_t)compar_size_t);
    assert_that(cc, is_equal_to_contents_of(cc_gt[i], ncc[i]*sizeof(size_t)));
  }

  TEAR_DOWN_MESH();
}

static int find_edge(size_t (*edges)[2], int num_edges, size_t i, size_t j) {
  int k;
  for (k = 0; k < num_edges; ++k) {
    if (edges[k][0] == i && edges[k][1] == j) {
      break;
    }
  }
  return k;
}

Ensure (mesh3, nec_works_for_cube) {
  SET_UP_CUBE_MESH();

  size_t edges[][2] = {
    {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, // (0, *)
    {1, 0}, {1, 3}, {1, 5}, // (1, *)
    {2, 0}, {2, 3}, {2, 6}, // (2, *)
    {3, 0}, {3, 1}, {3, 2}, {3, 5}, {3, 6}, {3, 7}, // (3, *)
    {4, 0}, {4, 5}, {4, 6}, // (4, *)
    {5, 0}, {5, 1}, {5, 3}, {5, 4}, {5, 6}, {5, 7}, // (5, *)
    {6, 0}, {6, 2}, {6, 3}, {6, 4}, {6, 5}, {6, 7}, // (6, *)
    {7, 3}, {7, 5}, {7, 6} // (7, *)
  };

  int nec_gt[] = {
    1, 1, 3, 1, 3, 3, // (0, *)
    1, 1, 1, // (1, *)
    1, 1, 1, // (2, *)
    3, 1, 1, 3, 3, 1, // (3, *)
    1, 1, 1, // (4, *)
    3, 1, 3, 1, 3, 1, // (5, *)
    3, 1, 3, 1, 3, 1, // (6, *)
    1, 1, 1 // (7, *)
  };

  int num_edges = sizeof(edges)/sizeof(edges[0]);

  int nec, k;
  for (size_t i = 0; i < 8; ++i) {
    for (size_t j = 0; j < 8; ++j) {
      k = find_edge(edges, num_edges, i, j);
      nec = mesh3_nec(mesh, i, j);
      if (k < num_edges) {
        assert_that(nec, is_equal_to(nec_gt[k]));
      } else {
        assert_that(nec, is_equal_to(0));
      }
    }
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, ec_works_for_cube) {
  SET_UP_CUBE_MESH();

  size_t edges[][2] = {
    {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, // (0, *)
    {1, 0}, {1, 3}, {1, 5}, // (1, *)
    {2, 0}, {2, 3}, {2, 6}, // (2, *)
    {3, 0}, {3, 1}, {3, 2}, {3, 5}, {3, 6}, {3, 7}, // (3, *)
    {4, 0}, {4, 5}, {4, 6}, // (4, *)
    {5, 0}, {5, 1}, {5, 3}, {5, 4}, {5, 6}, {5, 7}, // (5, *)
    {6, 0}, {6, 2}, {6, 3}, {6, 4}, {6, 5}, {6, 7}, // (6, *)
    {7, 3}, {7, 5}, {7, 6} // (7, *)
  };

  int nec[] = {
    1, 1, 3, 1, 3, 3, // (0, *)
    1, 1, 1, // (1, *)
    1, 1, 1, // (2, *)
    3, 1, 1, 3, 3, 1, // (3, *)
    1, 1, 1, // (4, *)
    3, 1, 3, 1, 3, 1, // (5, *)
    3, 1, 3, 1, 3, 1, // (6, *)
    1, 1, 1 // (7, *)
  };

  size_t ec_gt[][3] = {
    // (0, *)
    {3, NO_INDEX, NO_INDEX}, // 1
    {0, NO_INDEX, NO_INDEX}, // 2
    {0, 2, 3}, // 3
    {1, NO_INDEX, NO_INDEX}, // 4
    {1, 2, 3}, // 5
    {0, 1, 2}, // 6
    // (1, *)
    {3, NO_INDEX, NO_INDEX}, // 0
    {3, NO_INDEX, NO_INDEX}, // 3
    {3, NO_INDEX, NO_INDEX}, // 5
    // (2, *)
    {0, NO_INDEX, NO_INDEX}, // 0
    {0, NO_INDEX, NO_INDEX}, // 3
    {0, NO_INDEX, NO_INDEX}, // 6
    // (3, *)
    {0, 2, 3}, // 0
    {3, NO_INDEX, NO_INDEX}, // 1
    {0, NO_INDEX, NO_INDEX}, // 2
    {2, 3, 4}, // 5
    {0, 2, 4}, // 6
    {4, NO_INDEX, NO_INDEX}, // 7
    // (4, *)
    {1, NO_INDEX, NO_INDEX}, // 0
    {1, NO_INDEX, NO_INDEX}, // 5
    {1, NO_INDEX, NO_INDEX}, // 6
    // (5, *)
    {1, 2, 3}, // 0
    {3, NO_INDEX, NO_INDEX}, // 1
    {2, 3, 4}, // 3
    {1, NO_INDEX, NO_INDEX}, // 4
    {1, 2, 4}, // 6
    {4, NO_INDEX, NO_INDEX}, // 7
    // (6, *)
    {0, 1, 2}, // 0
    {0, NO_INDEX, NO_INDEX}, // 2
    {0, 2, 4}, // 3
    {1, NO_INDEX, NO_INDEX}, // 4
    {1, 2, 4}, // 5
    {4, NO_INDEX, NO_INDEX}, // 7
    // (7, *)
    {4, NO_INDEX, NO_INDEX}, // 3
    {4, NO_INDEX, NO_INDEX}, // 5
    {4, NO_INDEX, NO_INDEX}, // 7
  };

  int num_edges = sizeof(ec_gt)/sizeof(ec_gt[0]);

  size_t ec[3];

  for (int i = 0; i < num_edges; ++i) {
    mesh3_ec(mesh, edges[i][0], edges[i][1], ec);
    qsort(ec, nec[i], sizeof(size_t), (compar_t)compar_size_t);
    assert_that(ec, is_equal_to_contents_of(ec_gt[i], nec[i]*sizeof(size_t)));
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, nbde_works_for_cube) {
  SET_UP_CUBE_MESH();

  size_t nbde = mesh3_nbde(mesh);
  assert_that(nbde, is_equal_to(18));

  TEAR_DOWN_MESH();
}

Ensure (mesh3, nbdf_works_for_cube) {
  SET_UP_CUBE_MESH();

  size_t nbdf = mesh3_nbdf(mesh);
  assert_that(nbdf, is_equal_to(12));

  TEAR_DOWN_MESH();
}

Ensure (mesh3, bdc_works_for_cube) {
  SET_UP_CUBE_MESH();

  bool bdc_gt[5] = {true, true, false, true, true};

  for (int i = 0; i < 5; ++i) {
    bool bdc = mesh3_bdc(mesh, i);
    assert_that(bdc, is_equal_to(bdc_gt[i]));
  }

  TEAR_DOWN_MESH();
}

Ensure (mesh3, bdv_works_for_cube) {
  SET_UP_CUBE_MESH();

  bool bdv_gt[8] = {true, true, true, true, true, true, true, true};

  for (int i = 0; i < 8; ++i) {
    bool bdv = mesh3_bdv(mesh, i);
    assert_that(bdv, is_equal_to(bdv_gt[i]));
  }

  TEAR_DOWN_MESH();
}
