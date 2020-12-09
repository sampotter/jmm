#include "mesh3.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "index.h"

struct mesh3 {
  dvec3 *verts;
  size_t nverts;
  ind4 *cells;
  size_t ncells;
  size_t *vc;
  size_t *vc_offsets;
};

void mesh3_alloc(mesh3_s **mesh) {
  *mesh = malloc(sizeof(mesh3_s));
}

void mesh3_dealloc(mesh3_s **mesh) {
  free(*mesh);
  *mesh = NULL;
}

static void init_vc(mesh3_s *mesh) {
  // Allocate space to count the number of cells incident on each
  // vertex.
  int *nvc = calloc(mesh->nverts, sizeof(int));

  // Traverse the cells, incrementing the count of each incident
  // vertex.
  for (size_t i = 0; i < mesh->ncells; ++i) {
    for (int j = 0; j < 4; ++j) {
      size_t k = mesh->cells[i].data[j];
      assert(k >= 0);
      assert(k < mesh->nverts);
      ++nvc[k];
    }
  }

  // Compute the offsets into the array of vc's. The last entry
  // contains the length of the number of elements pointed to by vc.
  size_t *vc_offsets = malloc(sizeof(size_t)*(mesh->nverts + 1));
  vc_offsets[0] = 0;
  for (size_t i = 0; i < mesh->nverts; ++i) {
    vc_offsets[i + 1] = vc_offsets[i] + nvc[i];
  }

  // To avoid allocating a new block of memory, we zero out nvc here
  // and use it to store indices past the corresponding vc_offsets
  memset(nvc, 0x0, sizeof(int)*mesh->nverts);

  // Now that we've allocated some space and have the required offsets
  // into the array, we can traverse the cells again and fill the
  // array of vc's.
  size_t *vc = malloc(sizeof(size_t)*vc_offsets[mesh->nverts]);
  for (size_t i = 0; i < mesh->ncells; ++i) {
    for (int j = 0; j < 4; ++j) {
      size_t k = mesh->cells[i].data[j];
      vc[vc_offsets[k] + nvc[k]++] = i;
    }
  }

  mesh->vc = vc;
  mesh->vc_offsets = vc_offsets;

  free(nvc);
}

void mesh3_init(mesh3_s *mesh,
                dbl const *verts, size_t nverts,
                size_t const *cells, size_t ncells) {
  size_t k = 0;

  mesh->verts = malloc(sizeof(dvec3)*nverts);
  for (size_t i = 0; i < nverts; ++i) {
    for (int j = 0; j < 3; ++j) {
      mesh->verts[i].data[j] = verts[k++];
    }
  }
  mesh->nverts = nverts;

  k = 0;

  mesh->cells = malloc(sizeof(ind4)*ncells);
  for (size_t i = 0; i < ncells; ++i) {
    for (int j = 0; j < 4; ++j) {
      mesh->cells[i].data[j] = cells[k++];
    }
  }
  mesh->ncells = ncells;

  init_vc(mesh);
}

void mesh3_deinit(mesh3_s *mesh) {
  free(mesh->verts);
  free(mesh->cells);
  free(mesh->vc);
  free(mesh->vc_offsets);

  mesh->verts = NULL;
  mesh->cells = NULL;
  mesh->vc = NULL;
  mesh->vc_offsets = NULL;
}

void mesh3_get_vert(mesh3_s const *mesh, size_t i, dbl *v) {
  memcpy((void *)v, (void *)&mesh->verts[i], sizeof(dvec3));
}

size_t mesh3_nverts(mesh3_s const *mesh) {
  return mesh->nverts;
}

int mesh3_nvc(mesh3_s const *mesh, size_t i) {
  return mesh->vc_offsets[i + 1] - mesh->vc_offsets[i];
}

void mesh3_vc(mesh3_s const *mesh, size_t i, size_t *vc) {
  int nvc = mesh3_nvc(mesh, i);
  size_t *vci = &mesh->vc[mesh->vc_offsets[i]];
  memcpy((void *)vc, (void *)vci, sizeof(size_t)*nvc);
}

static bool contains(void *arr, int len, void *elt, size_t size) {
  char *ptr = (char *)arr;
  for (int i = 0; i < len; ++i) {
    if (!memcmp((void *)(ptr + size*i), elt, size)) {
      return true;
    }
  }
  return false;
}

/**
 * TODO: the functions below are simple, unoptimized
 * implementations. There are lots of ways to improve these, but we
 * want to wait on that until later.
 */

int mesh3_nvv(mesh3_s const *mesh, size_t i) {
  int nvc = mesh3_nvc(mesh, i);
  size_t *vc = malloc(sizeof(size_t)*nvc);
  mesh3_vc(mesh, i, vc);

  size_t *vv = malloc(3*sizeof(size_t)*nvc);

  int nvv = 0;
  for (int p = 0; p < nvc; ++p) {
    ind4 cell = mesh->cells[vc[p]];
    for (int q = 0; q < 4; ++q) {
      size_t j = cell.data[q];
      if (i == j || contains((void *)vv, nvv, &j, sizeof(size_t))) {
        continue;
      }
      vv[nvv++] = j;
    }
  }

  free(vv);
  free(vc);

  return nvv;
}

void mesh3_vv(mesh3_s const *mesh, size_t i, size_t *vv) {
  int nvc = mesh3_nvc(mesh, i);
  size_t *vc = malloc(sizeof(size_t)*nvc);
  mesh3_vc(mesh, i, vc);

  int k = 0;
  for (int p = 0; p < nvc; ++p) {
    ind4 cell = mesh->cells[vc[p]];
    for (int q = 0; q < 4; ++q) {
      size_t j = cell.data[q];
      if (i == j || contains((void *)vv, k, &j, sizeof(size_t))) {
        continue;
      }
      vv[k++] = j;
    }
  }

  free(vc);
}

int mesh3_nve(mesh3_s const *mesh, size_t i) {
  int nvc = mesh3_nvc(mesh, i);
  size_t *vc = malloc(sizeof(size_t)*nvc);
  mesh3_vc(mesh, i, vc);

  ind2 *ve = malloc(3*sizeof(ind2)*nvc);

  int nve = 0;
  for (int p = 0; p < nvc; ++p) {
    printf("p = %d\n", p);

    ind4 cell = mesh->cells[vc[p]];

    // Find the indices of cell vertices that aren't `i`.
    int r = 0, Q[4];
    for (int q = 0; q < 4; ++q) {
      if (cell.data[q] == i) continue;
      Q[r++] = q;
    }
    Q[r] = Q[0]; // Avoid having to use % below

    // Traverse distinct indices and add them to ve, incrementing nve
    // as we go.
    ind2 e;
    size_t j, k;
    for (int r = 0; r < 3; ++r) {
      j = cell.data[Q[r]];
      k = cell.data[Q[r + 1]];
      e.data[0] = j < k ? j : k;
      e.data[1] = j < k ? k : j;
      if (!contains(ve, nve, &e, sizeof(ind2))) {
        ve[nve++] = e;
      }
    }
  }

  free(ve);
  free(vc);

  return nve;
}

void mesh3_ve(mesh3_s const *mesh, size_t i, size_t *ve) {
  int nvc = mesh3_nvc(mesh, i);
  size_t *vc = malloc(sizeof(size_t)*nvc);
  mesh3_vc(mesh, i, vc);

  int nve = 0;
  for (int p = 0; p < nvc; ++p) {
    ind4 cell = mesh->cells[vc[p]];

    // Find the indices of cell vertices that aren't `i`.
    int r = 0, Q[4];
    for (int q = 0; q < 4; ++q) {
      if (cell.data[q] == i) continue;
      Q[r++] = q;
    }
    Q[r] = Q[0]; // Avoid having to use % below

    // Traverse distinct indices and add them to ve, incrementing nve
    // as we go.
    ind2 e;
    size_t j, k;
    for (int r = 0; r < 3; ++r) {
      j = cell.data[Q[r]];
      k = cell.data[Q[r + 1]];
      e.data[0] = j < k ? j : k;
      e.data[1] = j < k ? k : j;
      if (!contains(ve, nve, &e, sizeof(ind2))) {
        memcpy((void *)&ve[2*nve], &e.data[0], 2*sizeof(size_t));
        ++nve;
      }
    }
  }

  free(vc);
}

int mesh3_nvf(mesh3_s const *mesh, size_t i) {
  return mesh3_nvc(mesh, i);
}

void mesh3_vf(mesh3_s const *mesh, size_t i, size_t *vf) {
  int nvc = mesh3_nvc(mesh, i);
  size_t *vc = malloc(sizeof(size_t)*nvc);
  mesh3_vc(mesh, i, vc);

  int nvf = 0;
  for (int p = 0; p < nvc; ++p) {
    ind4 cell = mesh->cells[vc[p]];

    ind3 f;
    int r = 0;
    for (int q = 0; q < 4; ++q) {
      size_t j = cell.data[q];
      if (i == j) continue;
      f.data[r++] = j;
    }
    assert(r == 3);

    memcpy((void *)&vf[3*nvf], &f.data[0], 3*sizeof(size_t));
    ++nvf;
  }

  free(vc);
}

static int num_shared_verts(ind4 const *c1, ind4 const *c2) {
  // TODO: speed up using SIMD?
  int n = 0;
  for (int p = 0; p < 4; ++p) {
    for (int q = 0; q < 4; ++q) {
      if (c1->data[p] == c2->data[q]) {
        ++n;
        continue;
      }
    }
  }
  return n;
}

int mesh3_ncc(mesh3_s const *mesh, size_t i) {
  ind4 c = mesh->cells[i];

  int nvc[4], max_nvc = -1;
  for (int p = 0; p < 4; ++p) {
    nvc[p] = mesh3_nvc(mesh, c.data[p]);
    if (nvc[p] > max_nvc) max_nvc = nvc[p];
  }

  size_t *vc = malloc(sizeof(size_t)*max_nvc);

  size_t *cc = malloc(4*sizeof(size_t));

  int ncc = 0;

  size_t j, k;
  for (int p = 0; p < 4; ++p) {
    j = c.data[p];
    mesh3_vc(mesh, j, vc);
    for (int q = 0; q < nvc[p]; ++q) {
      k = vc[q];
      if (i != k
          && num_shared_verts(&c, &mesh->cells[k]) == 3
          && !contains(cc, ncc, &k, sizeof(size_t))) {
        cc[ncc++] = k;
      }
    }
  }

  free(cc);
  free(vc);

  return ncc;
}

void mesh3_cc(mesh3_s const *mesh, size_t i, size_t *cc) {
  ind4 c = mesh->cells[i];

  int nvc[4], max_nvc = -1;
  for (int p = 0; p < 4; ++p) {
    nvc[p] = mesh3_nvc(mesh, c.data[p]);
    if (nvc[p] > max_nvc) max_nvc = nvc[p];
  }

  size_t *vc = malloc(sizeof(size_t)*max_nvc);

  int ncc = 0;

  size_t j, k;
  for (int p = 0; p < 4; ++p) {
    j = c.data[p];
    mesh3_vc(mesh, j, vc);
    for (int q = 0; q < nvc[p]; ++q) {
      k = vc[q];
      if (i != k
          && num_shared_verts(&c, &mesh->cells[k]) == 3
          && !contains(cc, ncc, &k, sizeof(size_t))) {
        cc[ncc++] = k;
      }
    }
  }

  free(vc);
}

void mesh3_cv(mesh3_s const *mesh, size_t i, size_t *cv) {
  memcpy((void *)cv, (void *)&mesh->cells[i].data, 4*sizeof(size_t));
}

int mesh3_nec(mesh3_s const *mesh, size_t i, size_t j) {
  // TODO: really horrible implementation! :-(

  int nvci = mesh3_nvc(mesh, i);
  size_t *vci = malloc(sizeof(size_t)*nvci);
  mesh3_vc(mesh, i, vci);

  int nvcj = mesh3_nvc(mesh, j);
  size_t *vcj = malloc(sizeof(size_t)*nvcj);
  mesh3_vc(mesh, j, vcj);

  int nec = 0;

  size_t c;
  for (int a = 0; a < nvci; ++a) {
    c = vci[a];
    for (int b = 0; b < nvcj; ++b) {
      if (c == vcj[b]) {
        ++nec;
        continue;
      }
    }
  }

  free(vci);
  free(vcj);

  return nec;
}

void mesh3_ec(mesh3_s const *mesh, size_t i, size_t j, size_t *ec) {
  // TODO: really horrible implementation! :-(

  int nvci = mesh3_nvc(mesh, i);
  size_t *vci = malloc(sizeof(size_t)*nvci);
  mesh3_vc(mesh, i, vci);

  int nvcj = mesh3_nvc(mesh, j);
  size_t *vcj = malloc(sizeof(size_t)*nvcj);
  mesh3_vc(mesh, j, vcj);

  int nec = 0;

  size_t c;
  for (int a = 0; a < nvci; ++a) {
    c = vci[a];
    for (int b = 0; b < nvcj; ++b) {
      if (c == vcj[b]) {
        ec[nec++] = c;
        continue;
      }
    }
  }

  free(vci);
  free(vcj);
}

bool mesh3_bdc(mesh3_s const *mesh, size_t i) {
  return mesh3_ncc(mesh, i) < 4;
}
