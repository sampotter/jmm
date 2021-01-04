#include "mesh3.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "index.h"
#include "util.h"

struct mesh3 {
  dvec3 *verts;
  size_t nverts;
  ind4 *cells;
  size_t ncells;
  size_t *vc;
  size_t *vc_offsets;
  bool *bdc;
  bool *bdv;
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

typedef struct {
  size_t lf[3];
  size_t lc;
} tagged_face_s;

int compar_tagged_face(tagged_face_s const *f1, tagged_face_s const *f2) {
  for (int i = 0, cmp; i < 3; ++i) {
    cmp = compar_size_t(&f1->lf[i], &f2->lf[i]);
    if (cmp != 0) {
      return cmp;
    }
  }
  return 0;
}

/**
 * In this function we figure out which cells (tetrahedra) and
 * vertices are on the boundary. This is slightly arbitrary. We
 * stipulate that a vertex is on the boundary if every ball
 * surrounding the vertex intersects the exterior of the domain. A
 * cell is a boundary cell if it has a face that's incident on the
 * boundary of the domain.
 *
 * Using this information we find the unique faces in the mesh. A
 */
static void init_bd(mesh3_s *mesh) {
  mesh->bdc = calloc(mesh->ncells, sizeof(bool));
  mesh->bdv = calloc(mesh->nverts, sizeof(bool));

  // Allocate some space for the faces of each cell in the mesh.
  size_t nf = 4*mesh->ncells;
  tagged_face_s *f = malloc(nf*sizeof(tagged_face_s));

  // Traverse the cells in the mesh, and populate `f`. These faces are
  // "tagged", meaning that they have a backpointer to the originating
  // cell.
  size_t *C;
  for (size_t lc = 0; lc < mesh->ncells; ++lc) {
    C = &mesh->cells[lc].data[0];
    f[4*lc] = (tagged_face_s) {.lf = {C[0], C[1], C[2]}, .lc = lc};
    f[4*lc + 1] = (tagged_face_s) {.lf = {C[0], C[1], C[3]}, .lc = lc};
    f[4*lc + 2] = (tagged_face_s) {.lf = {C[0], C[2], C[3]}, .lc = lc};
    f[4*lc + 3] = (tagged_face_s) {.lf = {C[1], C[2], C[3]}, .lc = lc};
  }

  // Sort the components of each tagged face.
  for (size_t l = 0; l < nf; ++l) {
    qsort(f[l].lf, 3, sizeof(size_t), (compar_t)compar_size_t);
  }

  // Sort the tagged faces themselves into a dictionary order.
  qsort(f, nf, sizeof(tagged_face_s), (compar_t)compar_tagged_face);

  // After sorting, we can tell if a face is duplicated by checking
  // whether it's equal to the succeeding face in `f`. If there's a
  // duplicate, we update the relevant boundary information.
  for (size_t l = 0; l < nf - 1; ++l) {
    if (!compar_tagged_face(&f[l], &f[l + 1])) {
      ++l;
      continue;
    }
    mesh->bdc[f[l].lc] = true;
    mesh->bdv[f[l].lf[0]] = true;
    mesh->bdv[f[l].lf[1]] = true;
    mesh->bdv[f[l].lf[2]] = true;
  }

  free(f);
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
  init_bd(mesh);
}

void mesh3_deinit(mesh3_s *mesh) {
  free(mesh->verts);
  free(mesh->cells);
  free(mesh->vc);
  free(mesh->vc_offsets);
  free(mesh->bdc);
  free(mesh->bdv);

  mesh->verts = NULL;
  mesh->cells = NULL;
  mesh->vc = NULL;
  mesh->vc_offsets = NULL;
  mesh->bdc = NULL;
  mesh->bdv = NULL;
}

dvec3 mesh3_get_vert(mesh3_s const *mesh, size_t i) {
  return mesh->verts[i];
}

dbl const *mesh3_get_vert_ptr(mesh3_s const *mesh, size_t i) {
  return &mesh->verts[i].data[0];
}

void mesh3_copy_vert(mesh3_s const *mesh, size_t i, dbl *v) {
  memcpy(v, &mesh->verts[i].data[0], 3*sizeof(dbl));
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

int mesh3_nvf(mesh3_s const *mesh, size_t l) {
  return mesh3_nvc(mesh, l);
}

void mesh3_vf(mesh3_s const *mesh, size_t l, size_t (*vf)[3]) {
  // TODO: this is a very bad implementation!

  int nvc = mesh3_nvc(mesh, l);
  size_t *vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(mesh, l, vc);

  size_t cv[4];

  int k;
  for (int i = 0; i < nvc; ++i) {
    mesh3_cv(mesh, vc[i], cv);

    k = 0;
    for (int j = 0; j < 4; ++j) {
      if (cv[j] == l) {
        continue;
      }
      vf[i][k++] = cv[j];
    }
  }

  free(vc);
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

  if (i == j) {
    return 0;
  }

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

  if (i == j) {
    return;
  }

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
        continue; // TODO: should be break? test later...
      }
    }
  }

  free(vci);
  free(vcj);
}

bool *mesh3_get_bdc_ptr(mesh3_s *mesh) {
  return mesh->bdc;
}

bool mesh3_bdc(mesh3_s const *mesh, size_t i) {
  return mesh->bdc[i];
}

bool *mesh3_get_bdv_ptr(mesh3_s *mesh) {
  return mesh->bdv;
}

bool mesh3_bdv(mesh3_s const *mesh, size_t i) {
  return mesh->bdv[i];
}
