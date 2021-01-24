#include "mesh3.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "index.h"
#include "macros.h"
#include "mat.h"
#include "util.h"

bool point_in_face(size_t l, size_t const f[3]) {
  return f[0] == l || f[1] == l || f[2] == l;
}

bool point_in_cell(size_t l, size_t const c[4]) {
  return c[0] == l || c[1] == l || c[2] == l || c[3] == l;
}

typedef struct {
  size_t le[2];
  bool diff;
} diff_edge_s;

diff_edge_s make_diff_edge(size_t l0, size_t l1) {
  return (diff_edge_s) {.le = {MIN(l0, l1), MAX(l0, l1)}, .diff = false};
}

int diff_edge_cmp(diff_edge_s const *e1, diff_edge_s const *e2) {
  int cmp = compar_size_t(&e1->le[0], &e2->le[0]);
  if (cmp != 0) {
    return cmp;
  } else {
    return compar_size_t(&e1->le[1], &e2->le[1]);
  }
}

typedef struct {
  size_t lf[3];
  size_t lc;
} tagged_face_s;

tagged_face_s make_tagged_face(size_t l0, size_t l1, size_t l2, size_t lc) {
  tagged_face_s f = {.lf = {l0, l1, l2}, .lc = lc};
  qsort(f.lf, 3, sizeof(size_t), (compar_t)compar_size_t);
  return f;
}

void tagged_face_init(tagged_face_s *f, size_t const *lf, size_t lc) {
  memcpy(f->lf, lf, 3*sizeof(size_t));
  qsort(f->lf, 3, sizeof(size_t), (compar_t)compar_size_t);
  f->lc = lc;
}

int tagged_face_cmp(tagged_face_s const *f1, tagged_face_s const *f2) {
  for (int i = 0, cmp; i < 3; ++i) {
    cmp = compar_size_t(&f1->lf[i], &f2->lf[i]);
    if (cmp != 0) {
      return cmp;
    }
  }
  return 0;
}

struct mesh3 {
  dvec3 *verts;
  size_t nverts;
  ind4 *cells;
  size_t ncells;
  size_t *vc;
  size_t *vc_offsets;
  bool *bdc;
  bool *bdv;
  size_t nbdf;
  tagged_face_s *bdf;
  size_t nbde;
  diff_edge_s *bde;
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

static void get_op_edge(mesh3_s const *mesh, size_t lc, diff_edge_s const *e,
                        size_t l[2]) {
  size_t lv[4];
  mesh3_cv(mesh, lc, lv);
  for (int i = 0, j = 0; i < 4; ++i) {
    if (lv[i] == e->le[0] || lv[i] == e->le[1]) continue;
    l[j++] = lv[i];
  }
}

static dbl get_dihedral_angle(mesh3_s const *mesh, size_t lc, diff_edge_s const *e) {
  // Get edge endpoints
  dbl const *x0 = mesh3_get_vert_ptr(mesh, e->le[0]);
  dbl const *x1 = mesh3_get_vert_ptr(mesh, e->le[1]);

  // Compute normalize direction vector along edge
  dbl t[3];
  dbl3_sub(x1, x0, t);
  dbl3_normalize(t);

  // Get indices of opposite edge for tetrahedron lc
  size_t l[2];
  get_op_edge(mesh, lc, e, l);

  // Get corresponding vertices
  dbl const *x2 = mesh3_get_vert_ptr(mesh, l[0]);
  dbl const *x3 = mesh3_get_vert_ptr(mesh, l[1]);

  dbl tmp[3], s, xs[3], t2[3], t3[3];

  dbl3_sub(x2, x0, tmp);
  s = dbl3_dot(t, tmp);
  dbl3_saxpy(s, t, x0, xs);
  dbl3_sub(x2, xs, t2);
  dbl3_normalize(t2);

  dbl3_sub(x3, x0, tmp);
  s = dbl3_dot(t, tmp);
  dbl3_saxpy(s, t, x0, xs);
  dbl3_sub(x3, xs, t3);
  dbl3_normalize(t3);

  dbl t2_dot_t3 = dbl3_dot(t2, t3);

  return acos(t2_dot_t3);
}

/**
 * We want to check and see if this edge is a diffracting edge.  An
 * edge is diffracting if its interior dihedral angle is greater than
 * 180 degrees. So, we traverse the tetrahedra surrounding it, and sum of the angles they make with
 */
static bool edge_is_diff(mesh3_s const *mesh, diff_edge_s *e) {
  dbl const atol = 1e-14;

  int nec = mesh3_nec(mesh, e->le[0], e->le[1]);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, e->le[0], e->le[1], ec);

  dbl angle_sum = 0;
  for (int i = 0; i < nec; ++i) {
    angle_sum += get_dihedral_angle(mesh, ec[i], e);
  }

  free(ec);

  return angle_sum > PI + atol;
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
    f[4*lc] = make_tagged_face(C[0], C[1], C[2], lc);
    f[4*lc + 1] = make_tagged_face(C[0], C[1], C[3], lc);
    f[4*lc + 2] = make_tagged_face(C[0], C[2], C[3], lc);
    f[4*lc + 3] = make_tagged_face(C[1], C[2], C[3], lc);
  }

  // Sort the tagged faces themselves into a dictionary order.
  qsort(f, nf, sizeof(tagged_face_s), (compar_t)tagged_face_cmp);

  mesh->nbdf = 0;

  // After sorting, we can tell if a face is duplicated by checking
  // whether it's equal to the succeeding face in `f`. If there's a
  // duplicate, we update the relevant boundary information.
  for (size_t l = 0; l < nf - 1; ++l) {
    if (!tagged_face_cmp(&f[l], &f[l + 1])) {
      ++l;
      continue;
    }
    mesh->bdc[f[l].lc] = true;
    mesh->bdv[f[l].lf[0]] = true;
    mesh->bdv[f[l].lf[1]] = true;
    mesh->bdv[f[l].lf[2]] = true;
    ++mesh->nbdf;
  }

  // Traverse the sorted list again and pull out the boundary faces
  // now that we know how many there are
  mesh->bdf = malloc(mesh->nbdf*sizeof(tagged_face_s));
  for (size_t l = 0, lf = 0; l < nf - 1; ++l) {
    if (!tagged_face_cmp(&f[l], &f[l + 1])) {
      ++l; // Increment here to skip equal pairs
      continue;
    }
    tagged_face_init(&mesh->bdf[lf++], f[l].lf, f[l].lc);
  }

  // Sort the faces so that we can quickly query whether a face is a
  // boundary face or not.
  qsort(mesh->bdf, mesh->nbdf, sizeof(tagged_face_s), (compar_t)tagged_face_cmp);

  // Build a sorted array of all of the boundary edges, which are just
  // the edges incident on the boundary faces. This array will have
  // duplicates, so we'll have to prune these next.
  diff_edge_s *bde = malloc(3*mesh->nbdf*sizeof(diff_edge_s));
  for (size_t lf = 0, *l; lf < mesh->nbdf; ++lf) {
    l = mesh->bdf[lf].lf;
    bde[3*lf] = make_diff_edge(l[0], l[1]);
    bde[3*lf + 1] = make_diff_edge(l[1], l[2]);
    bde[3*lf + 2] = make_diff_edge(l[2], l[0]);
  }
  qsort(bde, 3*mesh->nbdf, sizeof(diff_edge_s), (compar_t)diff_edge_cmp);

  // Now, let's count the number of distinct boundary edges.
  mesh->nbde = 0;
  for (size_t l = 0; l < 3*mesh->nbdf - 1; ++l) {
    if  (!diff_edge_cmp(&bde[l], &bde[l + 1]))
      continue;
    ++mesh->nbde;
  }

  // Traverse the array again, copying over distinct boundary
  // edges. Note: there's no need to sort mesh->bde, since it will
  // already be sorted.
  mesh->bde = malloc(mesh->nbde*sizeof(diff_edge_s));
  for (size_t l = 0, l_ = 0; l < 3*mesh->nbdf - 1; ++l) {
    if (!diff_edge_cmp(&bde[l], &bde[l + 1]))
      continue;
    mesh->bde[l_++] = bde[l];
  }

  for (size_t l = 0; l < mesh->nbde; ++l) {
    mesh->bde[l].diff = edge_is_diff(mesh, &mesh->bde[l]);
  }

  free(bde);
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

size_t mesh3_ncells(mesh3_s const *mesh) {
  return mesh->ncells;
}

size_t mesh3_nverts(mesh3_s const *mesh) {
  return mesh->nverts;
}

void mesh3_get_bbox(mesh3_s const *mesh, rect3 *bbox) {
  dbl *min = bbox->min, *max = bbox->max, *v;

  min[0] = min[1] = min[2] = INFINITY;
  for (size_t l = 0; l < mesh->nverts; ++l) {
    v = &mesh->verts[l].data[0];
    min[0] = fmin(min[0], v[0]);
    min[1] = fmin(min[1], v[1]);
    min[2] = fmin(min[2], v[2]);
  }

  max[0] = max[1] = max[2] = -INFINITY;
  for (size_t l = 0; l < mesh->nverts; ++l) {
    v = &mesh->verts[l].data[0];
    max[0] = fmax(max[0], v[0]);
    max[1] = fmax(max[1], v[1]);
    max[2] = fmax(max[2], v[2]);
  }
}

void mesh3_get_cell_bbox(mesh3_s const *mesh, size_t i, rect3 *bbox) {
  size_t const *cell = &mesh->cells[i].data[0];

  dbl *min = bbox->min, *max = bbox->max;

  min[0] = min[1] = min[2] = INFINITY;
  for (int i = 0; i < 4; ++i) {
    min[0] = fmin(min[0], mesh->verts[cell[i]].data[0]);
    min[1] = fmin(min[1], mesh->verts[cell[i]].data[1]);
    min[2] = fmin(min[2], mesh->verts[cell[i]].data[2]);
  }

  max[0] = max[1] = max[2] = -INFINITY;
  for (int i = 0; i < 4; ++i) {
    max[0] = fmax(max[0], mesh->verts[cell[i]].data[0]);
    max[1] = fmax(max[1], mesh->verts[cell[i]].data[1]);
    max[2] = fmax(max[2], mesh->verts[cell[i]].data[2]);
  }
}

bool mesh3_dbl3_in_cell(mesh3_s const *mesh, size_t lc, dbl const x[3],
                        dbl b[4]) {
  dbl const atol = 1e-13;

  assert(lc < mesh->ncells);

  size_t const *c = &mesh->cells[lc].data[0];

  dbl lhs[4][4];
  for (int j = 0; j < 4; ++j) {
    lhs[0][j] = 1;
    for (int i = 1; i < 4; ++i) {
      lhs[i][j] = mesh->verts[c[j]].data[i - 1];
    }
  }

  dbl vol = dbl44_det(lhs);

  dbl rhs[4];
  rhs[0] = 1;
  for (int i = 1; i < 4; ++i) {
    rhs[i] = x[i - 1];
  }

  dbl tmp[4];
  for (int j = 0; j < 4; ++j) {
    dbl44_get_col(lhs, j, tmp);
    dbl44_set_col(lhs, j, rhs);
    b[j] = dbl44_det(lhs)/vol;
    dbl44_set_col(lhs, j, tmp);
  }

  // After solving this system, the components of b will sum to unity,
  // so we just need to check if b is nonnegative.
  return b[0] >= -atol && b[1] >= -atol && b[2] >= -atol && b[3] >= -atol;
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

bool mesh3_cee(mesh3_s const *mesh, size_t c, size_t const e[2],
               size_t e_out[2]) {
  size_t cv[4];
  mesh3_cv(mesh, c, cv);
  int k = 0;
  for (int i = 0; i < 4; ++i) {
    if (cv[i] == e[0] || cv[i] == e[1])
      continue;
    if (k < 2)
      e_out[k++] = cv[i];
    else
      return false;
  }
  return true;
}

/**
 * This function takes two chains of vertices corresponding to the
 * edges of a loop and orients them so that one can traverse the loop
 * in order.
 */
static void orient_chains(size_t (*e)[2], int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      // We don't want any duplicate entries among the l0 or l1
      // indices, but this can easily happen. This is the "orient"
      // part of this function.
      if (e[j][0] == e[i][0] || e[j][1] == e[i][1]) {
        SWAP(e[j][0], e[j][1]);
      }
      if (e[j][1] == e[i][0]) {
        SWAP(e[i + 1][0], e[j][0]);
        SWAP(e[i + 1][1], e[j][1]);
      }
    }
  }
}

int mesh3_nee(mesh3_s const *mesh, size_t const e[2]) {
  return mesh3_nec(mesh, e[0], e[1]);
}

void mesh3_ee(mesh3_s const *mesh, size_t const e[2], size_t (*ee)[2]) {
  int nec = mesh3_nec(mesh, e[0], e[1]);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, e[0], e[1], ec);

  for (int i = 0; i < nec; ++i)
    mesh3_cee(mesh, ec[i], e, ee[i]);
  orient_chains(ee, nec);

  free(ec);
}

bool mesh3_cfv(mesh3_s const *mesh, size_t lc, size_t const lf[3], size_t *lv) {
  size_t cv[4];
  mesh3_cv(mesh, lc, cv);
  // First, check if the verts in lf actually belong to cv. If they
  // don't, return false, since cfv no longer makes any sense.
  for (int i = 0; i < 3; ++i)
    if (!point_in_cell(lf[i], cv))
      return false;
  // Now, figure out which vertex isn't one of cv isn't one of lf, and
  // set *lv to that vertex.
  *lv = NO_PARENT;
  for (int i = 0; i < 4; ++i) {
    if (cv[i] == lf[0] || cv[i] == lf[1] || cv[i] == lf[2])
      continue;
    *lv = cv[i];
    break;
  }
  assert(*lv != NO_PARENT);
  return true;
}

bool mesh3_ccfv(mesh3_s const *mesh, size_t lc, size_t const lf[3],
                size_t *lv_out) {
  int ncc = mesh3_ncc(mesh, lc);
  size_t *cc = malloc(ncc*sizeof(size_t));
  mesh3_cc(mesh, lc, cc);
  for (int i = 0; i < ncc; ++i)
    if (mesh3_cfv(mesh, cc[i], lf, lv_out))
      return true;
  return false;
}

bool mesh3_cvf(mesh3_s const *mesh, size_t lc, size_t lv, size_t lf[3]) {
  size_t cv[4];
  mesh3_cv(mesh, lc, cv);
  // First, check if lv is actually incident on lc
  if (!point_in_cell(lv, cv))
    return false;
  // Now, copy the other vertices over to lf
  int j = 0;
  for (int i = 0; i < 4; ++i)
    if (cv[i] != lv)
      lf[j++] = cv[i];
  assert(j == 3);
  return true;
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

bool mesh3_bde(mesh3_s const *mesh, size_t const le[2]) {
  diff_edge_s e = make_diff_edge(le[0], le[1]);
  return bsearch(&e, mesh->bde, mesh->nbde, sizeof(diff_edge_s),
                 (compar_t)diff_edge_cmp);
}

bool mesh3_bdf(mesh3_s const *mesh, size_t const lf[3]) {
  tagged_face_s f = make_tagged_face(lf[0], lf[1], lf[2], NO_PARENT);
  return bsearch(&f, mesh->bdf, mesh->nbdf, sizeof(tagged_face_s),
                 (compar_t)tagged_face_cmp);
}

bool mesh3_is_diff_edge(mesh3_s const *mesh, size_t const le[2]) {
  diff_edge_s q = make_diff_edge(le[0], le[1]);
  diff_edge_s const *e = bsearch(
    &q, mesh->bde, mesh->nbde, sizeof(diff_edge_s), (compar_t)diff_edge_cmp);
  return e && e->diff;
}

bool mesh3_vert_incident_on_diff_edge(mesh3_s const *mesh, size_t l) {
  int nvv = mesh3_nvv(mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l, vv);

  bool is_incident = false;
  for (int i = 0; i < nvv; ++i)
    if ((is_incident = mesh3_is_diff_edge(mesh, (size_t[2]) {l, vv[i]})))
      break;

  free(vv);

  return is_incident;
}
