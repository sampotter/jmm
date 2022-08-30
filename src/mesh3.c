#include <jmm/mesh3.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <jmm/array.h>
#include <jmm/edge.h>
#include <jmm/index.h>
#include <jmm/mat.h>
#include <jmm/mesh1.h>
#include <jmm/mesh2.h>
#include <jmm/util.h>

#include "macros.h"
#include "mesh_util.h"

typedef struct {
  size_t le[2];
  bool diff;
} bde_s;

bde_s make_bde(size_t l0, size_t l1) {
  return (bde_s) {.le = {MIN(l0, l1), MAX(l0, l1)}, .diff = false};
}

int bde_cmp(bde_s const *e1, bde_s const *e2) {
  return edge_cmp(e1->le, e2->le);
}

typedef struct {
  size_t lf[3];
  size_t lc;
} bdf_s;

bdf_s make_bdf(size_t l0, size_t l1, size_t l2, size_t lc) {
  bdf_s f = {.lf = {l0, l1, l2}, .lc = lc};
  qsort(f.lf, 3, sizeof(size_t), (compar_t)compar_size_t);
  return f;
}

void bdf_init(bdf_s *f, size_t const *lf, size_t lc) {
  memcpy(f->lf, lf, 3*sizeof(size_t));
  qsort(f->lf, 3, sizeof(size_t), (compar_t)compar_size_t);
  f->lc = lc;
}

int bdf_cmp(bdf_s const *f1, bdf_s const *f2) {
  for (int i = 0, cmp; i < 3; ++i) {
    cmp = compar_size_t(&f1->lf[i], &f2->lf[i]);
    if (cmp != 0)
      return cmp;
  }
  return 0;
}

struct mesh3 {
  size_t nverts;
  dbl3 *verts;

  size_t ncells;
  uint4 *cells;

  size_t *vc;
  size_t *vc_offsets;

  size_t (*edges)[2];
  size_t nedges;

  bool has_bd_info;
  bool *bdc;
  bool *bdv;
  size_t nbdf;
  bdf_s *bdf;
  size_t nbde;
  bde_s *bde;

  size_t num_bdf_labels;
  size_t *bdf_label; // Labels for distinct reflecting surfaces

  size_t num_bde_labels;
  size_t *bde_label;

  /* "Mesh epsilon": a small parameter derived from the mesh, used to
   * make geometric calculations a bit more robust. */
  dbl eps;

  // Geometric quantities
  dbl min_tetra_alt; // The minimum of all tetrahedron altitudes in the mesh.
  dbl min_edge_length;
  dbl mean_edge_length;
  dbl diam;
};

tri3 mesh3_tetra_get_face(mesh3_tetra_s const *tetra, int f[3]) {
  size_t *cv = tetra->mesh->cells[tetra->l];
  dbl (*verts)[3] = tetra->mesh->verts;
  tri3 tri;
  memcpy(tri.v[0], verts[cv[f[0]]], sizeof(dbl[3]));
  memcpy(tri.v[1], verts[cv[f[1]]], sizeof(dbl[3]));
  memcpy(tri.v[2], verts[cv[f[2]]], sizeof(dbl[3]));
  return tri;
}

bool mesh3_tetra_equal(mesh3_tetra_s const *t1, mesh3_tetra_s const *t2) {
  return t1->mesh == t2->mesh && t1->l == t2->l;
}

size_t mesh3_data_append_vert(mesh3_data_s *data, dbl3 const x) {
  data->verts = reallocarray(data->verts, data->nverts + 1, sizeof(dbl3));
  dbl3_copy(x, data->verts[data->nverts]);
  return data->nverts++;
}

void mesh3_data_append_cells(mesh3_data_s *data, size_t n, uint4 const C[4]) {
  data->cells = reallocarray(data->cells, data->ncells + n, sizeof(uint4));
  memcpy(data->cells[data->ncells], C, n*sizeof(uint4));
  data->ncells += n;
}

void mesh3_data_delete_cell(mesh3_data_s *data, size_t lc0) {
  --data->ncells;
  for (size_t lc = lc0; lc < data->ncells; ++lc)
    memcpy(data->cells[lc], data->cells[lc + 1], sizeof(uint4));
  data->cells = reallocarray(data->cells, data->ncells, sizeof(uint4));
}

bool mesh3_data_cell_contains_point(mesh3_data_s const *data, size_t lc,
                                    dbl3 const x, dbl eps) {
  tetra3 tetra;
  for (size_t i = 0; i < 4; ++i)
    dbl3_copy(data->verts[data->cells[lc][i]], tetra.v[i]);
  return tetra3_contains_point(&tetra, x, &eps);
}

size_t mesh3_data_find_cell_containing_point(mesh3_data_s const *data, dbl3 const x, dbl eps) {
  for (size_t lc = 0; lc < data->ncells; ++lc)
    if (mesh3_data_cell_contains_point(data, lc, x, eps))
      return lc;
  return (size_t)NO_INDEX;
}

error_e mesh3_data_insert_vert(mesh3_data_s *data, dbl3 const x, dbl eps) {
  size_t lc = mesh3_data_find_cell_containing_point(data, x, eps);
  if (lc == (size_t)NO_INDEX)
    return BAD_ARGUMENT;

  size_t lv = mesh3_data_append_vert(data, x);

  uint4 cells[4];
  for (size_t i = 0; i < 4; ++i) {
    memcpy(cells[i], data->cells[lc], sizeof(uint4));
    cells[i][i] = lv;
  }

  mesh3_data_delete_cell(data, lc);
  mesh3_data_append_cells(data, 4, cells);

  return SUCCESS;
}

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
      size_t k = mesh->cells[i][j];
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
      size_t k = mesh->cells[i][j];
      vc[vc_offsets[k] + nvc[k]++] = i;
    }
  }

  mesh->vc = vc;
  mesh->vc_offsets = vc_offsets;

  free(nvc);
}

static void init_edges(mesh3_s *mesh) {
  array_s *edge_arr;
  array_alloc(&edge_arr);
  array_init(edge_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);

  /* Initially accumulate all the cell edges into one big array */
  size_t ie[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  for (size_t lc = 0; lc < mesh->ncells; ++lc) {
    size_t const *cell = mesh->cells[lc];
    for (size_t i = 0; i < 6; ++i) {
      size_t le[2] = {cell[ie[i][0]], cell[ie[i][1]]};
      SORT2(le[0], le[1]);
      array_append(edge_arr, &le);
    }
  }

  /* Sort the array */
  array_sort(edge_arr, (compar_t)edge_cmp);

  array_s *unique_edge_arr;
  array_alloc(&unique_edge_arr);
  array_init(unique_edge_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);

  /* Pull out the unique edges into a new array */
  size_t n = array_size(edge_arr);
  for (size_t i = 0, le[2]; i < n; ) {
    array_get(edge_arr, i++, &le);
    array_append(unique_edge_arr, &le);
    while (i < n && !edge_cmp(le, array_get_ptr(edge_arr, i)))
      ++i;
  }

  /* Copy the edges over */
  mesh->nedges = array_size(unique_edge_arr);
  size_t nbytes = mesh->nedges*sizeof(size_t[2]);
  mesh->edges = malloc(nbytes);
  memcpy(mesh->edges, array_get_ptr(unique_edge_arr, 0), nbytes);

  /* Sanity check */
#if JMM_DEBUG
  for (size_t i = 1; i < mesh->nedges; ++i)
    assert(edge_cmp(mesh->edges[i - 1], mesh->edges[i]));
#endif

  array_deinit(unique_edge_arr);
  array_dealloc(&unique_edge_arr);

  array_deinit(edge_arr);
  array_dealloc(&edge_arr);
}

static void get_op_edge(mesh3_s const *mesh, size_t lc, size_t const le[2],
                        size_t l[2]) {
  size_t lv[4];
  mesh3_cv(mesh, lc, lv);
  for (int i = 0, j = 0; i < 4; ++i) {
    if (lv[i] == le[0] || lv[i] == le[1]) continue;
    l[j++] = lv[i];
  }
}

static dbl get_dihedral_angle(mesh3_s const *mesh, size_t lc, size_t const le[2]) {
  // Get edge endpoints
  dbl const *x0 = mesh3_get_vert_ptr(mesh, le[0]);
  dbl const *x1 = mesh3_get_vert_ptr(mesh, le[1]);

  // Compute normalize direction vector along edge
  dbl t[3];
  dbl3_sub(x1, x0, t);
  dbl3_normalize(t);

  // Get indices of opposite edge for tetrahedron lc
  size_t l[2];
  get_op_edge(mesh, lc, le, l);

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
 * 180 degrees. So, we traverse the tetrahedra surrounding it, and sum
 * the corresponding dihedral angles.
 */
static bool edge_is_diff(mesh3_s const *mesh, size_t const le[2]) {
  int nec = mesh3_nec(mesh, le);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, le, ec);

  dbl angle_sum = 0;
  for (int i = 0; i < nec; ++i)
    angle_sum += get_dihedral_angle(mesh, ec[i], le);

  free(ec);

  return angle_sum > JMM_PI + mesh->eps;
}

static size_t find_bdf(mesh3_s const *mesh, bdf_s const *bdf) {
  bdf_s const *found = bsearch(
    bdf, mesh->bdf, mesh->nbdf, sizeof(bdf_s), (compar_t)bdf_cmp);
  return (size_t)(found ? found - mesh->bdf : NO_INDEX);
}

/* Find all of the boundary faces that are adjacent to `mesh->bdf[l]`
 * and append them to `nb`. This doesn't assume that `nb` is empty. */
static void get_bdf_nbs(mesh3_s const *mesh, size_t l, array_s *nb) {
  bdf_s const *bdf = &mesh->bdf[l];
  bdf_s nb_bdf;

  for (size_t i = 0, lv, nve, (*ve)[2]; i < 3; ++i) {
    lv = bdf->lf[i];

    nve = mesh3_nve(mesh, lv);
    ve = malloc(nve*sizeof(size_t[2]));
    mesh3_ve(mesh, lv, ve);

    for (size_t j = 0, l_nb; j < nve; ++j) {
      nb_bdf = make_bdf(lv, ve[j][0], ve[j][1], (size_t)NO_INDEX);
      if (bdf_cmp(bdf, &nb_bdf) == 0)
        continue;
      l_nb = find_bdf(mesh, &nb_bdf);
      if (l_nb != (size_t)NO_INDEX && !array_contains(nb, &l_nb))
        array_append(nb, &l_nb);
    }

    free(ve);
  }
}

static bool bdfs_are_coplanar(mesh3_s const *mesh, size_t l0, size_t l1) {
  tri3 const tri0 = mesh3_get_tri(mesh, mesh->bdf[l0].lf);
  tri3 const tri1 = mesh3_get_tri(mesh, mesh->bdf[l1].lf);
  return tri3_coplanar(&tri0, &tri1, &mesh->eps);
}

static bool label_reflector(mesh3_s *mesh) {
  bool done = false;

  /* Set up a queue for the breadth-first search */
  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Set up a queue for streaming the neighbors of each face */
  array_s *nb;
  array_alloc(&nb);
  array_init(nb, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Find the first unlabeled face. */
  size_t lf = 0;
  while (lf < mesh->nbdf
         && mesh->bdf_label[lf] != NO_LABEL)
    ++lf;

  /* If we couldn't find an unlabeled face, we're done, so we can
   * return `false`, signaling no more reflectors to process. */
  if (lf == mesh->nbdf) {
    done = true;
    goto cleanup;
  }

  /* Push `lf` onto `queue` to start the BFS. */
  array_append(queue, &lf);

  /* Traverse the faces of the current reflector and label them. */
  do {
    array_pop_front(queue, &lf);

    /* Set the label of the face to the current label. */
    mesh->bdf_label[lf] = mesh->num_bdf_labels;

    /* Get `lf`'s neighbors, and store them in `nb`. */
    get_bdf_nbs(mesh, lf, nb);

    /* Traverse the neighboring boundary faces. If the neighboring
     * faces are coplanar, append them to `queue`. */
    size_t lf_nb;
    while (!array_is_empty(nb)) {
      array_pop_front(nb, &lf_nb);

      /* If face `lf_nb` is already labeled, skip it. */
      if (mesh->bdf_label[lf_nb] != NO_LABEL)
        continue;

      /* If `lf_nb` and `lf` are coplanar, try to add it to
       * `queue`. Either way, continue the iteration. We only want
       * to set `lf_next = `lf_nb` if they aren't coplanar. */
      if (bdfs_are_coplanar(mesh, lf, lf_nb)
          && !array_contains(queue, &lf_nb))
        array_append(queue, &lf_nb);
    }
  } while (!array_is_empty(queue));

cleanup:
  array_deinit(nb);
  array_dealloc(&nb);

  array_deinit(queue);
  array_dealloc(&queue);

  return done;
}

static void init_bdf_labels(mesh3_s *mesh) {
  /* Allocate and initialize all labels with `NO_LABEL` */
  mesh->bdf_label = malloc(mesh->nbdf*sizeof(size_t));
  for (size_t i = 0; i < mesh->nbdf; ++i)
    mesh->bdf_label[i] = NO_LABEL;

  /* No labels initially */
  mesh->num_bdf_labels = 0;

  /* Label each reflecting component, incrementing the label count as
   * we go */
  while (!label_reflector(mesh))
    ++mesh->num_bdf_labels;
}

static size_t find_bde(mesh3_s const *mesh, bde_s const *bde) {
  bde_s const *found = bsearch(
    bde, mesh->bde, mesh->nbde, sizeof(bde_s), (compar_t)bde_cmp);
  return (size_t)(found ? found - mesh->bde : NO_INDEX);
}

static void get_diff_bde_nbs(mesh3_s const *mesh, size_t le, array_s *nb) {
  bde_s const *bde = &mesh->bde[le];
  bde_s nb_bde;

  assert(bde->diff);

  for (size_t i = 0, l, nvv, *vv; i < 2; ++i) {
    l = bde->le[i];

    nvv = mesh3_nvv(mesh, l);
    vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(mesh, l, vv);

    for (size_t j = 0, le_nb; j < nvv; ++j) {
      nb_bde = make_bde(l, vv[j]);
      if (bde_cmp(bde, &nb_bde) == 0)
        continue;
      le_nb = find_bde(mesh, &nb_bde);
      if (le_nb != (size_t)NO_INDEX
          && mesh->bde[le_nb].diff /* only append diffracting edges! */
          && !array_contains(nb, &le_nb))
        array_append(nb, &le_nb);
    }

    free(vv);
  }
}

static bool bdes_are_colinear(mesh3_s const *mesh, size_t l0, size_t l1) {
  size_t const *le[2] = {mesh->bde[l0].le, mesh->bde[l1].le};

  line3 line;
  mesh3_copy_vert(mesh, le[0][0], line.x);
  mesh3_copy_vert(mesh, le[0][1], line.y);

  dbl const *x1[2];
  x1[0] = mesh3_get_vert_ptr(mesh, le[1][0]);
  x1[1] = mesh3_get_vert_ptr(mesh, le[1][1]);

  dbl const atol = mesh->eps;
  return line3_point_colinear(&line, x1[0], atol)
    && line3_point_colinear(&line, x1[1], atol);
}

static bool label_diffractor(mesh3_s *mesh) {
  bool done = false;

  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  array_s *nb;
  array_alloc(&nb);
  array_init(nb, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  size_t le = 0;
  while (le < mesh->nbde
         && (mesh->bde_label[le] != NO_LABEL
             || !mesh->bde[le].diff))
    ++le;

  if (le == mesh->nbde) {
    done = true;
    goto cleanup;
  }

  array_append(queue, &le);

  do {
    array_pop_front(queue, &le);

    mesh->bde_label[le] = mesh->num_bde_labels;

    get_diff_bde_nbs(mesh, le, nb);

    size_t le_nb;
    while (!array_is_empty(nb)) {
      array_pop_front(nb, &le_nb);

      if (mesh->bde_label[le_nb] != NO_LABEL)
        continue;

      if (bdes_are_colinear(mesh, le, le_nb)
          && !array_contains(queue, &le_nb))
        array_append(queue, &le_nb);
    }
  } while (!array_is_empty(queue));

cleanup:
  array_deinit(nb);
  array_dealloc(&nb);

  array_deinit(queue);
  array_dealloc(&queue);

  return done;
}

static void init_bde_labels(mesh3_s *mesh) {
  /* Allocate and initialize all labels with `NO_LABEL` */
  mesh->bde_label = malloc(mesh->nbde*sizeof(size_t));
  for (size_t i = 0; i < mesh->nbde; ++i)
    mesh->bde_label[i] = NO_LABEL;

  /* No labels initially */
  mesh->num_bde_labels = 0;

  /* Label each diffractor */
  while (!label_diffractor(mesh))
    ++mesh->num_bde_labels;
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
  bdf_s *f = malloc(nf*sizeof(bdf_s));

  // Traverse the cells in the mesh, and populate `f`. These faces are
  // "tagged", meaning that they have a backpointer to the originating
  // cell.
  size_t *C;
  for (size_t lc = 0; lc < mesh->ncells; ++lc) {
    C = mesh->cells[lc];
    f[4*lc] = make_bdf(C[0], C[1], C[2], lc);
    f[4*lc + 1] = make_bdf(C[0], C[1], C[3], lc);
    f[4*lc + 2] = make_bdf(C[0], C[2], C[3], lc);
    f[4*lc + 3] = make_bdf(C[1], C[2], C[3], lc);
  }

  // Sort the tagged faces themselves into a dictionary order.
  qsort(f, nf, sizeof(bdf_s), (compar_t)bdf_cmp);

  /**
   * Set up the boundary vertex, cell, face data structures (stored in
   * `mesh->bdv`, `mesh->bdv`, and `mesh->bdf`, respectively).
   */

  mesh->nbdf = 0;

  // After sorting, we can tell if a face is duplicated by checking
  // whether it's equal to the succeeding face in `f`. If there's a
  // duplicate, we update the relevant boundary information.
  for (size_t l = 0; l < nf - 1; ++l) {
    if (!bdf_cmp(&f[l], &f[l + 1])) {
      ++l;
      continue;
    }
    mesh->bdc[f[l].lc] = true;
    mesh->bdv[f[l].lf[0]] = true;
    mesh->bdv[f[l].lf[1]] = true;
    mesh->bdv[f[l].lf[2]] = true;
    ++mesh->nbdf;
  }

  // Don't forget to check the last face! It could be a boundary face,
  // too, in which case the foregoing loop will miss it
  if (bdf_cmp(&f[nf - 2], &f[nf - 1])) {
    mesh->bdc[f[nf - 1].lc] = true;
    mesh->bdv[f[nf - 1].lf[0]] = true;
    mesh->bdv[f[nf - 1].lf[1]] = true;
    mesh->bdv[f[nf - 1].lf[2]] = true;
    ++mesh->nbdf;
  }

  size_t lf = 0;

  // Traverse the sorted list again and pull out the boundary faces
  // now that we know how many there are
  mesh->bdf = malloc(mesh->nbdf*sizeof(bdf_s));
  for (size_t l = 0; l < nf - 1; ++l) {
    if (!bdf_cmp(&f[l], &f[l + 1])) {
      ++l; // Increment here to skip equal pairs
      continue;
    }
    bdf_init(&mesh->bdf[lf++], f[l].lf, f[l].lc);
  }

  // ... and the last face
  if (bdf_cmp(&f[nf - 2], &f[nf - 1])) {
    bdf_init(&mesh->bdf[lf++], f[nf - 1].lf, f[nf - 1].lc);
  }

  assert(lf == mesh->nbdf); // sanity
  assert(lf > 0);           // check

  // Sort the faces so that we can quickly query whether a face is a
  // boundary face or not.
  qsort(mesh->bdf, mesh->nbdf, sizeof(bdf_s), (compar_t)bdf_cmp);

  /**
   * Set up the boundary edge data structure (stored in `mesh->bde`)
   */

  // Build a sorted array of all of the boundary edges, which are just
  // the edges incident on the boundary faces. This array will have
  // duplicates, so we'll have to prune these next.
  bde_s *bde = malloc(3*mesh->nbdf*sizeof(bde_s));
  for (size_t lf = 0, *l; lf < mesh->nbdf; ++lf) {
    l = mesh->bdf[lf].lf;
    bde[3*lf] = make_bde(l[0], l[1]);
    bde[3*lf + 1] = make_bde(l[1], l[2]);
    bde[3*lf + 2] = make_bde(l[2], l[0]);
  }
  qsort(bde, 3*mesh->nbdf, sizeof(bde_s), (compar_t)bde_cmp);

  // Now, let's count the number of distinct boundary edges.
  mesh->nbde = 1; // count the first edge (we assume nbdf > 0)
  for (size_t l = 0; l < 3*mesh->nbdf - 1; ++l)
    if  (bde_cmp(&bde[l], &bde[l + 1]))
      ++mesh->nbde;

  size_t le = 0;

  // Traverse the array again, copying over distinct boundary
  // edges. Note: there's no need to sort mesh->bde, since it will
  // already be sorted.
  mesh->bde = malloc(mesh->nbde*sizeof(bde_s));
  mesh->bde[le++] = bde[0];
  for (size_t l = 0; l < 3*mesh->nbdf - 1; ++l)
    if (bde_cmp(&bde[l], &bde[l + 1]))
      mesh->bde[le++] = bde[l + 1];

  assert(le == mesh->nbde); // sanity
  assert(le > 0);           // check

  /* Check whether each boundary edge is a diffracting edge */
  for (size_t l = 0; l < mesh->nbde; ++l)
    mesh->bde[l].diff = edge_is_diff(mesh, mesh->bde[l].le);

  // Cleanup
  free(bde);
  free(f);
}

static void compute_geometric_quantities(mesh3_s *mesh) {
  // Compute minimum tetrahedron altitude
  mesh->min_tetra_alt = INFINITY;
  for (size_t l = 0; l < mesh->ncells; ++l) {
    dbl x[4][3];
    for (int i = 0; i < 4; ++i)
      mesh3_copy_vert(mesh, mesh->cells[l][i], x[i]);
    dbl h = min_tetra_altitude(x);
    mesh->min_tetra_alt = fmin(mesh->min_tetra_alt, h);
  }

  /* Compute the minimum edge length */
  mesh->min_edge_length = INFINITY;
  for (size_t i = 0; i < mesh->nedges; ++i) {
    size_t const *le = mesh->edges[i];
    dbl h = dbl3_dist(mesh->verts[le[0]], mesh->verts[le[1]]);
    mesh->min_edge_length = fmin(h, mesh->min_edge_length);
  }

  /* Compute mean edge length */
  mesh->mean_edge_length = 0;
  for (size_t i = 0; i < mesh->nedges; ++i) {
    size_t const *le = mesh->edges[i];
    dbl h = dbl3_dist(mesh->verts[le[0]], mesh->verts[le[1]]);
    mesh->mean_edge_length += h;
  }
  mesh->mean_edge_length /= mesh->nedges;

  /* Approximate the diameter */
  mesh->diam = mesh3_diam_2approx_rand(mesh, 100, NULL);
}

void mesh3_init(mesh3_s *mesh, mesh3_data_s const *data,
                bool compute_bd_info, dbl const *eps) {
  mesh->verts = malloc(data->nverts*sizeof(dbl3));
  memcpy(mesh->verts, data->verts, data->nverts*sizeof(dbl3));
  mesh->nverts = data->nverts;

  mesh->cells = malloc(data->ncells*sizeof(uint4));
  memcpy(mesh->cells, data->cells, data->ncells*sizeof(uint4));
  mesh->ncells = data->ncells;

  init_vc(mesh);

  init_edges(mesh);

  compute_geometric_quantities(mesh);

  mesh->eps = eps ? *eps : EPS;

  mesh->has_bd_info = compute_bd_info;
  if (compute_bd_info) {
    init_bd(mesh);
    init_bdf_labels(mesh);
    init_bde_labels(mesh);
  }
}

void mesh3_deinit(mesh3_s *mesh) {
  free(mesh->verts);
  free(mesh->cells);
  free(mesh->edges);
  free(mesh->vc);
  free(mesh->vc_offsets);

  mesh->verts = NULL;
  mesh->cells = NULL;
  mesh->edges = NULL;
  mesh->vc = NULL;
  mesh->vc_offsets = NULL;

  if (mesh->has_bd_info) {
    free(mesh->bdc);
    free(mesh->bdv);
    free(mesh->bdf);
    free(mesh->bde);
    free(mesh->bdf_label);
    free(mesh->bde_label);

    mesh->bdc = NULL;
    mesh->bdv = NULL;
    mesh->bdf = NULL;
    mesh->bde = NULL;
    mesh->bdf_label = NULL;
    mesh->bde_label = NULL;
  }
}

dbl3 const *mesh3_get_verts_ptr(mesh3_s const *mesh) {
  return mesh->verts;
}

size_t const *mesh3_get_cells_ptr(mesh3_s const *mesh) {
  return mesh->cells[0];
}

dbl const *mesh3_get_vert_ptr(mesh3_s const *mesh, size_t i) {
  return mesh->verts[i];
}

void mesh3_get_vert_ptrs(mesh3_s const *mesh, size_t const *l, int n,
                         dbl const **x) {
  for (int i = 0; i < n; ++i)
    x[i] = mesh->verts[l[i]];
}

void mesh3_copy_vert(mesh3_s const *mesh, size_t i, dbl *v) {
  memcpy(v, mesh->verts[i], 3*sizeof(dbl));
}

tetra3 mesh3_get_tetra(mesh3_s const *mesh, size_t lc) {
  tetra3 tetra;
  for (int i = 0; i < 4; ++i)
    mesh3_copy_vert(mesh, mesh->cells[lc][i], tetra.v[i]);
  return tetra;
}

tri3 mesh3_get_tri(mesh3_s const *mesh, size_t const lf[3]) {
  tri3 tri;
  for (size_t i = 0; i < 3; ++i)
    mesh3_copy_vert(mesh, lf[i], tri.v[i]);
  return tri;
}

void mesh3_get_centroid(mesh3_s const *mesh, size_t lc, dbl c[3]) {
  tetra3 tetra = mesh3_get_tetra(mesh, lc);
  for (int j = 0; j < 3; ++j) {
    c[j] = 0;
    for (int i = 0; i < 4; ++i)
      c[j] += tetra.v[i][j];
    c[j] /= 4;
  }
}

void mesh3_get_edge_centroid(mesh3_s const *mesh, size_t e[2], dbl c[3]) {
  assert(mesh3_is_edge(mesh, e));
  dbl3_cc(mesh->verts[e[0]], mesh->verts[e[1]], 0.5, c);
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
    v = mesh->verts[l];
    min[0] = fmin(min[0], v[0]);
    min[1] = fmin(min[1], v[1]);
    min[2] = fmin(min[2], v[2]);
  }

  max[0] = max[1] = max[2] = -INFINITY;
  for (size_t l = 0; l < mesh->nverts; ++l) {
    v = mesh->verts[l];
    max[0] = fmax(max[0], v[0]);
    max[1] = fmax(max[1], v[1]);
    max[2] = fmax(max[2], v[2]);
  }
}

void mesh3_get_cell_bbox(mesh3_s const *mesh, size_t i, rect3 *bbox) {
  size_t const *cell = mesh->cells[i];

  dbl *min = bbox->min, *max = bbox->max, *v;

  min[0] = min[1] = min[2] = INFINITY;
  for (int i = 0; i < 4; ++i) {
    v = mesh->verts[cell[i]];
    min[0] = fmin(min[0], v[0]);
    min[1] = fmin(min[1], v[1]);
    min[2] = fmin(min[2], v[2]);
  }

  max[0] = max[1] = max[2] = -INFINITY;
  for (int i = 0; i < 4; ++i) {
    v = mesh->verts[cell[i]];
    max[0] = fmax(max[0], v[0]);
    max[1] = fmax(max[1], v[1]);
    max[2] = fmax(max[2], v[2]);
  }
}

bool mesh3_cell_contains_point(mesh3_s const *mesh, size_t lc, dbl const x[3]) {
  tetra3 tetra = mesh3_get_tetra(mesh, lc);
  return tetra3_contains_point(&tetra, x, &mesh->eps);
}

size_t mesh3_find_cell_containing_point(mesh3_s const *mesh, dbl const x[3],
                                        size_t lc) {
  /* First, try the passed guess... */
  if (lc != (size_t)NO_INDEX && mesh3_cell_contains_point(mesh, lc, x))
    return lc;

  /* ... otherwise, proceed as usual */
  for (size_t lc = 0; lc < mesh->ncells; ++lc)
    if (mesh3_cell_contains_point(mesh, lc, x))
      return lc;
  return (size_t)NO_INDEX;
}

bool mesh3_contains_point(mesh3_s const *mesh, dbl3 const x) {
  size_t lc = mesh3_find_cell_containing_point(mesh, x, (size_t)NO_INDEX);
  return lc != (size_t)NO_INDEX;
}

/* This function is inexact. Instead of checking containment to
 * machine precision, we check the distance from `x` to each boundary
 * vertex. If `x` is contained in the mesh and the distance to the
 * nearest boundary vertex is greater than or equal to `r`, then we
 * declare victory. */
bool mesh3_contains_ball(mesh3_s const *mesh, dbl3 const x, dbl r) {
  if (!mesh3_contains_point(mesh, x))
    return false;
  for (size_t l = 0; l < mesh->nverts; ++l)
    if (mesh3_bdv(mesh, l) && dbl3_dist(x, mesh->verts[l]) < r)
      return false;
  return true;
}

int mesh3_nvc(mesh3_s const *mesh, size_t i) {
  assert(i < mesh->nverts);
  return mesh->vc_offsets[i + 1] - mesh->vc_offsets[i];
}

void mesh3_vc(mesh3_s const *mesh, size_t i, size_t *vc) {
  int nvc = mesh3_nvc(mesh, i);
  size_t *vci = &mesh->vc[mesh->vc_offsets[i]];
  memcpy((void *)vc, (void *)vci, sizeof(size_t)*nvc);
}

static void get_opposite_edges(size_t const cv[4], size_t lv, edge_s edge[3]) {
  size_t l[3];
  for (int i = 0, j = 0; i < 4; ++i) {
    if (cv[i] == lv)
      continue;
    l[j++] = cv[i];
  }
  edge[0] = make_edge(l[0], l[1]);
  edge[1] = make_edge(l[1], l[2]);
  edge[2] = make_edge(l[2], l[0]);
}

int mesh3_nve(mesh3_s const *mesh, size_t lv) {
  array_s *edges;
  array_alloc(&edges);
  array_init(edges, sizeof(edge_s), /* capacity */ 8);

  int nvc = mesh3_nvc(mesh, lv);
  size_t *vc = malloc(sizeof(size_t)*nvc);
  mesh3_vc(mesh, lv, vc);

  size_t cv[4];

  edge_s new_edges[3];
  for (int i = 0; i < nvc; ++i) {
    mesh3_cv(mesh, vc[i], cv);
    get_opposite_edges(cv, lv, new_edges);
    for (int j = 0; j < 3; ++j) {
      if (array_contains(edges, &new_edges[j]))
        continue;
      array_append(edges, &new_edges[j]);
    }
  }

  int nve = array_size(edges);

  free(vc);

  array_deinit(edges);
  array_dealloc(&edges);

  return nve;
}

void mesh3_ve(mesh3_s const *mesh, size_t lv, size_t (*ve)[2]) {
  array_s *edges;
  array_alloc(&edges);
  array_init(edges, sizeof(edge_s), /* capacity */ 8);

  int nvc = mesh3_nvc(mesh, lv);
  size_t *vc = malloc(sizeof(size_t)*nvc);
  mesh3_vc(mesh, lv, vc);

  size_t cv[4];

  edge_s new_edges[3];
  for (int i = 0; i < nvc; ++i) {
    mesh3_cv(mesh, vc[i], cv);
    get_opposite_edges(cv, lv, new_edges);
    for (int j = 0; j < 3; ++j) {
      if (array_contains(edges, &new_edges[j]))
        continue;
      array_append(edges, &new_edges[j]);
    }
  }

  int nve = array_size(edges);

  edge_s edge;
  for (int i = 0; i < nve; ++i) {
    array_get(edges, i, &edge);
    ve[i][0] = edge.l[0];
    ve[i][1] = edge.l[1];
  }

  free(vc);

  array_deinit(edges);
  array_dealloc(&edges);
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

/**
 * TODO: the functions below are simple, unoptimized
 * implementations. There are lots of ways to improve these, but we
 * want to wait on that until later.
 */

int mesh3_nvv(mesh3_s const *mesh, size_t i) {
  int nvc = mesh3_nvc(mesh, i);
  size_t *vc = malloc(sizeof(size_t)*nvc);
  mesh3_vc(mesh, i, vc);

  /* trivial upper bound: nvv <= 3*nvc */
  size_t *vv = malloc(3*sizeof(size_t)*nvc);

  int nvv = 0;
  for (int p = 0; p < nvc; ++p) {
    size_t const *cell = mesh->cells[vc[p]];
    for (int q = 0; q < 4; ++q) {
      size_t j = cell[q];
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
    size_t const *cell = mesh->cells[vc[p]];
    for (int q = 0; q < 4; ++q) {
      size_t j = cell[q];
      if (i == j || contains((void *)vv, k, &j, sizeof(size_t))) {
        continue;
      }
      vv[k++] = j;
    }
  }

  free(vc);
}

static int num_shared_verts(size_t const *cell1, size_t const *cell2) {
  // TODO: speed up using SIMD?
  int n = 0;
  for (int p = 0; p < 4; ++p) {
    for (int q = 0; q < 4; ++q) {
      if (cell1[p] == cell2[q]) {
        ++n;
        continue;
      }
    }
  }
  return n;
}

int mesh3_ncc(mesh3_s const *mesh, size_t i) {
  size_t const *cell = mesh->cells[i];

  int nvc[4], max_nvc = -1;
  for (int p = 0; p < 4; ++p) {
    nvc[p] = mesh3_nvc(mesh, cell[p]);
    if (nvc[p] > max_nvc) max_nvc = nvc[p];
  }

  size_t *vc = malloc(sizeof(size_t)*max_nvc);

  size_t *cc = malloc(4*sizeof(size_t));

  int ncc = 0;

  size_t j, k;
  for (int p = 0; p < 4; ++p) {
    j = cell[p];
    mesh3_vc(mesh, j, vc);
    for (int q = 0; q < nvc[p]; ++q) {
      k = vc[q];
      if (i != k
          && num_shared_verts(cell, mesh->cells[k]) == 3
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
  size_t const *cell = mesh->cells[i];

  int nvc[4], max_nvc = -1;
  for (int p = 0; p < 4; ++p) {
    nvc[p] = mesh3_nvc(mesh, cell[p]);
    if (nvc[p] > max_nvc) max_nvc = nvc[p];
  }

  size_t *vc = malloc(sizeof(size_t)*max_nvc);

  int ncc = 0;

  size_t j, k;
  for (int p = 0; p < 4; ++p) {
    j = cell[p];
    mesh3_vc(mesh, j, vc);
    for (int q = 0; q < nvc[p]; ++q) {
      k = vc[q];
      if (i != k
          && num_shared_verts(cell, mesh->cells[k]) == 3
          && !contains(cc, ncc, &k, sizeof(size_t))) {
        cc[ncc++] = k;
      }
    }
  }

  free(vc);
}

/* Fill `lf` with the faces incident on `lc`. These are the four
 * different sets of vertices comprising the faces of the cell indexed
 * by `lc`, but returned in sorted order. */
void mesh3_cf(mesh3_s const *mesh, size_t lc, size_t lf[4][3]) {
  size_t const *cv = mesh->cells[lc];
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0, k = 0; k < 4; ++k) {
      if (i == k)
        continue;
      lf[i][j++] = cv[k];
    }
    SORT3(lf[i][0], lf[i][1], lf[i][2]);
  }
}

void mesh3_cv(mesh3_s const *mesh, size_t i, size_t *cv) {
  memcpy(cv, mesh->cells[i], 4*sizeof(size_t));
}

int mesh3_nec(mesh3_s const *mesh, size_t const le[2]) {
  assert(le[0] != le[1]);

  size_t i = le[0], j = le[1];

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
        break;
      }
    }
  }

  free(vci);
  free(vcj);

  return nec;
}

void mesh3_ec(mesh3_s const *mesh, size_t const le[2], size_t *lc) {
  assert(le[0] != le[1]);

  size_t i = le[0], j = le[1];

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
        lc[nec++] = c;
        break;
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

int mesh3_nee(mesh3_s const *mesh, size_t const le[2]) {
  return mesh3_nec(mesh, le);
}

void mesh3_ee(mesh3_s const *mesh, size_t const le[2], size_t (*ee)[2]) {
  int nec = mesh3_nec(mesh, le);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, le, ec);

  for (int i = 0; i < nec; ++i)
    mesh3_cee(mesh, ec[i], le, ee[i]);
  orient_chains(ee, nec);

  free(ec);
}

size_t mesh3_nev(mesh3_s const *mesh, size_t const le[2]) {
  size_t nec = mesh3_nec(mesh, le);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, le, ec);

  array_s *ev;
  array_alloc(&ev);
  array_init(ev, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  size_t ee[2];

  /* Fill `fev` with the vertices opposite `e` in each face incident
   * on `e`. We'll count the number of elements in `fev` to calculate
   * `nef`. */
  for (size_t i = 0; i < nec; ++i) {
    mesh3_cee(mesh, ec[i], le, ee);
    for (size_t j = 0; j < 2; ++j) {
      if (array_contains(ev, &ee[j]))
        continue;
      array_append(ev, &ee[j]);
    }
  }

  size_t nev = array_size(ev);

  array_deinit(ev);
  array_dealloc(&ev);

  free(ec);

  return nev;
}

void mesh3_ev(mesh3_s const *mesh, size_t const e[2], size_t *v) {
  size_t nec = mesh3_nec(mesh, e);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, e, ec);

  array_s *ev;
  array_alloc(&ev);
  array_init(ev, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  size_t ee[2];

  /* Fill `fev` with the vertices opposite `e` in each face incident
   * on `e`. We'll count the number of elements in `fev` to calculate
   * `nef`. */
  size_t nev = 0;
  for (size_t i = 0; i < nec; ++i) {
    mesh3_cee(mesh, ec[i], e, ee);
    for (size_t j = 0; j < 2; ++j) {
      if (array_contains(ev, &ee[j]))
        continue;
      array_append(ev, &ee[j]);
      v[nev++] = ee[j];
    }
  }

  array_deinit(ev);
  array_dealloc(&ev);

  free(ec);
}

int mesh3_nfc(mesh3_s const *mesh, size_t const f[3]) {
  // We find the cells neighboring one of the vertices of the face and
  // then determine which ones are adjacent to the rest of the
  // vertices. It doesn't matter which vertex of `f` we use to do
  // this.

  int nvc = mesh3_nvc(mesh, f[0]);
  size_t *vc = malloc(nvc*sizeof(size_t)), cv[4];
  mesh3_vc(mesh, f[0], vc);

  int nfc = 0;
  for (int i = 0; i < nvc; ++i) {
    mesh3_cv(mesh, vc[i], cv);
    nfc += face_in_cell(f, cv);
  }
  assert(nfc == 1 || nfc == 2);

  free(vc);

  return nfc;
}

void mesh3_fc(mesh3_s const *mesh, size_t const f[3], size_t *fc) {
  int nvc = mesh3_nvc(mesh, f[0]);
  size_t *vc = malloc(nvc*sizeof(size_t)), cv[4];
  mesh3_vc(mesh, f[0], vc);

  int nfc = 0;
  for (int i = 0; i < nvc; ++i) {
    mesh3_cv(mesh, vc[i], cv);
    if (face_in_cell(f, cv))
      fc[nfc++] = vc[i];
  }

  free(vc);
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

bool mesh3_has_bd_info(mesh3_s const *mesh) {
  return mesh->has_bd_info;
}

bool *mesh3_get_bdc_ptr(mesh3_s *mesh) {
  assert(mesh->has_bd_info);
  return mesh->bdc;
}

bool mesh3_bdc(mesh3_s const *mesh, size_t i) {
  assert(mesh->has_bd_info);
  return mesh->bdc[i];
}

bool *mesh3_get_bdv_ptr(mesh3_s *mesh) {
  assert(mesh->has_bd_info);
  return mesh->bdv;
}

bool mesh3_bdv(mesh3_s const *mesh, size_t i) {
  assert(mesh->has_bd_info);
  return mesh->bdv[i];
}

size_t mesh3_nbde(mesh3_s const *mesh) {
  return mesh->nbde;
}

bool mesh3_bde(mesh3_s const *mesh, size_t const le[2]) {
  assert(mesh->has_bd_info);
  bde_s e = make_bde(le[0], le[1]);
  return bsearch(&e, mesh->bde, mesh->nbde, sizeof(bde_s),
                 (compar_t)bde_cmp);
}

size_t mesh3_nbdf(mesh3_s const *mesh) {
  return mesh->nbdf;
}

void mesh3_get_bdf_inds(mesh3_s const *mesh, size_t l, size_t lf[3]) {
  assert(mesh->has_bd_info);
  assert(l < mesh->nbdf);
  memcpy(lf, mesh->bdf[l].lf, sizeof(size_t[3]));
}

void mesh3_set_bdf(mesh3_s *mesh, size_t lf[3]) {
  // TODO: this is very inefficient!

  mesh->bdv[lf[0]] = mesh->bdv[lf[1]] = mesh->bdv[lf[2]] = true;

  assert(mesh->has_bd_info);
  bdf_s bdf = make_bdf(lf[0], lf[1], lf[2], NO_PARENT);

  int cmp;
  size_t l = 0;
  while ((cmp = bdf_cmp(&bdf, &mesh->bdf[l])) >= 0)
    ++l;

  if (cmp == 0)
    return;

  mesh->bdf = realloc(mesh->bdf, (mesh->nbdf + 1)*sizeof(bdf_s));
  memmove(&mesh->bdf[l + 1], &mesh->bdf[l],
          (mesh->nbdf - l)*sizeof(bdf_s));
  mesh->bdf[l] = bdf;
  ++mesh->nbdf;
}

bool mesh3_is_bdf(mesh3_s const *mesh, size_t const lf[3]) {
  assert(mesh->has_bd_info);
  bdf_s f = make_bdf(lf[0], lf[1], lf[2], NO_PARENT);
  return bsearch(
    &f, mesh->bdf, mesh->nbdf, sizeof(bdf_s), (compar_t)bdf_cmp) != NULL;
}

bool mesh3_is_edge(mesh3_s const *mesh, size_t const l[2]) {
  int nvv = mesh3_nvv(mesh, l[0]);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l[0], vv);

  bool is_edge = false;
  for (int i = 0; i < nvv; ++i) {
    if (vv[i] == l[1]) {
      is_edge = true;
      break;
    }
  }

  free(vv);

  return is_edge;
}

bool mesh3_is_diff_edge(mesh3_s const *mesh, size_t const le[2]) {
  assert(mesh->has_bd_info);
  bde_s q = make_bde(le[0], le[1]);
  bde_s const *e = bsearch(
    &q, mesh->bde, mesh->nbde, sizeof(bde_s), (compar_t)bde_cmp);
  return e && e->diff;
}

bool mesh3_vert_incident_on_diff_edge(mesh3_s const *mesh, size_t l) {
  assert(mesh->has_bd_info);

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

bool mesh3_vert_is_terminal_diff_edge_vert(mesh3_s const *mesh, size_t l) {
  assert(mesh->has_bd_info);

  if (!mesh3_vert_incident_on_diff_edge(mesh, l))
    return false;

  array_s *diff_labels_seen;
  array_alloc(&diff_labels_seen);
  array_init(diff_labels_seen, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  size_t nvv = mesh3_nvv(mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l, vv);

  bool terminal = true;

  for (size_t i = 0; i < nvv; ++i) {
    if (!mesh3_is_diff_edge(mesh, (size_t[2]) {l, vv[i]}))
      continue;
    bde_s bde = make_bde(l, vv[i]);
    size_t j = find_bde(mesh, &bde);
    size_t bde_label = mesh->bde_label[j];
    assert(bde_label != NO_LABEL);
    if (array_contains(diff_labels_seen, &bde_label)) {
      terminal = false;
      break;
    } else {
      array_append(diff_labels_seen, &bde_label);
    }
  }

  free(vv);

  array_deinit(diff_labels_seen);
  array_dealloc(&diff_labels_seen);

  return terminal;
}

bool mesh3_cell_incident_on_diff_edge(mesh3_s const *mesh, size_t lc) {
  for (size_t i = 0; i < 4; ++i)
    if (mesh3_vert_incident_on_diff_edge(mesh, mesh->cells[lc][i]))
      return true;
  return false;
}

static bool local_ray_in_tetra_cone(mesh3_s const *mesh, dbl3 const p, size_t lc, size_t lv) {
  dbl3 xhat;
  mesh3_copy_vert(mesh, lv, xhat);

  size_t l[3];
  assert(mesh3_cvf(mesh, lc, lv, l));

  dbl3 x[3];
  for (size_t i = 0; i < 3; ++i)
    mesh3_copy_vert(mesh, l[i], x[i]);

  dbl33 dX;
  for (size_t i = 0; i < 3; ++i) {
    dbl3 dx;
    dbl3_sub(x[i], xhat, dx);
    dbl33_set_column(dX, i, dx);
  }

  dbl3 alpha;
  dbl33_dbl3_solve(dX, p, alpha);

  dbl const atol = 1e-13;
  return alpha[0] >= -atol && alpha[1] >= -atol && alpha[2] >= -atol;
}

bool mesh3_local_ray_in_vertex_cone(mesh3_s const *mesh, dbl3 const p, size_t lv) {
  size_t nvc = mesh3_nvc(mesh, lv);
  size_t *vc = malloc(mesh->nverts*sizeof(size_t));
  mesh3_vc(mesh, lv, vc);

  bool in_cone = false;
  for (size_t i = 0; i < nvc; ++i)
    if ((in_cone = local_ray_in_tetra_cone(mesh, p, vc[i], lv)))
      break;

  free(vc);

  return in_cone;
}

bool mesh3_ray_prop_from_edge_is_occluded(mesh3_s const *mesh, dbl3 const t,
                                          uint2 const l) {
  bool occluded = true;

  /* gets cells incident on active edge */
  size_t nec = mesh3_nec(mesh, l);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, l, ec);

  /* check whether `dxhat` points into a tetrahedron incident on the
   * base of the active edge */
  for (size_t i = 0, l_op[2]; i < nec; ++i) {
    /* get the opposite edge */
    mesh3_cee(mesh, ec[i], l, l_op);

    dbl3 x[2];
    mesh3_copy_vert(mesh, l[0], x[0]);
    mesh3_copy_vert(mesh, l[1], x[1]);

    dbl3 y[2];
    mesh3_copy_vert(mesh, l_op[0], y[0]);
    mesh3_copy_vert(mesh, l_op[1], y[1]);

    dbl3 te;
    dbl3_sub(x[1], x[0], te);
    dbl3_normalize(te);

    dbl3 yproj[2];
    dbl3_saxpy(dbl3_dot(te, y[0]) - dbl3_dot(te, x[0]), te, x[0], yproj[0]);
    dbl3_saxpy(dbl3_dot(te, y[1]) - dbl3_dot(te, x[0]), te, x[0], yproj[1]);

    /* Compute the x-axis for the arctan2 computation */
    dbl3 q1;
    dbl3_sub(y[0], yproj[0], q1);
    dbl3_normalize(q1);

    /* Compute the y-axis for the arctan2 computation */
    dbl3 q2;
    dbl3_cross(te, q1, q2);

    /* Check that q2 has the correct orientation (to ensure that we
     * measure the angle using arctan2 consistently) */
    dbl3 u;
    dbl3_sub(y[1], yproj[1], u);
    dbl3_normalize(u);
    if (dbl3_dot(u, q2) < 0)
      dbl3_negate(q2);

    /* Compute the maximum angle */
    dbl thetamax = atan2(dbl3_dot(q2, u), dbl3_dot(q1, u));

    /* ... and compute the angle of interest */
    dbl theta = atan2(dbl3_dot(q2, t), dbl3_dot(q1, t));

    /* Check whether it's in the feasible range */
    dbl const atol = 1e-13;
    if (-atol <= theta && theta <= thetamax + atol) {
      occluded = false;
      break;
    }
  }

  free(ec);

  return occluded;
}

bool mesh3_local_ray_is_occluded(mesh3_s const *mesh, size_t lhat, par3_s const *par) {
  dbl3 xb;
  par3_get_xb(par, mesh, xb);

  dbl const *xhat = mesh->verts[lhat];

  dbl3 dxhat; /* xlam -> xhat */
  dbl3_sub(xhat, xb, dxhat);

  /* Check whether the update ray is contained in the cone spanned by
   * the cells incident on `utetra->lhat`. */

  dbl3 dxlam; /* xhat -> xlam */
  dbl3_sub(xb, xhat, dxlam);

  if (!mesh3_local_ray_in_vertex_cone(mesh, dxlam, lhat))
    return true;

  /* Make sure the start of the update ray doesn't exit the domain */

  bool ray_start_is_feasible = false;

  size_t l_active[3];
  size_t num_active_constraints = par3_get_active_inds(par, l_active);
  num_active_constraints = 3 - num_active_constraints; // TODO: fix this...

  /* interior minimizer */
  if (num_active_constraints == 0) {
    /* get cells incident on base of update */
    size_t nfc = mesh3_nfc(mesh, l_active);
    size_t *fc = malloc(nfc*sizeof(size_t));
    mesh3_fc(mesh, l_active, fc);

    /* check whether the `dxhat` points into a tetrahedron incident on
     * the base of the update */
    for (size_t i = 0, lv; i < nfc; ++i) {
      mesh3_cfv(mesh, fc[i], l_active, &lv);

      size_t lf[3];
      mesh3_cvf(mesh, fc[i], lv, lf);

      dbl3 x[3];
      for (size_t j = 0; j < 3; ++j)
        mesh3_copy_vert(mesh, lf[j], x[j]);

      dbl3 dx[3];
      dbl3_sub(mesh3_get_vert_ptr(mesh, lv), x[0], dx[0]);
      dbl3_sub(x[1], x[0], dx[1]);
      dbl3_sub(x[2], x[0], dx[2]);

      /* compute the face normal and orient it so that it points
       * outside the tetrahedron */
      dbl3 n;
      dbl3_cross(dx[1], dx[2], n);
      if (dbl3_dot(n, dx[0]) < 0)
        dbl3_negate(n);

      if (dbl3_dot(n, dxhat) > 0) {
        ray_start_is_feasible = true;
        break;
      }
    }

    free(fc);
  }

  /* edge minimizer */
  else if (num_active_constraints == 1) {
    dbl3 x[2];
    mesh3_copy_vert(mesh, l_active[0], x[0]);
    mesh3_copy_vert(mesh, l_active[1], x[1]);

    dbl3 te;
    dbl3_sub(x[1], x[0], te);
    dbl3_normalize(te);

    dbl3 xproj;
    dbl3_saxpy(dbl3_dot(te, xhat) - dbl3_dot(te, x[0]), te, x[0], xproj);

    dbl3 t;
    dbl3_sub(xhat, xproj, t);
    dbl3_normalize(t);

    ray_start_is_feasible = !mesh3_ray_prop_from_edge_is_occluded(mesh, t, l_active);
  }

  /* corner minimizer */
  else if (num_active_constraints == 2)
    ray_start_is_feasible = mesh3_local_ray_in_vertex_cone(mesh, dxhat, l_active[0]);

  else assert(false);

  return !ray_start_is_feasible;
}

dbl mesh3_get_min_tetra_alt(mesh3_s const *mesh) {
  return mesh->min_tetra_alt;
}

dbl mesh3_get_min_edge_length(mesh3_s const *mesh) {
  return mesh->min_edge_length;
}

dbl mesh3_get_mean_edge_length(mesh3_s const *mesh) {
  return mesh->mean_edge_length;
}

/**
 * Create and return a new mesh2 consisting of the boundary of the
 * tetrahedron mesh stored in `mesh`. The returned mesh must be
 * destroyed by the caller.
 */
mesh2_s *mesh3_get_surface_mesh(mesh3_s const *mesh) {
  assert(mesh->has_bd_info);

  // TODO: Right now we're just splatting the relevant vertices and
  // faces into a new mesh2. It wouldn't be too difficult to compress
  // the vertex array. An even cooler thing to do would be to only
  // create a faces array which references the vertex array in the
  // generating mesh3. But at that point, we're probably looking at
  // smart pointers, so we'll hold off on that for now.

  // Grab all of the referenced vertices and splat them into a new
  // array. There will be duplicate vertices!
  dbl3 *verts = malloc(3*mesh->nbdf*sizeof(dbl3));
  for (size_t lf = 0; lf < mesh->nbdf; ++lf) {
    bdf_s const *face = &mesh->bdf[lf];
    for (int i = 0; i < 3; ++i) {
      size_t lv = face->lf[i];
      for (int j = 0; j < 3; ++j)
        verts[3*lf + i][j] = mesh->verts[lv][j];
    }
  }

  // The faces array is very easy to construct at this point---it's
  // just an array of 3*nbdf size_t's, filled with the values 0, ...,
  // 3*mesh->nbdf - 1.
  uint3 *faces = malloc(mesh->nbdf*sizeof(uint3));
  for (size_t lf = 0; lf < mesh->nbdf; ++lf)
    for (size_t i = 0; i < 3; ++i)
      faces[lf][i] = 3*lf + i;

  /* Precompute the array of face normals. The face normal for the
   * tetrahedron mesh is assumed to point *into* the interior of the
   * domain, so to get an outward-facing normal, we need to negate it
   * here. */
  dbl3 *face_normals = malloc(mesh->nbdf*sizeof(dbl3));
  for (size_t lf = 0; lf < mesh->nbdf; ++lf) {
    mesh3_get_face_normal(mesh, mesh->bdf[lf].lf, face_normals[lf]);
    dbl3_negate(face_normals[lf]);
  }

  mesh2_s *surface_mesh;
  mesh2_alloc(&surface_mesh);
  mesh2_init(surface_mesh,
             verts, 3*mesh->nbdf, /* verts_policy: */ POLICY_XFER,
             faces, mesh->nbdf, /* faces_policy: */ POLICY_XFER,
             face_normals, /* face_normals_policy: */ POLICY_XFER);

#if JMM_DEBUG
  /* Make sure no faces have duplicate vertices */
  for (size_t lf = 0; lf < mesh->nbdf; ++lf)
    assert(faces[lf][0] != faces[lf][1] &&
           faces[lf][1] != faces[lf][2] &&
           faces[lf][2] != faces[lf][0]);
#endif

  /* Don't need to free `verts`, `faces`, or `face_normals` here since
   * we've passed ownership to `surface_mesh`. They will be freed when
   * `surface_mesh` is destructed. */

  return surface_mesh;
}

/* Count the number of diffracting edges that have `l` as as
 * endpoint. */
size_t mesh3_get_num_inc_diff_edges(mesh3_s const *mesh, size_t l) {
  assert(mesh->has_bd_info);

  int nvv = mesh3_nvv(mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l, vv);

  size_t num_inc_diff_edges = 0;
  for (int i = 0; i < nvv; ++i)
    if (mesh3_is_diff_edge(mesh, (size_t[2]) {l, vv[i]}))
      ++num_inc_diff_edges;

  free(vv);

  return num_inc_diff_edges;
}

/* Get the diffracting edges incident on `l`. This assumes that `le`
 * has space enough for all edges, the number of which can be found by
 * calling `mesh3_get_num_inc_diff_edges`. */
void mesh3_get_inc_diff_edges(mesh3_s const *mesh, size_t l, size_t (*le)[2]) {
  assert(mesh->has_bd_info);

  int nvv = mesh3_nvv(mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l, vv);

  int j = 0;
  for (int i = 0; i < nvv; ++i) {
    if (mesh3_is_diff_edge(mesh, (size_t[2]) {l, vv[i]})) {
      le[j][0] = l;
      le[j][1] = vv[i];
      ++j;
    }
  }

  free(vv);
}

size_t mesh3_get_num_inc_bdf(mesh3_s const *mesh, size_t l) {
  assert(mesh->has_bd_info);

  if (!mesh3_bdv(mesh, l))
    return 0;

  size_t nve = mesh3_nve(mesh, l);
  size_t (*ve)[2] = malloc(nve*sizeof(size_t[2]));
  mesh3_ve(mesh, l, ve);

  size_t nbdf = 0;
  for (size_t i = 0; i < nve; ++i)
    nbdf += mesh3_is_bdf(mesh, (size_t[3]) {l, ve[i][0], ve[i][1]});

  free(ve);

  return nbdf;
}

void mesh3_get_inc_bdf(mesh3_s const *mesh, size_t l, size_t (*lf)[3]) {
  assert(mesh->has_bd_info);

  if (!mesh3_bdv(mesh, l))
    return;

  size_t nve = mesh3_nve(mesh, l);
  size_t (*ve)[2] = malloc(nve*sizeof(size_t[2]));
  mesh3_ve(mesh, l, ve);

  size_t next_lf[3] = {[0] = l};

  size_t j = 0;
  for (size_t i = 0; i < nve; ++i) {
    memcpy(&next_lf[1], ve[i], sizeof(size_t[2]));
    if (mesh3_is_bdf(mesh, next_lf))
      memcpy(lf[j++], next_lf, sizeof(size_t[3]));
  }

  free(ve);
}

size_t mesh3_get_num_reflectors(mesh3_s const *mesh) {
  return mesh->num_bdf_labels;
}

size_t mesh3_get_reflector_size(mesh3_s const *mesh, size_t i) {
  size_t count = 0;
  for (size_t l = 0; l < mesh->nbdf; ++l)
    count += mesh->bdf_label[l] == i;
  return count;
}

void mesh3_get_reflector(mesh3_s const *mesh, size_t i, size_t (*lf)[3]) {
  size_t nf = 0;
  for (size_t l = 0; l < mesh->nbdf; ++l)
    if (mesh->bdf_label[l] == i)
      memcpy(lf[nf++], mesh->bdf[l].lf, sizeof(size_t[3]));
}

/* Create a `mesh2_s` corresponding to the `i`th reflector. It is the
 * caller's responsibility to clean up the returned mesh. */
mesh2_s *mesh3_get_refl_mesh(mesh3_s const *mesh, size_t i) {
  size_t nf = mesh3_get_reflector_size(mesh, i);
  uint3 *lf = malloc(nf*sizeof(uint3));
  mesh3_get_reflector(mesh, i, lf);

  dbl3 *face_normals = malloc(nf*sizeof(dbl3));
  for (size_t i = 0; i < nf; ++i)
    mesh3_get_face_normal(mesh, lf[i], face_normals[i]);

  mesh2_s *reflector_mesh;
  mesh2_alloc(&reflector_mesh);
  mesh2_init(
    reflector_mesh,
    mesh->verts, mesh->nverts, POLICY_VIEW,
    lf, nf, POLICY_XFER,
    face_normals, POLICY_XFER);

  /* Don't need to free `lf` or `face_normals` here since we transfer
   * ownership to `reflector_mesh`. */

  return reflector_mesh;
}

size_t mesh3_get_num_diffractors(mesh3_s const *mesh) {
  return mesh->num_bde_labels;
}

size_t mesh3_get_diffractor_size(mesh3_s const *mesh, size_t i) {
  size_t count = 0;
  for (size_t l = 0; l < mesh->nbde; ++l)
    count += mesh->bde_label[l] == i;
  return count;
}

void mesh3_get_diffractor(mesh3_s const *mesh, size_t i, size_t (*le)[2]) {
  size_t k = 0;
  for (size_t l = 0; l < mesh->nbde; ++l)
    if (mesh->bde_label[l] == i)
      memcpy(le[k++], mesh->bde[l].le, sizeof(size_t[2]));
}

mesh1_s *mesh3_get_diff_mesh(mesh3_s const *mesh, size_t diff_index) {
  size_t nedges = mesh3_get_diffractor_size(mesh, diff_index);
  uint2 *le = malloc(nedges*sizeof(uint2));
  mesh3_get_diffractor(mesh, diff_index, le);

  mesh1_s *diff_mesh;
  mesh1_alloc(&diff_mesh);
  mesh1_init(diff_mesh,
             mesh->verts, mesh->nverts, POLICY_XFER,
             le, nedges, POLICY_XFER);

  /* Don't need to free `le` since we transfer ownership to
   * `diff_mesh`. */

  return diff_mesh;
}

size_t mesh3_get_num_diffs_inc_on_refl(mesh3_s const *mesh, size_t refl_index) {
  size_t nf = mesh3_get_reflector_size(mesh, refl_index);
  size_t (*lf)[3] = malloc(nf*sizeof(size_t[3]));
  mesh3_get_reflector(mesh, refl_index, lf);

  array_s *diff_index_arr;
  array_alloc(&diff_index_arr);
  array_init(diff_index_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < nf; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      bde_s bde = make_bde(lf[i][j], lf[i][(j + 1) % 3]);
      size_t l = find_bde(mesh, &bde);
      if (l != (size_t)NO_INDEX &&
          mesh->bde[l].diff &&
          !array_contains(diff_index_arr, &mesh->bde_label[l]))
        array_append(diff_index_arr, &mesh->bde_label[l]);
    }
  }

  size_t num_diffs_inc_on_refl = array_size(diff_index_arr);

  array_deinit(diff_index_arr);
  array_dealloc(&diff_index_arr);

  free(lf);

  return num_diffs_inc_on_refl;
}

void mesh3_get_diffs_inc_on_refl(mesh3_s const *mesh, size_t refl_index, size_t *diff_index) {
  size_t nf = mesh3_get_reflector_size(mesh, refl_index);
  size_t (*lf)[3] = malloc(nf*sizeof(size_t[3]));
  mesh3_get_reflector(mesh, refl_index, lf);

  array_s *diff_index_arr;
  array_alloc(&diff_index_arr);
  array_init(diff_index_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < nf; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      bde_s bde = make_bde(lf[i][j], lf[i][(j + 1) % 3]);
      size_t l = find_bde(mesh, &bde);
      if (l != (size_t)NO_INDEX &&
          mesh->bde[l].diff &&
          !array_contains(diff_index_arr, &mesh->bde_label[l]))
        array_append(diff_index_arr, &mesh->bde_label[l]);
    }
  }

  for (size_t i = 0; i < array_size(diff_index_arr); ++i)
    array_get(diff_index_arr, i, &diff_index[i]);

  array_deinit(diff_index_arr);
  array_dealloc(&diff_index_arr);

  free(lf);
}

size_t mesh3_get_diff_index_for_bde(mesh3_s const *mesh, size_t const le[2]) {
  bde_s bde = make_bde(le[0], le[1]);
  size_t l = find_bde(mesh, &bde);
  return mesh->bde_label[l];
}

void mesh3_get_bde_inds(mesh3_s const *mesh, size_t l, size_t le[2]) {
  memcpy(le, mesh->bde[l].le, sizeof(size_t[2]));
}

dbl mesh3_get_eps(mesh3_s const *mesh) {
  return mesh->eps;
}

void mesh3_get_face_normal(mesh3_s const *mesh, size_t const lf[3], dbl normal[3]) {
  assert(mesh3_is_bdf(mesh, lf));

  dbl const *x0 = mesh3_get_vert_ptr(mesh, lf[0]);
  dbl const *x1 = mesh3_get_vert_ptr(mesh, lf[1]);
  dbl const *x2 = mesh3_get_vert_ptr(mesh, lf[2]);

  dbl dx1[3], dx2[3];
  dbl3_sub(x1, x0, dx1);
  dbl3_sub(x2, x0, dx2);

  dbl3_cross(dx1, dx2, normal);
  dbl3_normalize(normal);

  // Get the cell which `lf` is incident on. There should only be one...
  assert(mesh3_nfc(mesh, lf) == 1);
  size_t lc;
  mesh3_fc(mesh, lf, &lc);

  // Get the cell vertex which isn't one of the face vertices
  size_t l3;
  assert(mesh3_cfv(mesh, lc, lf, &l3));

  dbl const *x3 = mesh3_get_vert_ptr(mesh, l3);
  dbl dx3[3]; dbl3_sub(x3, x0, dx3);

  // If the normal isn't pointing into the interior of the domain, flip it
  if (dbl3_dot(normal, dx3) < 0)
    dbl3_negate(normal);
}

void mesh3_get_R_for_face(mesh3_s const *mesh, size_t const lf[3], dbl33 R) {
  dbl3 n;
  mesh3_get_face_normal(mesh, lf, n);

  R_from_n(n, R);
}

void mesh3_get_R_for_reflector(mesh3_s const *mesh, size_t refl_index, dbl33 R) {
  /* Find a face that lies on the reflector, doesn't matter which */
  size_t i = 0;
  while (i < mesh->nbdf && mesh->bdf_label[i] != refl_index)
    ++i;
  assert(mesh->bdf_label[i] == refl_index);

  /* Get the reflection matrix for that face */
  mesh3_get_R_for_face(mesh, mesh->bdf[i].lf, R);
}

void mesh3_get_R_for_interior_reflector_vertex(mesh3_s const *mesh, size_t l, dbl33 R) {
  /* Verify that this is a boundary vertex which is in the interior of
   * a facet. */
  assert(mesh3_bdv(mesh, l));
  assert(!mesh3_vert_incident_on_diff_edge(mesh, l));

  /* Get the boundary faces incident on `l`. */
  size_t nbdf  = mesh3_get_num_inc_bdf(mesh, l);
  assert(nbdf > 0);
  size_t (*lf)[3] = malloc(nbdf*sizeof(size_t[3]));
  mesh3_get_inc_bdf(mesh, l, lf);

  /* All of the incident boundary faces lie on the same planar
   * reflector, so which one we use to compute the reflection matrix
   * `R` doesn't matter. */
  mesh3_get_R_for_face(mesh, lf[0], R);

  free(lf);
}

void mesh3_get_diff_edge_tangent(mesh3_s const *mesh, size_t const le[2],
                                 dbl t[3]) {
  dbl3_sub(mesh->verts[le[1]], mesh->verts[le[0]], t);
  dbl3_normalize(t);
}

dbl mesh3_get_edge_ext_angle(mesh3_s const *mesh, size_t const le[2]) {
  int nec = mesh3_nec(mesh, le);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, le, ec);

  dbl ext_angle = 2*JMM_PI;
  for (int i = 0; i < nec; ++i)
    ext_angle -= get_dihedral_angle(mesh, ec[i], le);

  free(ec);

  return ext_angle;
}

bool mesh3_edge_contains_point(mesh3_s const *mesh, size_t const le[2], dbl const x[3]) {
  line3 line;
  mesh3_copy_vert(mesh, le[0], line.x);
  mesh3_copy_vert(mesh, le[1], line.x);
  return line3_point_in_interval(&line, x, 1e-14);
}

size_t mesh3_get_num_bdf_inc_on_edge(mesh3_s const *mesh, size_t const le[2]) {
  array_s *lf_arr;
  array_alloc(&lf_arr);
  array_init(lf_arr, sizeof(size_t[3]), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < 2; ++i) {
    size_t nlf = mesh3_get_num_inc_bdf(mesh, le[i]);
    size_t (*lf)[3] = malloc(nlf*sizeof(size_t[3]));
    mesh3_get_inc_bdf(mesh, le[i], lf);

    for (size_t j = 0; j < nlf; ++j)
      SORT3(lf[j][0], lf[j][1], lf[j][2]);

    for (size_t j = 0; j < nlf; ++j)
      if (point_in_face(le[0], lf[j]) &&
          point_in_face(le[1], lf[j]) &&
          !array_contains(lf_arr, &lf[j]))
        array_append(lf_arr, &lf[j]);

    free(lf);
  }

  size_t nbdf = array_size(lf_arr);

  array_deinit(lf_arr);
  array_dealloc(&lf_arr);

  return nbdf;
}

void mesh3_get_bdf_inc_on_edge(mesh3_s const *mesh, size_t const le[2],
                                 size_t (*lf)[3]) {
  array_s *lf_arr;
  array_alloc(&lf_arr);
  array_init(lf_arr, sizeof(size_t[3]), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < 2; ++i) {
    size_t nlf = mesh3_get_num_inc_bdf(mesh, le[i]);
    size_t (*lf)[3] = malloc(nlf*sizeof(size_t[3]));
    mesh3_get_inc_bdf(mesh, le[i], lf);

    for (size_t j = 0; j < nlf; ++j)
      SORT3(lf[j][0], lf[j][1], lf[j][2]);

    for (size_t j = 0; j < nlf; ++j)
      if (point_in_face(le[0], lf[j]) &&
          point_in_face(le[1], lf[j]) &&
          !array_contains(lf_arr, &lf[j]))
        array_append(lf_arr, &lf[j]);

    free(lf);
  }

  for (size_t i = 0; i < array_size(lf_arr); ++i)
    memcpy(lf[i], array_get_ptr(lf_arr, i), sizeof(size_t[3]));

  array_deinit(lf_arr);
  array_dealloc(&lf_arr);
}

void mesh3_get_edge_tangent(mesh3_s const *mesh, size_t const le[2], dbl t[3]) {
  dbl3_sub(mesh->verts[le[1]], mesh->verts[le[0]], t);
  dbl3_normalize(t);
}

void mesh3_get_edge_midpoint(mesh3_s const *mesh, size_t const le[2], dbl p[3]) {
  dbl3_cc(mesh->verts[le[0]], mesh->verts[le[1]], 0.5, p);
}

void mesh3_get_face_centroid(mesh3_s const *mesh, size_t const lf[3], dbl p[3]) {
  dbl3_zero(p);
  for (size_t i = 0; i < 3; ++i)
    dbl3_add_inplace(p, mesh->verts[lf[i]]);
  dbl3_dbl_div_inplace(p, 3);
}

void mesh3_dump_verts(mesh3_s const *mesh, char const *path) {
  FILE *fp = fopen(path, "wb");
  fwrite(mesh->verts, sizeof(dbl[3]), mesh->nverts, fp);
  fclose(fp);
}

void mesh3_dump_cells(mesh3_s const *mesh, char const *path) {
  FILE *fp = fopen(path, "wb");
  fwrite(mesh->cells, sizeof(size_t[4]), mesh->ncells, fp);
  fclose(fp);
}

bool mesh3_has_vertex(mesh3_s const *mesh, dbl3 const x) {
  for (size_t l = 0; l < mesh->nverts; ++l)
    if (dbl3_dist(x, mesh->verts[l]) < 1e-13)
      return true;
  return false;
}

size_t mesh3_get_vert_index(mesh3_s const *mesh, dbl3 const x) {
  for (size_t l = 0; l < mesh->nverts; ++l)
    if (dbl3_dist(x, mesh->verts[l]) < 1e-13)
      return l;
  return (size_t)NO_INDEX;
}

dbl mesh3_linterp(mesh3_s const *mesh, dbl const *values, dbl3 const x) {
  size_t lc = mesh3_find_cell_containing_point(mesh, x, (size_t)NO_INDEX);
  tetra3 tetra = mesh3_get_tetra(mesh, lc);
  dbl4 b; tetra3_get_bary_coords(&tetra, x, b);
  size_t lv[4]; mesh3_cv(mesh, lc, lv);
  dbl value = 0;
  for (size_t i = 0; i < 4; ++i)
    value += b[i]*values[lv[i]];
  return value;
}

/* Approximate the diameter of the mesh to within a factor of two. For
 * the vertex with index `l`, find the vertex with the maximum
 * distance to `l`. The true diameter of the mesh is within a factor
 * of two of one half this quantity. */
dbl mesh3_diam_2approx(mesh3_s const *mesh, size_t l) {
  dbl const *x = mesh3_get_vert_ptr(mesh, l);
  dbl rmax = 0;
  for (size_t m = 0; m < mesh3_nverts(mesh); ++m)
    rmax = fmax(rmax, dbl3_dist(x, mesh3_get_vert_ptr(mesh, m)));
  return 2*rmax;
}

/* Randomized algorithm for approximating the diameter of
 * `mesh`. Applies the standard 2-approximation `trials` times and
 * returns the minimum over the trials. If `seed` is passed, it is
 * used to seed the random number generator first. */
dbl mesh3_diam_2approx_rand(mesh3_s const *mesh, size_t trials, size_t const *seed) {
  if (seed)
    srandom(*seed);
  dbl diam_min = INFINITY;
  for (size_t _ = 0; _ < trials; ++_) {
    size_t l = random() % mesh->nverts;
    diam_min = fmin(diam_min, mesh3_diam_2approx(mesh, l));
  }
  return diam_min;
}

dbl mesh3_get_diam(mesh3_s const *mesh) {
  return mesh->diam;
}

dbl mesh3_get_edge_tol(mesh3_s const *mesh, uint2 const le) {
  assert(mesh3_is_edge(mesh, le));
  dbl h = dbl3_dist(mesh->verts[le[0]], mesh->verts[le[1]]);
  return pow(h/mesh->diam, 2);
}

dbl mesh3_get_face_tol(mesh3_s const *mesh, uint3 const lf) {
  dbl h[3] = {
    dbl3_dist(mesh->verts[lf[0]], mesh->verts[lf[1]]),
    dbl3_dist(mesh->verts[lf[1]], mesh->verts[lf[2]]),
    dbl3_dist(mesh->verts[lf[2]], mesh->verts[lf[0]])
  };
  dbl hmin = fmin(h[0], fmin(h[1], h[2]));
  return pow(hmin/mesh->diam, 2);
}
