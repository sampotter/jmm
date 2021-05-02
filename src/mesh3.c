#include "mesh3.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "array.h"
#include "edge.h"
#include "index.h"
#include "macros.h"
#include "mat.h"
#include "mesh2.h"
#include "util.h"

bool face_in_cell(size_t const f[3], size_t const c[4]) {
  return point_in_cell(f[0], c) && point_in_cell(f[1], c) &&
    point_in_cell(f[2], c);
}

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

  bool has_bd_info;
  bool *bdc;
  bool *bdv;
  size_t nbdf;
  tagged_face_s *bdf;
  size_t nbde;
  diff_edge_s *bde;

  size_t nlabels;
  int *bdf_label; // Labels for distinct reflecting surfaces

  // Geometric quantities
  dbl min_tetra_alt; // The minimum of all tetrahedron altitudes in the mesh.
  dbl min_edge_length;
};

tri3 mesh3_tetra_get_face(mesh3_tetra_s const *tetra, int f[3]) {
  size_t *cv = &tetra->mesh->cells[tetra->l].data[0];
  dvec3 *verts = tetra->mesh->verts;
  tri3 tri;
  memcpy(tri.v[0], &verts[cv[f[0]]].data[0], sizeof(dbl[3]));
  memcpy(tri.v[1], &verts[cv[f[1]]].data[0], sizeof(dbl[3]));
  memcpy(tri.v[2], &verts[cv[f[2]]].data[0], sizeof(dbl[3]));
  return tri;
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
      size_t k = mesh->cells[i].data[j];
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

static size_t find_bdf(mesh3_s const *mesh, tagged_face_s const *bdf) {
  tagged_face_s const *found = bsearch(
    bdf, mesh->bdf, mesh->nbdf, sizeof(tagged_face_s),
    (compar_t)tagged_face_cmp);
  return found ? found - mesh->bdf : (size_t)NO_INDEX;
}

/* Find all of the boundary faces that are adjacent to `mesh->bdf[l]`
 * and append them to `nb`. This doesn't assume that `nb` is empty. */
static void get_bdf_nbs(mesh3_s const *mesh, size_t l, array_s *nb) {
  tagged_face_s const *bdf = &mesh->bdf[l];
  tagged_face_s nb_bdf;

  for (size_t i = 0, lv, nve, (*ve)[2]; i < 3; ++i) {
    lv = bdf->lf[i];

    nve = mesh3_nve(mesh, lv);
    ve = malloc(nve*sizeof(size_t[2]));
    mesh3_ve(mesh, lv, ve);

    for (size_t j = 0, l_nb; j < nve; ++j) {
      nb_bdf = make_tagged_face(lv, ve[j][0], ve[j][1], (size_t)NO_INDEX);
      if (tagged_face_cmp(bdf, &nb_bdf) == 0)
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
  dbl const atol = 1e-10;
  return tri3_coplanar(&tri0, &tri1, &atol);
}

static bool label_reflector(mesh3_s *mesh) {
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
  if (lf == mesh->nbdf)
    return false;

  /* Push `lf` onto `queue` to start the BFS. */
  array_append(queue, &lf);

  /* Traverse the faces of the current reflector and label them. */
  do {
    array_pop_front(queue, &lf);

    /* Set the label of the face to the current label. */
    mesh->bdf_label[lf] = mesh->nlabels;

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

  array_deinit(nb);
  array_dealloc(&nb);

  array_deinit(queue);
  array_dealloc(&queue);

  /* We successfully labeled a component. */
  return true;
}

static void init_bdf_labels(mesh3_s *mesh) {
  /* Allocate and initialize all labels with `NO_LABEL` */
  mesh->bdf_label = malloc(mesh->nbdf*sizeof(int));
  for (size_t i = 0; i < mesh->nbdf; ++i)
    mesh->bdf_label[i] = NO_LABEL;

  /* No labels initially */
  mesh->nlabels = 0;

  /* Label each reflecting component, incrementing the label count as
   * we go */
  while (label_reflector(mesh))
    ++mesh->nlabels;
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

  /**
   * Set up the boundary vertex, cell, face data structures (stored in
   * `mesh->bdv`, `mesh->bdv`, and `mesh->bdf`, respectively).
   */

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

  // Don't forget to check the last face! It could be a boundary face,
  // too, in which case the foregoing loop will miss it
  if (tagged_face_cmp(&f[nf - 2], &f[nf - 1])) {
    mesh->bdc[f[nf - 1].lc] = true;
    mesh->bdv[f[nf - 1].lf[0]] = true;
    mesh->bdv[f[nf - 1].lf[1]] = true;
    mesh->bdv[f[nf - 1].lf[2]] = true;
    ++mesh->nbdf;
  }

  size_t lf = 0;

  // Traverse the sorted list again and pull out the boundary faces
  // now that we know how many there are
  mesh->bdf = malloc(mesh->nbdf*sizeof(tagged_face_s));
  for (size_t l = 0; l < nf - 1; ++l) {
    if (!tagged_face_cmp(&f[l], &f[l + 1])) {
      ++l; // Increment here to skip equal pairs
      continue;
    }
    tagged_face_init(&mesh->bdf[lf++], f[l].lf, f[l].lc);
  }

  // ... and the last face
  if (tagged_face_cmp(&f[nf - 2], &f[nf - 1])) {
    tagged_face_init(&mesh->bdf[lf++], f[nf - 1].lf, f[nf - 1].lc);
  }

  assert(lf == mesh->nbdf); // sanity
  assert(lf > 0);           // check

  // Sort the faces so that we can quickly query whether a face is a
  // boundary face or not.
  qsort(mesh->bdf, mesh->nbdf, sizeof(tagged_face_s), (compar_t)tagged_face_cmp);

  /**
   * Set up the boundary edge data structure (stored in `mesh->bde`)
   */

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
  mesh->nbde = 1; // count the first edge (we assume nbdf > 0)
  for (size_t l = 0; l < 3*mesh->nbdf - 1; ++l)
    if  (diff_edge_cmp(&bde[l], &bde[l + 1]))
      ++mesh->nbde;

  size_t le = 0;

  // Traverse the array again, copying over distinct boundary
  // edges. Note: there's no need to sort mesh->bde, since it will
  // already be sorted.
  mesh->bde = malloc(mesh->nbde*sizeof(diff_edge_s));
  mesh->bde[le++] = bde[0];
  for (size_t l = 0; l < 3*mesh->nbdf - 1; ++l)
    if (diff_edge_cmp(&bde[l], &bde[l + 1]))
      mesh->bde[le++] = bde[l + 1];

  assert(le == mesh->nbde); // sanity
  assert(le > 0);           // check

  /* Check whether each boundary edge is a diffracting edge */
  for (size_t l = 0; l < mesh->nbde; ++l)
    mesh->bde[l].diff = edge_is_diff(mesh, &mesh->bde[l]);

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
      mesh3_copy_vert(mesh, mesh->cells->data[i], x[i]);
    dbl h = min_tetra_altitude(x);
    mesh->min_tetra_alt = fmin(mesh->min_tetra_alt, h);
  }

  mesh->min_edge_length = INFINITY;
  for (size_t l = 0, nvv, *vv; l < mesh->nverts; ++l) {
    nvv = mesh3_nvv(mesh, l);
    vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(mesh, l, vv);
    for (size_t i = 0; i < nvv; ++i) {
      dbl h = dbl3_dist(&mesh->verts[l].data[0], &mesh->verts[vv[i]].data[0]);
      mesh->min_edge_length = fmin(mesh->min_edge_length, h);
    }
    free(vv);
  }
}

void mesh3_init(mesh3_s *mesh,
                dbl const *verts, size_t nverts,
                size_t const *cells, size_t ncells,
                bool compute_bd_info) {
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

  mesh->has_bd_info = compute_bd_info;
  if (compute_bd_info) {
    init_bd(mesh);
    init_bdf_labels(mesh);
  }

  compute_geometric_quantities(mesh);
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

  if (mesh->has_bd_info) {
    free(mesh->bdc);
    free(mesh->bdv);
    free(mesh->bdf);
    free(mesh->bde);
    free(mesh->bdf_label);

    mesh->bdc = NULL;
    mesh->bdv = NULL;
    mesh->bdf = NULL;
    mesh->bde = NULL;
    mesh->bdf_label = NULL;
  }
}

dbl const *mesh3_get_verts_ptr(mesh3_s const *mesh) {
  return &mesh->verts[0].data[0];
}

size_t const *mesh3_get_cells_ptr(mesh3_s const *mesh) {
  return &mesh->cells[0].data[0];
}

dvec3 mesh3_get_vert(mesh3_s const *mesh, size_t i) {
  return mesh->verts[i];
}

dbl const *mesh3_get_vert_ptr(mesh3_s const *mesh, size_t i) {
  return &mesh->verts[i].data[0];
}

void mesh3_get_vert_ptrs(mesh3_s const *mesh, size_t const *l, int n,
                         dbl const **x) {
  for (int i = 0; i < n; ++i)
    x[i] = &mesh->verts[l[i]].data[0];
}

void mesh3_copy_vert(mesh3_s const *mesh, size_t i, dbl *v) {
  memcpy(v, &mesh->verts[i].data[0], 3*sizeof(dbl));
}

tetra3 mesh3_get_tetra(mesh3_s const *mesh, size_t lc) {
  tetra3 tetra;
  for (int i = 0; i < 4; ++i)
    mesh3_copy_vert(mesh, mesh->cells[lc].data[i], tetra.v[i]);
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
  dbl3_cc(mesh->verts[e[0]].data, mesh->verts[e[1]].data, 0.5, c);
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

bool mesh3_cell_contains_point(mesh3_s const *mesh, size_t lc, dbl const x[3]) {
  tetra3 tetra = mesh3_get_tetra(mesh, lc);
  return tetra3_contains_point(&tetra, x);
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

/* Fill `lf` with the faces incident on `lc`. These are the four
 * different sets of vertices comprising the faces of the cell indexed
 * by `lc`, but returned in sorted order. */
void mesh3_cf(mesh3_s const *mesh, size_t lc, size_t lf[4][3]) {
  size_t const *cv = mesh->cells[lc].data;
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

size_t mesh3_nev(mesh3_s const *mesh, size_t const e[2]) {
  size_t nec = mesh3_nec(mesh, e[0], e[1]);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, e[0], e[1], ec);

  array_s *ev;
  array_alloc(&ev);
  array_init(ev, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  size_t ee[2];

  /* Fill `fev` with the vertices opposite `e` in each face incident
   * on `e`. We'll count the number of elements in `fev` to calculate
   * `nef`. */
  for (size_t i = 0; i < nec; ++i) {
    mesh3_cee(mesh, ec[i], e, ee);
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
  size_t nec = mesh3_nec(mesh, e[0], e[1]);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(mesh, e[0], e[1], ec);

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
  diff_edge_s e = make_diff_edge(le[0], le[1]);
  return bsearch(&e, mesh->bde, mesh->nbde, sizeof(diff_edge_s),
                 (compar_t)diff_edge_cmp);
}

size_t mesh3_nbdf(mesh3_s const *mesh) {
  return mesh->nbdf;
}

bool mesh3_bdf(mesh3_s const *mesh, size_t const lf[3]) {
  assert(mesh->has_bd_info);
  tagged_face_s f = make_tagged_face(lf[0], lf[1], lf[2], NO_PARENT);
  return bsearch(&f, mesh->bdf, mesh->nbdf, sizeof(tagged_face_s),
                 (compar_t)tagged_face_cmp);
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
  diff_edge_s q = make_diff_edge(le[0], le[1]);
  diff_edge_s const *e = bsearch(
    &q, mesh->bde, mesh->nbde, sizeof(diff_edge_s), (compar_t)diff_edge_cmp);
  return e && e->diff;
}

bool mesh3_is_nondiff_boundary_edge(mesh3_s const *mesh, size_t const le[2]) {
  assert(mesh->has_bd_info);
  diff_edge_s q = make_diff_edge(le[0], le[1]);
  diff_edge_s const *e = bsearch(
    &q, mesh->bde, mesh->nbde, sizeof(diff_edge_s), (compar_t)diff_edge_cmp);
  return e && !e->diff;
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

dbl mesh3_get_min_tetra_alt(mesh3_s const *mesh) {
  return mesh->min_tetra_alt;
}

dbl mesh3_get_min_edge_length(mesh3_s const *mesh) {
  return mesh->min_edge_length;
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
  dbl *verts = malloc(3*mesh->nbdf*sizeof(dbl[3]));
  for (size_t lf = 0; lf < mesh->nbdf; ++lf)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        verts[3*(3*lf + i) + j] = mesh->verts[mesh->bdf[lf].lf[i]].data[j];

  // The faces array is very easy to construct at this point---it's
  // just an array of 3*nbdf size_t's, filled with the values 0, ...,
  // 3*mesh->nbdf - 1.
  size_t *faces = malloc(mesh->nbdf*sizeof(size_t[3]));
  for (size_t lf = 0; lf < 3*mesh->nbdf; ++lf)
    faces[lf] = lf;

  mesh2_s *surface_mesh;
  mesh2_alloc(&surface_mesh);
  mesh2_init(surface_mesh, verts, 3*mesh->nbdf, faces, mesh->nbdf);

  free(faces);
  free(verts);

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
    nbdf += mesh3_bdf(mesh, (size_t[3]) {l, ve[i][0], ve[i][1]});

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
    if (mesh3_bdf(mesh, next_lf))
      memcpy(lf[j++], next_lf, sizeof(size_t[3]));
  }

  free(ve);
}

size_t mesh3_get_num_reflectors(mesh3_s const *mesh) {
  return mesh->nlabels;
}

size_t mesh3_get_reflector_size(mesh3_s const *mesh, int r) {
  size_t count = 0;
  for (size_t l = 0; l < mesh->nbdf; ++l)
    count += mesh->bdf_label[l] == r;
  return count;
}

void mesh3_get_reflector(mesh3_s const *mesh, int r, size_t (*lf)[3]) {
  size_t nf = 0;
  for (size_t l = 0; l < mesh->nbdf; ++l)
    if (mesh->bdf_label[l] == r)
      memcpy(lf[nf++], mesh->bdf[l].lf, sizeof(size_t[3]));
}
