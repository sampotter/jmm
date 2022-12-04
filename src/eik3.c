#include <jmm/eik3.h>

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include <jmm/alist.h>
#include <jmm/array.h>
#include <jmm/bb.h>
#include <jmm/edge.h>
#include <jmm/eik3_transport.h>
#include <jmm/heap.h>
#include <jmm/mat.h>
#include <jmm/mesh1.h>
#include <jmm/mesh2.h>
#include <jmm/mesh3.h>
#include <jmm/utetra.h>
#include <jmm/util.h>
#include <jmm/utri.h>
#include <jmm/utri_cache.h>
#include <jmm/vec.h>

#include "log.h"
#include "macros.h"

static bool l_OK(size_t l) {
  return l != (size_t)NO_INDEX;
}

static bool is_edge(uint3 const l) {
  return l_OK(l[0]) && l_OK(l[1]) && !l_OK(l[2]);
}

static bool is_face(uint3 const l) {
  return l_OK(l[0]) && l_OK(l[1]) && l_OK(l[2]);
}

/* A structure managing a jet marching method solving the eikonal
 * equation in 3D on an unstructured tetrahedron mesh.
 *
 * NOTE: this is just for s = 1 at the moment. Will extend this to
 * handle s != later. */
struct eik3 {
  mesh3_s const *mesh;
  sfunc_s const *sfunc;

  jet31t *jet;
  state_e *state;
  int *pos;
  par3_s *par;
  heap_s *heap;

  /* In some cases, we'll skip old updates that might be useful at a
   * later stage. We keep track of them here. */
  array_s *old_utetra;
  utri_cache_s *old_bd_utri; // old two-point boundary `utri`
  utri_cache_s *old_diff_utri; // old two-point updates from diff. edges

  alist_s *T_diff;

  array_s *trial_inds, *bc_inds;

  /* Useful statistics for debugging */
  size_t num_accepted; /* number of nodes fixed by `eik3_step` */

  /* An array containing the order in which the individual nodes were
   * accepted. That is, `accepted[i] == l` means that `eik3_step()`
   * returned `l` when it was called for the `i`th time. */
  size_t *accepted;

  bool is_initialized;
};

static char const *state_str[] = {
  [FAR] = "FAR",
  [TRIAL] = "TRIAL",
  [VALID] = "VALID"
};

static void print_node(eik3_s const *eik, size_t l) {
  printf("l = %lu: state = %s, T = %.3g, DT = {%.3g, %.3g, %.3g}\n",
         l, state_str[eik->state[l]], eik->jet[l].f,
         eik->jet[l].Df[0], eik->jet[l].Df[1], eik->jet[l].Df[2]);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
static void print_node_and_nbs(eik3_s const *eik, size_t l) {
  print_node(eik, l);

  size_t nvv = mesh3_nvv(eik->mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(eik->mesh, l, vv);

  printf("neighbors:\n");
  for (size_t i = 0; i < nvv; ++i) {
    printf("  ");
    print_node(eik, vv[i]);
  }

  free(vv);
}
#pragma GCC diagnostic pop

void eik3_alloc(eik3_s **eik) {
  *eik = malloc(sizeof(eik3_s));

  (*eik)->is_initialized = false;
}

void eik3_dealloc(eik3_s **eik) {
  free(*eik);
  *eik = NULL;
}

static dbl value(void *ptr, int l) {
  eik3_s *eik = (eik3_s *)ptr;
  assert(l >= 0);
  assert(l < (int)mesh3_nverts(eik->mesh));
  dbl T = eik->jet[l].f;
  return T;
}

static void setpos(void *ptr, int l, int pos) {
  eik3_s *eik = (eik3_s *)ptr;
  eik->pos[l] = pos;
}

void eik3_init(eik3_s *eik, mesh3_s const *mesh, sfunc_s const *sfunc) {
  eik->mesh = mesh;

  eik->sfunc = sfunc;

  size_t nverts = mesh3_nverts(mesh);

  eik->jet = malloc(nverts*sizeof(jet31t));
  for (size_t l = 0; l < nverts; ++l)
    eik->jet[l] = jet31t_make_empty();

  eik->state = malloc(nverts*sizeof(state_e));
  for (size_t l = 0; l < nverts; ++l) {
    eik->state[l] = FAR;
  }

  eik->pos = malloc(nverts*sizeof(int));
  for (size_t l = 0; l < nverts; ++l) {
    eik->pos[l] = NO_INDEX;
  }

  eik->par = malloc(nverts*sizeof(par3_s));
  for (size_t l = 0; l < nverts; ++l) {
    par3_init_empty(&eik->par[l]);
  }

  /**
   * When we compute the initial heap capacity, we want to estimate
   * the number of nodes that could comprise the expanding numerical
   * front at any one time. We can't know this ahead of time, so we
   * set it to a constant multiple times (# nodes)^(1/d). In this
   * case, d=3. Even if this is an underestimate, well still reduce
   * the number of times the heap needs to be expanded at solve time.
   */
  int capacity = (int)3*cbrt(nverts);

  heap_alloc(&eik->heap);
  heap_init(eik->heap, capacity, value, setpos, (void *)eik);

  eik->num_accepted = 0;

  eik->accepted = malloc(nverts*sizeof(size_t));
  for (size_t i = 0; i < nverts; ++i)
    eik->accepted[i] = (size_t)NO_INDEX;

  array_alloc(&eik->old_utetra);
  array_init(eik->old_utetra, sizeof(utetra_s *), 16);

  utri_cache_alloc(&eik->old_bd_utri);
  utri_cache_init(eik->old_bd_utri);

  utri_cache_alloc(&eik->old_diff_utri);
  utri_cache_init(eik->old_diff_utri);

  array_alloc(&eik->bc_inds);
  array_init(eik->bc_inds, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  array_alloc(&eik->trial_inds);
  array_init(eik->trial_inds, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  alist_alloc(&eik->T_diff);
  alist_init(eik->T_diff, sizeof(size_t[2]), sizeof(bb31), ARRAY_DEFAULT_CAPACITY);

  eik->is_initialized = true;
}

void eik3_deinit(eik3_s *eik) {
  free(eik->jet);
  eik->jet = NULL;

  free(eik->state);
  eik->state = NULL;

  free(eik->pos);
  eik->pos = NULL;

  free(eik->par);
  eik->par = NULL;

  free(eik->accepted);
  eik->accepted = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);

  /* release the "old update" data structures */
  {
    utetra_s *utetra;
    for (size_t i = 0; i < array_size(eik->old_utetra); ++i) {
      array_get(eik->old_utetra, i, &utetra);
      utetra_dealloc(&utetra);
    }
    array_deinit(eik->old_utetra);
    array_dealloc(&eik->old_utetra);
  }

  utri_cache_deinit(eik->old_bd_utri);
  utri_cache_dealloc(&eik->old_bd_utri);

  utri_cache_deinit(eik->old_diff_utri);
  utri_cache_dealloc(&eik->old_diff_utri);

  array_deinit(eik->bc_inds);
  array_dealloc(&eik->bc_inds);

  array_deinit(eik->trial_inds);
  array_dealloc(&eik->trial_inds);

  alist_deinit(eik->T_diff);
  alist_dealloc(&eik->T_diff);

  eik->is_initialized = false;
}

bool eik3_is_initialized(eik3_s const *eik) {
  return eik->is_initialized;
}

void eik3_dump_jet(eik3_s const *eik, char const *path) {
  FILE *fp = fopen(path, "wb");
  fwrite(eik->jet, sizeof(eik->jet[0]), mesh3_nverts(eik->mesh), fp);
  fclose(fp);
}

void eik3_dump_state(eik3_s const *eik, char const *path) {
  FILE *fp = fopen(path, "wb");
  fwrite(eik->state, sizeof(eik->state[0]), mesh3_nverts(eik->mesh), fp);
  fclose(fp);
}

void eik3_dump_par_l(eik3_s const *eik, char const *path) {
  FILE *fp = fopen(path, "wb");
  for (size_t i = 0; i < mesh3_nverts(eik->mesh); ++i)
    fwrite(eik->par[i].l, sizeof(eik->par[i].l[0]), 3, fp);
  fclose(fp);
}

void eik3_dump_par_b(eik3_s const *eik, char const *path) {
  FILE *fp = fopen(path, "wb");
  for (size_t i = 0; i < mesh3_nverts(eik->mesh); ++i)
    fwrite(eik->par[i].b, sizeof(eik->par[i].b[0]), 3, fp);
  fclose(fp);
}

void eik3_dump_accepted(eik3_s const *eik, char const *path) {
  FILE *fp = fopen(path, "wb");
  fwrite(eik->accepted, sizeof(eik->accepted[0]), mesh3_nverts(eik->mesh), fp);
  fclose(fp);
}

size_t eik3_peek(eik3_s const *eik) {
  return heap_front(eik->heap);
}

static void adjust(eik3_s *eik, size_t l) {
  assert(eik->state[l] == TRIAL);
  assert(l < mesh3_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l]);
}

/** Functions for `do_utri`: */

static void commit_utri(eik3_s *eik, size_t lhat, utri_s const *utri) {
  /* TODO: see comment about caustics in `commit_utetra` */
  assert(utri_get_value(utri) < eik->jet[lhat].f);

  utri_get_jet31t(utri, &eik->jet[lhat]);

  eik3_set_par(eik, lhat, utri_get_par(utri));
}

/* Check whether the edge indexed by `l` is:
 *   1) a boundary edge (a mesh edge that's  immersed in the boundary)
 *   2) on the `VALID` front (in this case, since we know that `l[0]`
 *   is newly `VALID`, this means that `l[1]` has at least one
 *   boundary vertex neighbor that's `TRIAL`) */
static bool is_valid_front_bde(eik3_s const *eik, size_t const l[2]) {
  mesh3_s const *mesh = eik->mesh;

  if (!mesh3_bde(mesh, l))
    return false;

  size_t nvv = mesh3_nvv(mesh, l[1]);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l[1], vv);

  bool has_trial_nb = false;
  for (size_t i = 0; i < nvv; ++i) {
    if (vv[i] != l[0] &&
        mesh3_bdv(mesh, vv[i]) &&
        eik->state[vv[i]] == TRIAL) {
      has_trial_nb = true;
      break;
    }
  }

  free(vv);

  return has_trial_nb;
}

static bool is_diff_edge(eik3_s const *eik, size_t const l[2]) {
  return mesh3_is_diff_edge(eik->mesh, l);
}

static void
get_valid_inc_edges(eik3_s const *eik, size_t l0, array_s *l1,
                    bool (*pred)(eik3_s const *, size_t const[2])) {
mesh3_s const *mesh = eik3_get_mesh(eik);

  int nvv = mesh3_nvv(mesh, l0);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l0, vv);

  size_t le[2] = {[0] = l0};
  for (int i = 0; i < nvv; ++i) {
    le[1] = vv[i];

    if (!pred(eik, le))
      continue;

    if (!eik3_is_valid(eik, le[1]))
      continue;

    array_append(l1, &le[1]);
  }

  free(vv);
}

static void
do_utri(eik3_s *eik, size_t l, size_t l0, size_t l1, utri_cache_s *utri_cache,
        par3_s *par) {
  if (utri_cache_contains_inds(utri_cache, l, (uint2) {l0, l1}))
    return;

  if (par != NULL)
    par3_init_empty(par);

  utri_s *utri;
  utri_alloc(&utri);
  utri_init(utri, eik, l, (uint2) {l0, l1});

  if (utri_is_backwards(utri, eik))
    goto cleanup;

  if (utri_is_degenerate(utri))
    goto cleanup;

  utri_solve(utri);

  if (par != NULL && utri_get_value(utri) < eik->jet[l].f)
    *par = utri_get_par(utri);

  if (utri_get_value(utri) >= eik->jet[l].f)
    goto cleanup;

  if (utri_ray_is_occluded(utri, eik))
    goto cleanup;

  if (utri_has_interior_point_solution(utri) ||
      utri_active_vert_is_terminal_diff_vert(utri, eik)) {
    commit_utri(eik, l, utri);
    adjust(eik, l);
    goto cleanup;
  }

  /* see if we can find an old triangle update which matches the
   * current one, and commit the update if we can
   *
   * (if successful, delete the old one!) */
  utri_s *utri_other = utri_cache_pop(utri_cache, utri);
  if (utri_other) {
    commit_utri(eik, l, utri);
    adjust(eik, l);
    utri_dealloc(&utri_other);
    goto cleanup;
  }

  /* if we failed, we cache this update for later (... if we
   * haven't already) */
  if (utri_cache_try_add_unique(utri_cache, utri))
    return; /* don't dealloc in this case! */

cleanup:
  utri_dealloc(&utri);
}

static void do_utris_if(eik3_s *eik, size_t l, size_t l0, utri_cache_s *utri_cache,
                        bool (*pred)(eik3_s const *, size_t const[2])) {
  assert(l != l0);

  /* find the diffracting edges incident on l0 with VALID indices */
  array_s *l1_arr;
  array_alloc(&l1_arr);
  array_init(l1_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);
  get_valid_inc_edges(eik, l0, l1_arr, pred);

  for (size_t i = 0, l1; i < array_size(l1_arr); ++i) {
    array_get(l1_arr, i, &l1);
    if (l == l1)
      continue;
    do_utri(eik, l, l0, l1, utri_cache, /* par: */ NULL);
  }

  /* release array of VALID diff edges */
  array_deinit(l1_arr);
  array_dealloc(&l1_arr);
}

/** Functions for `do_freespace_utetra`: */

static void get_update_fan(eik3_s const *eik, size_t l0, array_s *l_arr) {
  /* Find all of the cells incident on `l0` */
  size_t nvc = mesh3_nvc(eik->mesh, l0);
  size_t *vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(eik->mesh, l0, vc);

  /* Iterate over each cell incident on `l0` */
  for (size_t i = 0; i < nvc; ++i) {
    size_t cv[4];
    mesh3_cv(eik->mesh, vc[i], cv);

    /* Find the current cell's `VALID` vertices */
    size_t num_valid = 0, l_valid[4];
    for (size_t j = 0; j < 4; ++j)
      if (eik->state[cv[j]] == VALID)
        l_valid[num_valid++] = cv[j];

    /* Only interested in cells on the `VALID` front */
    if (num_valid != 3)
      continue;

    /* Grab the two non-`l0` `VALID` vertices in sorted order */
    size_t l[2];
    for (size_t j = 0, k = 0; j < num_valid; ++j)
      if (l_valid[j] != l0)
        l[k++] = l_valid[j];
    SORT2(l[0], l[1]);

    if (array_contains(l_arr, &l))
      continue;

    array_append(l_arr, &l);
  }

  /* Release the cell array */
  free(vc);
}

static void commit_utetra(eik3_s *eik, size_t l, utetra_s const *utetra) {
  /* TODO: at some point, we may want to look into dealing with the
   * case where the new value is approximately equal to the current
   * value. This is a sign that a caustic has formed. */
  assert(utetra_get_value(utetra) < eik->jet[l].f);

  utetra_get_jet31t(utetra, &eik->jet[l]);

  eik3_set_par(eik, l, utetra_get_parent(utetra));
}

static bool
utetras_bracket_ray_1(utetra_s const *utetra, size_t la,
                      utetra_s **utetra_other, size_t num_utetra_other) {
  (void)utetra;

  size_t i_current = 0;

  size_t l_other[2] = {NO_INDEX, NO_INDEX};
  utetra_get_other_inds(utetra_other[i_current], la, l_other);

  size_t l_start = l_other[0], l_current = l_other[1];

  /* Next, we walk around the ring of utetra until we've visited all
   * of them... */
  size_t num_visited = 1;
  while (num_visited++ < num_utetra_other) {
    /* Search through the utetra for the one we should step to next */
    for (size_t i = 0; i < num_utetra_other; ++i) {
      if (i == i_current)
        continue;

      utetra_get_other_inds(utetra_other[i], la, l_other);
      if (l_current != l_other[0] && l_current != l_other[1])
        continue;

      /* Update l_current and i_current and break */
      l_current = l_current == l_other[0] ? l_other[1] : l_other[0];
      i_current = i;
      break;

      /* TODO: actually check whether these vertices form a ring
       * around the active index after projecting onto the normal
       * plane */
    }
  }

  /* We're done walking... check if we're back where we started */
  return l_start == l_current;
}

static void get_unit_normal_for_verts_by_index(mesh3_s const *mesh,
                                               uint3 const l, dbl3 n) {
  dbl3 x0;
  mesh3_copy_vert(mesh, l[0], x0);

  dbl3 dx[2];
  for (size_t i = 0; i < 2; ++i)
    dbl3_sub(mesh3_get_vert_ptr(mesh, l[i + 1]), x0, dx[i]);

  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
}

static bool
utetras_bracket_ray_2(eik3_s const *eik, utetra_s const *utetra,
                      utetra_s **utetra_other, size_t num_utetra_other) {
  mesh3_s const *mesh = eik->mesh;

  size_t l = utetra_get_l(utetra);

  par3_s par = utetra_get_parent(utetra);
  uint3 la, li;
  size_t na = par3_get_active_and_inactive_inds(&par, la, li);
  assert(na == 2);
  SORT_UINT3(la);
  SORT_UINT3(li);

  dbl3 xa, xi;
  mesh3_copy_vert(mesh, la[0], xa);
  mesh3_copy_vert(mesh, li[0], xi);

  for (size_t i = 0; i < num_utetra_other; ++i) {
    assert(l == utetra_get_l(utetra_other[i]));

    par3_s par_ = utetra_get_parent(utetra_other[i]);

    uint3 la_, li_;
    size_t na_ = par3_get_active_and_inactive_inds(&par_, la_, li_);
    SORT_UINT3(la_);
    SORT_UINT3(li_);

    /* Check if both `utetra` have the same active indices */
    if (na_ != 2 && la_[0] != la[0] && la_[1] != la[1])
      continue;

    /* Get the unit normal for the plane determined by the update
     * index and the active edge */
    dbl3 n;
    get_unit_normal_for_verts_by_index(mesh, (uint3) {l, la[0], la[1]}, n);

    /* Get the inactive vertex for the current other utetra */
    dbl3 xi_;
    mesh3_copy_vert(mesh, li_[0], xi_);

    /* Check whether the two inactive vertices lie on either side of
     * the plane */
    dbl n_dot_xa = dbl3_dot(n, xa);
    dbl dp = dbl3_dot(n, xi) - n_dot_xa;
    dbl dp_ = dbl3_dot(n, xi_) - n_dot_xa;
    if (sgn(dp) != sgn(dp_))
      return true;
  }

  return false;
}

static bool
utetras_bracket_ray(eik3_s const *eik, utetra_s const *utetra,
                    utetra_s **utetra_other, size_t num_utetra_other) {
  if (num_utetra_other == 0)
    return false;

  par3_s par = utetra_get_parent(utetra);
  uint3 la;
  size_t na = par3_get_active_inds(&par, la);
  assert(na == 1 || na == 2);

  if (na == 1)
    return utetras_bracket_ray_1(utetra, la[0], utetra_other, num_utetra_other);
  else if (na == 2)
    return utetras_bracket_ray_2(eik, utetra, utetra_other, num_utetra_other);
  else
    assert(false);
}


static utetra_s **
find_and_delete_cached_utetra(eik3_s *eik, utetra_s const *utetra,
                              size_t *num_utetra_other) {
  size_t l = utetra_get_l(utetra);

  /* First, find the indices of the cached utetra which share the same
   * target node and have the same active indices as `utetra`. */
  array_s *i_arr;
  array_alloc(&i_arr);
  array_init(i_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);
  for (size_t i = 0; i < array_size(eik->old_utetra); ++i) {
    utetra_s const *utetra_other;
    array_get(eik->old_utetra, i, &utetra_other);

    if (l != utetra_get_l(utetra_other))
      continue;

    if (!utetras_have_same_minimizer(utetra, utetra_other, eik))
      continue;

    /* Initially just append here so we can check whether we actually
     * want to return these utetras! */
    array_append(i_arr, &i);
  }

  /* Get the pointers to the utetra we found */
  *num_utetra_other = array_size(i_arr);
  utetra_s **utetra_found = malloc(*num_utetra_other*sizeof(utetra_s *));
  for (size_t j = 0, i; j < array_size(i_arr); ++j) {
    array_get(i_arr, j, &i);
    array_get(eik->old_utetra, i, &utetra_found[j]);
    assert(!utetras_have_same_inds(utetra, utetra_found[j]));
  }

  /* We're good if we've found at least one utetra, and either: 1)
   * we're dealing with an optimum on the interior of an edge, or 2)
   * we have an optimum on a corner and we can "walk" the ring of
   * surrounding utetra */
  bool found = *num_utetra_other > 0 &&
    utetras_bracket_ray(eik, utetra, utetra_found, *num_utetra_other);
  if (found) {
    array_delete_all(eik->old_utetra, i_arr);
  } else {
    free(utetra_found);
    *num_utetra_other = 0;
    utetra_found = NULL;
  }

  array_deinit(i_arr);
  array_dealloc(&i_arr);

  return utetra_found;
}

static bool utetra_cached_already(eik3_s const *eik, utetra_s const *utetra) {
  utetra_s const *utetra_other = NULL;
  for (size_t i = 0; i < array_size(eik->old_utetra); ++i) {
    array_get(eik->old_utetra, i, &utetra_other);
    if (utetras_have_same_inds(utetra, utetra_other))
      return true;
  }
  return false;
}

static bool did_utetra_already(eik3_s const *eik, size_t lhat, uint3 const l) {
  for (size_t j = 0; j < array_size(eik->old_utetra); ++j) {
    utetra_s *utetra;
    array_get(eik->old_utetra, j, &utetra);
    if (utetra_has_inds(utetra, lhat, l))
      return true;
  }
  return false;
}

static bool find_and_do_matched_utetra(eik3_s *eik, size_t lhat,
                                       utetra_s const *utetra) {
  size_t num_interior = utetra_get_num_interior_coefs(utetra);
  assert(num_interior == 0 || num_interior == 2);

  size_t num_utetra_other = 0;
  utetra_s **utetra_other =
    find_and_delete_cached_utetra(eik, utetra, &num_utetra_other);

  if (!utetra_other) {
    assert(num_utetra_other == 0);
    return false;
  }

  commit_utetra(eik, lhat, utetra);
  adjust(eik, lhat);

  for (size_t j = 0; j < num_utetra_other; ++j)
    utetra_dealloc(&utetra_other[j]);
  free(utetra_other);

  return true;
}

void do_utetra(eik3_s *eik, size_t lhat, uint3 const l, par3_s *par) {
  if (did_utetra_already(eik, lhat, l))
    return;

  utetra_spec_s spec = utetra_spec_from_eik_and_inds(eik, lhat, l[0],l[1],l[2]);

  if (par != NULL)
    par3_init_empty(par);

  utetra_s *utetra;
  utetra_alloc(&utetra);
  utetra_init(utetra, &spec);

  if (utetra_is_backwards(utetra, eik))
    goto cleanup;

  if (utetra_is_degenerate(utetra))
    goto cleanup;

  utetra_solve(utetra, /* warm start: */ NULL);

  if (par != NULL)
    *par = utetra_get_parent(utetra);

  if (utetra_get_value(utetra) >= eik->jet[lhat].f)
    goto cleanup;

  if (utetra_ray_is_occluded(utetra, eik))
    goto cleanup;

  if (utetra_has_interior_point_solution(utetra)) {
    commit_utetra(eik, lhat, utetra);
    adjust(eik, lhat);
    goto cleanup;
  }

  if (find_and_do_matched_utetra(eik, lhat, utetra))
    goto cleanup;

  if (!utetra_cached_already(eik, utetra)) {
    array_append(eik->old_utetra, &utetra);
    return; /* don't dealloc! */
  }

cleanup:
  utetra_dealloc(&utetra);
}

static void do_1pt_update(eik3_s *eik, size_t l, size_t l0) {
  mesh3_s const *mesh = eik->mesh;

  dbl3 x, x0;
  mesh3_copy_vert(mesh, l, x);
  mesh3_copy_vert(mesh, l0, x0);

  jet31t jet;
  dbl3_sub(x, x0, jet.Df);
  jet.f = eik->jet[l0].f + dbl3_normalize(jet.Df);

  if (jet.f >= eik->jet[l].f)
    return;

  eik->jet[l] = jet;

  par3_init_empty(&eik->par[l]);
  eik->par[l].l[0] = l0;
  eik->par[l].b[0] = 1;

  adjust(eik, l);
}

static void do_utetra_fan(eik3_s *eik, size_t lhat, size_t l0) {
  assert(lhat != l0);

  /* Get the fan of `VALID` triangles incident on `l0`. */
  array_s *le_arr;
  array_alloc(&le_arr);
  array_init(le_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);
  get_update_fan(eik, l0, le_arr);

  size_t l[3] = {l0, (size_t)NO_INDEX, (size_t)NO_INDEX};

  /* Array to track which vertices on the rim of the update fan are
   * incident on `VALID` diffracting edges. */
  array_s *l_diff;
  array_alloc(&l_diff);
  array_init(l_diff, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Do them all */
  for (size_t i = 0; i < array_size(le_arr); ++i) {
    array_get(le_arr, i, &l[1]);

    /* Skip degenerate updates */
    if (lhat == l[1] || lhat == l[2])
      continue;

    /* Prospectively do 1-point updates from point sources */
    for (size_t j = 1; j < 3; ++j) {
      if (jet31t_is_point_source(&eik->jet[l[j]])) {
        assert(array_contains(eik->bc_inds, &l[j]));
        do_1pt_update(eik, lhat, l[j]);
        goto cleanup;
      }
    }

    /* Accumulate `VALID` vertices incident on diffracting edges */
    for (size_t j = 1; j < 3; ++j)
      if (mesh3_vert_incident_on_diff_edge(eik->mesh, l[j])
          && !array_contains(l_diff, &l[j]))
        array_append(l_diff, &l[j]);

    do_utetra(eik, lhat, l, /* par: */ NULL);
  }

  /* Do 2-point diffraction updates */
  for (size_t i = 0; i < array_size(l_diff); ++i) {
    size_t l0;
    array_get(l_diff, i, &l0);
    do_utris_if(eik, lhat, l0, eik->old_diff_utri, is_diff_edge);
  }

cleanup:

  array_deinit(l_diff);
  array_dealloc(&l_diff);

  array_deinit(le_arr);
  array_dealloc(&le_arr);
}

static void update(eik3_s *eik, size_t l, size_t l0) {
  if (jet31t_is_point_source(&eik->jet[l0])) {
    assert(array_contains(eik->bc_inds, &l0));
    do_1pt_update(eik, l, l0);
    return;
  }

  bool l0_is_on_diff_edge = mesh3_vert_incident_on_diff_edge(eik->mesh, l0);
  bool l_is_on_diff_edge = mesh3_vert_incident_on_diff_edge(eik->mesh, l);

  /* If `l0` is incident on a diffracting edge, look for corresponding
   * two-point updates to do. Do not do any other types of updates! */
  if (l0_is_on_diff_edge && !l_is_on_diff_edge)
    do_utris_if(eik, l, l0, eik->old_diff_utri, is_diff_edge);

  /* Check whether l0 is a boundary vertex */
  bool l0_is_bdv = l0_is_on_diff_edge || mesh3_bdv(eik->mesh, l0);

  /* If `l` is a boundary point, do two-point updates that are
   * immersed in the boundary. These are updates that can yield
   * "creeping rays". */
  if (l0_is_bdv && mesh3_bdv(eik->mesh, l))
    do_utris_if(eik, l, l0, eik->old_bd_utri, is_valid_front_bde);

  /* Finally, do the fan of tetrahedron updates. */
  do_utetra_fan(eik, l, l0);
}

void update_neighbors(eik3_s *eik, size_t l0) {
  size_t l; // Node l is a neighbor of node l0

  // Get i0's neighboring nodes.
  int nnb = mesh3_nvv(eik->mesh, l0);
  size_t *nb = malloc(nnb*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, nb);

  // Set FAR nodes to TRIAL and insert them into the heap.
  for (int i = 0; i < nnb; ++i) {
    if (eik->state[l = nb[i]] == FAR) {
      eik->state[l] = TRIAL;
      heap_insert(eik->heap, l);
    }
  }

  // Update neighboring nodes.
  for (int i = 0; i < nnb; ++i) {
    if (eik->state[l = nb[i]] == TRIAL) {
      update(eik, l, l0);
      adjust(eik, l); // TODO: we should avoid calling adjust
                      // repeatedly here and above. Instead we should
                      // use a flag to track which nodes actually had
                      // their values change and then adjust before
                      // finally returning from this function.
    }
  }

  free(nb);
}

/* Remove tetrahedron updates targeting `l0` from `utetra_cache`. */
static void purge_utetra(array_s *utetra_cache, size_t l0) {
  utetra_s *utetra;
  for (size_t i = array_size(utetra_cache); i > 0; --i) {
    array_get(utetra_cache, i - 1, &utetra);
    if (utetra_get_l(utetra) == l0) {
      utetra_dealloc(&utetra);
      array_delete(utetra_cache, i - 1);
    }
  }
}

jmm_error_e eik3_step(eik3_s *eik, size_t *l0) {
  /* Get the first node in the heap. It should be `TRIAL`. */
  *l0 = heap_front(eik->heap);
  assert(eik->state[*l0] == TRIAL);

  /* If the eikonal of the newly VALID node isn't finite, something
   * bad happened. Bail. */
  if (!isfinite(eik->jet[*l0].f))
    return JMM_ERROR_RUNTIME_ERROR;

  /* Otherwise, we pop `l0` from the heap and mark it `VALID`. */
  heap_pop(eik->heap);
  eik->state[*l0] = VALID;

  /* Purge cached updates to keep the cache size under control */
  purge_utetra(eik->old_utetra, *l0);
  utri_cache_purge(eik->old_bd_utri, *l0);
  utri_cache_purge(eik->old_diff_utri, *l0);

  update_neighbors(eik, *l0);

  /* Increment the number of nodes that have been accepted, and mark
   * that the `eik->num_accepted`th node was `l0`. */
  eik->accepted[eik->num_accepted++] = *l0;

  return JMM_ERROR_NONE;
}

jmm_error_e eik3_solve(eik3_s *eik) {
  jmm_error_e error;
  size_t l0;
  while (heap_size(eik->heap) > 0)
    if ((error = eik3_step(eik, &l0)) != JMM_ERROR_NONE)
      break;
  return error;
}

bool eik3_brute_force_remaining(eik3_s *eik) {
  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l)
    if (!eik3_is_valid(eik, l))
      array_append(queue, &l);

  size_t max_num_iter = array_size(queue);
  max_num_iter *= max_num_iter;

  size_t it = 0;
  while (!array_is_empty(queue)) {
    size_t l;
    array_pop_front(queue, &l);

    size_t nvf = mesh3_nvf(eik->mesh, l);
    uint3 *vf = malloc(nvf*sizeof(uint3));
    mesh3_vf(eik->mesh, l, vf);

    purge_utetra(eik->old_utetra, l);

    for (size_t i = 0; i < nvf; ++i)
      if (eik->state[vf[i][0]] == VALID &&
          eik->state[vf[i][1]] == VALID &&
          eik->state[vf[i][2]] == VALID)
        do_utetra(eik, l, vf[i], /* par: */ NULL);

    if (isfinite(eik->jet[l].f)) {
      eik->state[l] = VALID;
      eik->accepted[eik->num_accepted++] = l;
    } else {
      array_append(queue, &l);
    }

    free(vf);

    if (++it == max_num_iter)
      break;
  }

  array_deinit(queue);
  array_dealloc(&queue);

  return eik->num_accepted == mesh3_nverts(eik->mesh);
}

bool eik3_is_solved(eik3_s const *eik) {
  return eik->num_accepted == mesh3_nverts(eik->mesh);
}

static void unaccept_nodes(eik3_s *eik, array_s const *l_arr) {
  size_t j = 0;
  for (size_t i = 0; i < eik->num_accepted; ++i) {
    if (array_contains(l_arr, &eik->accepted[i]))
      continue;
    eik->accepted[j++] = eik->accepted[i];
  }
  for (; j < eik->num_accepted; ++j)
    eik->accepted[j] = (size_t)NO_INDEX;
  eik->num_accepted -= array_size(l_arr);
}

static void reset_nodes(eik3_s *eik, array_s const *l_arr) {
  for (size_t i = 0, l; i < array_size(l_arr); ++i) {
    array_get(l_arr, i, &l);
    assert(eik->state[l] == VALID);

    eik->jet[l] = jet31t_make_empty();
    eik->state[l] = FAR;
    par3_init_empty(&eik->par[l]);

    purge_utetra(eik->old_utetra, l);
    utri_cache_purge(eik->old_bd_utri, l);
    utri_cache_purge(eik->old_diff_utri, l);
  }

  unaccept_nodes(eik, l_arr);
}

static void fix_valid_front(eik3_s *eik) {
  array_s *l_arr;
  array_alloc(&l_arr);
  array_init(l_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l) {
    if (eik->state[l] != FAR)
      continue;

    size_t nvv = mesh3_nvv(eik->mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(eik->mesh, l, vv);

    for (size_t i = 0; i < nvv; ++i)
      if (eik->state[vv[i]] == VALID && !array_contains(l_arr, &vv[i]))
        array_append(l_arr, &vv[i]);

    free(vv);
  }

  /* Reinsert these nodes into the heap */
  for (size_t i = 0, l; i < array_size(l_arr); ++i) {
    array_get(l_arr, i, &l);
    eik->state[l] = TRIAL;
    heap_insert(eik->heap, l);
  }

  unaccept_nodes(eik, l_arr);

  array_deinit(l_arr);
  array_dealloc(&l_arr);
}

/* Sort `eik->accepted` topologically, otherwise can't propagate. */
static void fix_accepted(eik3_s *eik) {
  size_t nverts = mesh3_nverts(eik->mesh);

  /* TODO: Make sure all nodes *have* an accepted value... could relax
   * this, but not a high priority... just shift all nodes w/o an
   * accepted value to the end of array first */
  for (size_t l = 0; l < nverts; ++l)
    assert(l_OK(eik->accepted[l]));

  /* Go through and count the children of each node.
   *
   * Note: this is well-defined, since by construction a node can have
   * at most one parent. */
  size_t *num_ch = calloc(nverts, sizeof(size_t));
  for (size_t l_ch = 0; l_ch < nverts; ++l_ch) {
    uint3 l;
    size_t n = par3_get_active_inds(&eik->par[l_ch], l);
    for (size_t i = 0; i < n; ++i)
      ++num_ch[l[i]];
  }

  /* Figure out the offsets for each subrange of children which will
   * be contained by `ch` */
  size_t *ch_offset = malloc((nverts + 1)*sizeof(size_t));
  ch_offset[0] = 0;
  for (size_t l = 0; l < nverts; ++l)
    ch_offset[l + 1] = ch_offset[l] + num_ch[l];

  /* Reset `num_ch`---we'll count them again below */
  memset(num_ch, 0x0, nverts*sizeof(size_t));

  /* We need to keep track of the indegree of each vertex. We'll use
   * this to control when we insert new nodes into the queue. */
  size_t *indegree = malloc(nverts*sizeof(size_t));

  /* Simultaneously recount the children and populate the subranges of
   * `ch` where the child indices will be kept */
  size_t *ch = malloc(ch_offset[nverts]*sizeof(size_t));
  for (size_t l_ch = 0; l_ch < nverts; ++l_ch) {
    uint3 l;
    size_t n = par3_get_active_inds(&eik->par[l_ch], l);
    for (size_t i = 0; i < n; ++i)
      ch[ch_offset[l[i]] + num_ch[l[i]]++] = l_ch;

    /* Set indegree of child node. */
    indegree[l_ch] = n;
  }

  /* Prepare an array which we'll use for a queue backing the
   * topological sort's BFS */
  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Reset the number of accepted nodes. We'll use this as the
   * accumulation index into `eik->accepted`. */
  eik->num_accepted = 0;

  /* States indicating whether our BFS has visited each node */
  state_e *marked = malloc(nverts*sizeof(state_e));
  for (size_t l = 0; l < nverts; ++l)
    marked[l] = FAR;

  /* Initially populate the queue with nodes without a parent and
   * which also have BCs. We also mark these nodes as `TRIAL` so we
   * can quickly determine whether or not they've been queued. */
  for (size_t l = 0; l < nverts; ++l) {
    if (par3_is_empty(&eik->par[l]) && array_contains(eik->bc_inds, &l)) {
      array_append(queue, &l);
      marked[l] = TRIAL;
    }
  }

  /* Use a BFS to visit each node in topological order */
  while (!array_is_empty(queue)) {
    size_t l;
    array_pop_front(queue, &l);
    assert(indegree[l] == 0);

    /* Mark and accept the current node */
    assert(marked[l] == TRIAL);
    marked[l] = VALID;
    eik->accepted[eik->num_accepted++] = l;

    /* Make sure each of the newly marked node's parents have already
     * been visited */
    uint3 la;
    size_t na = par3_get_active_inds(&eik->par[l], la);
    for (size_t i = 0; i < na; ++i)
      assert(marked[la[i]] == VALID);

    /* Enqueue each child node which hasn't been marked and which
     * isn't already queued */
    for (size_t i = 0, l_ch; i < num_ch[l]; ++i) {
      l_ch = ch[ch_offset[l] + i];

      /* Skip each node with zero indegree... */
      if (indegree[l_ch] == 0)
        continue;

      /* ... and if we don't skip it, it shouldn't be marked yet */
      assert(marked[l_ch] == FAR);

      /* Decrement the indegree of the child node and insert it if it
       * reached zero */
      if (--indegree[l_ch] == 0) {
        array_append(queue, &l_ch);
        marked[l_ch] = TRIAL;
      }
    }
  }

  /* Make sure we visited every node. */
  assert(eik->num_accepted == nverts);
  for (size_t l = 0; l < nverts; ++l) {
    assert(indegree[l] == 0);
    assert(marked[l]);
  }

  /** Clean everything up: */

  free(marked);

  free(indegree);

  array_deinit(queue);
  array_dealloc(&queue);

  free(ch);
  free(ch_offset);
  free(num_ch);
}

void eik3_resolve_downwind_from_diff(eik3_s *eik, size_t diff_index, dbl rfac) {
  mesh3_s const *mesh = eik->mesh;

  /** First, find all of the nodes which are downwind of the
   * diffractor and whose origin value is less than 0.5. These are the
   * nodes that we'll reset, since their error will be polluted. */

  /* Get the current diffractor */
  size_t diff_size = mesh3_get_diffractor_size(mesh, diff_index);
  size_t (*le)[2] = malloc(diff_size*sizeof(size_t[2]));
  mesh3_get_diffractor(mesh, diff_index, le);

  /* Array of unique diffractor node indices */
  array_s *l_diff;
  array_alloc(&l_diff);
  array_init(l_diff, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* ... fill it */
  for (size_t i = 0; i < diff_size; ++i)
    for (size_t j = 0; j < 2; ++j)
      if (!array_contains(l_diff, &le[i][j]))
        array_append(l_diff, &le[i][j]);

  /* Free the diffractor edges */
  free(le);

  /* Queue backing the BFS */
  array_s *l_queue;
  array_alloc(&l_queue);
  array_init(l_queue, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* ... fill it */
  for (size_t i = 0; i < array_size(l_diff); ++i)
    array_append(l_queue, array_get_ptr(l_diff, i));

  /* Get the origins */
  dbl *org = malloc(mesh3_nverts(mesh)*sizeof(dbl));
  eik3_init_org_from_BCs(eik, org);
  eik3_prop_org(eik, org);

  /* Array used to accumulate indices of nodes to be reset */
  array_s *l_reset;
  array_alloc(&l_reset);
  array_init(l_reset, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Use a BFS to find each of the nodes downwind from the diffracting
   * edge which should be reset. */
  while (!array_is_empty(l_queue)) {
    size_t l;
    array_pop_front(l_queue, &l);
    if (array_contains(l_reset, &l))
      continue;

    if (!array_contains(l_diff, &l))
      array_append(l_reset, &l);

    /* Get the neighbors of the current node */
    size_t nvv = mesh3_nvv(mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(mesh, l, vv);

    /* For each neighbor... */
    for (size_t j = 0; j < nvv; ++j) {
      /* ... skip this node if its origin is too big or if we're
       * already resetting it. */
      if (org[vv[j]] >= 0.5 || array_contains(l_reset, &vv[j]))
        continue;

      /* Get the active parent indices of the current node. */
      size_t la[3];
      size_t na = par3_get_active_inds(&eik->par[vv[j]], la);

      for (size_t k = 0; k < na; ++k) {
        if (array_contains(l_diff, &la[k]) || array_contains(l_reset, &la[k])) {
          array_append(l_queue, &vv[j]);
          break;
        }
      }
    }

    free(vv);
  }

  reset_nodes(eik, l_reset);

  /** Now, add BCs for the diffractor and solve again */

  eik3_add_diff_bcs(eik, /* eik_in: */ eik, diff_index, rfac);
  fix_valid_front(eik);
  eik3_solve(eik);
  fix_accepted(eik);

  /** Cleanup */

  array_deinit(l_reset);
  array_dealloc(&l_reset);

  free(org);

  array_deinit(l_queue);
  array_dealloc(&l_queue);

  array_deinit(l_diff);
  array_dealloc(&l_diff);
}

void eik3_add_trial(eik3_s *eik, size_t l, jet31t jet) {
  if (eik->state[l] == VALID) {
    log_warn("failed to add TRIAL node %lu (already VALID)", l);
    return;
  }

  if (isfinite(eik->jet[l].f)) {
    log_warn("failed to add TRIAL node %lu (finite jet)", l);
    return;
  }

  assert(eik->pos[l] == NO_INDEX);
  assert(eik->state[l] == FAR);

  eik->jet[l] = jet;
  eik->state[l] = TRIAL;
  heap_insert(eik->heap, l);

  array_append(eik->trial_inds, &l);
}

size_t eik3_num_trial(eik3_s const *eik) {
  size_t num_trial = 0;
  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l)
    num_trial += eik->state[l] == TRIAL;
  return num_trial;
}

size_t eik3_num_valid(eik3_s const *eik) {
  return eik->num_accepted;
}

void eik3_add_bc(eik3_s *eik, size_t l, jet31t jet) {
  assert(!array_contains(eik->bc_inds, &l));

  eik->jet[l] = jet;
  eik->state[l] = VALID;
  eik->accepted[eik->num_accepted++] = l;

  array_append(eik->bc_inds, &l);
}

static bool has_nb_with_state(eik3_s const *eik, size_t l, state_e state) {
  size_t nvv = mesh3_nvv(eik->mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(eik->mesh, l, vv);

  bool has_nb = false;

  for (size_t i = 0; i < nvv; ++i) {
    if (eik->state[vv[i]] == state) {
      has_nb = true;
      break;
    }
  }

  free(vv);

  return has_nb;
}

static void freeze_bc_layer(eik3_s *eik) {
  array_s *l_arr;
  array_alloc(&l_arr);
  array_init(l_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  while (heap_size(eik->heap) > 0) {
    size_t l = heap_front(eik->heap);
    heap_pop(eik->heap);
    array_append(l_arr, &l);
  }

  for (size_t i = 0, l; i < array_size(l_arr); ++i) {
    array_get(l_arr, i, &l);

    if (has_nb_with_state(eik, l, FAR)) {
      heap_insert(eik->heap, l);
      continue;
    }

    eik->state[l] = VALID;
    eik->accepted[eik->num_accepted++] = l;
  }

  array_deinit(l_arr);
  array_dealloc(&l_arr);
}

void eik3_add_pt_src_bcs(eik3_s *eik, dbl3 const xsrc, dbl rfac) {
  mesh3_s const *mesh = eik->mesh;

  /* We assume that `xsrc` is a mesh vertex. */
  assert(mesh3_has_vertex(mesh, xsrc));

  /* Get the index of the vertex corresponding to `xsrc`. */
  size_t lsrc = mesh3_get_vert_index(mesh, xsrc);

  eik3_add_bc(eik, lsrc, jet31t_make_point_source(0));

  /** BFS starting from `xsrc`: */

  /* The queue driving the BFS */
  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Put `xsrc` into `queue` initially. */
  array_append(queue, &lsrc);

  /* Main BFS loop: */
  while (!array_is_empty(queue)) {
    size_t l;
    array_pop_front(queue, &l);

    /* Set the node to trial and insert it into the heap */
    if (eik3_is_far(eik, l))
      eik3_add_trial(eik, l, jet31t_make_empty());

    /* Do one-point updates from `xsrc`. */
    if (l != lsrc)
      do_1pt_update(eik, l, lsrc);
    assert(isfinite(eik->jet[l].f));

    /* Check if we're in the factoring ball... */
    if (eik->jet[l].f > rfac)
      continue;

    /* ... if we are, then add this node's neighbors to the queue. */
    size_t nvv = mesh3_nvv(mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(mesh, l, vv);
    for (size_t i = 0; i < nvv; ++i)
      if (vv[i] != lsrc && isinf(eik->jet[vv[i]].f))
        array_append(queue, &vv[i]);
    free(vv);
  }

  freeze_bc_layer(eik);

  /* Make sure we added some boundary data: */
  assert(!array_is_empty(eik->bc_inds));
  assert(!array_is_empty(eik->trial_inds));
}

/* Store the cubic polynomial `T` approximating the eikonal over this
 * diffracting edge. */
static void store_diff_bc_T(eik3_s *eik, size_t const le[2], bb31 const *T) {
  size_t key[2] = {le[0], le[1]};
  SORT2(key[0], key[1]);
  alist_append(eik->T_diff, key, T);
}

static void add_diff_bc_for_edge_from_bb31(eik3_s *eik, size_t const le[2],
                                           bb31 const *T){
  /* Add BCs for each node with values taken from T... but with
   * a singular gradient! */
  jet31t jet[2] = {
    jet31t_make_point_source(T->c[0]),
    jet31t_make_point_source(T->c[3])
  };

  for (size_t i = 0; i < 2; ++i)
    if (eik->state[le[i]] == FAR)
      eik3_add_bc(eik, le[i], jet[i]);

  store_diff_bc_T(eik, le, T);
}

bool eik3_has_diff_bc(eik3_s const *eik, size_t const le[2]) {
  size_t key[2] = {le[0], le[1]};
  SORT2(key[0], key[1]);
  return alist_contains(eik->T_diff, key);
}

void eik3_get_diff_bc(eik3_s const *eik, size_t const le[2], bb31 *T) {
  size_t key[2] = {le[0], le[1]};
  SORT2(key[0], key[1]);
  alist_get_by_key(eik->T_diff, key, T);
}

static dbl get_par_eik_value(eik3_s const *eik, size_t l) {
  par3_s const *par = &eik->par[l];

  if (par3_is_empty(par))
    return eik->jet[l].f;

  size_t la[3];
  dbl3 ba;
  size_t na = par3_get_active(par, la, ba);

  dbl T = 0;
  for (size_t i = 0; i < na; ++i)
    T += ba[i]*eik->jet[la[i]].f;

  return T;
}

typedef struct update_inds {
  size_t lhat;
  size_t l[3];
} update_inds_s;

update_inds_s make_empty_update_inds(size_t lhat) {
  return (update_inds_s) {
    .lhat = lhat,
    .l = {(size_t)NO_INDEX, (size_t)NO_INDEX, (size_t)NO_INDEX}
  };
}

/* Set up BCs for each diffracting edge for this diffractor */
static void add_diff_bcs_for_diffractor(eik3_s *eik, eik3_s const *eik_in,
                                        size_t diff_index) {
  assert(eik->mesh == eik_in->mesh);
  mesh3_s const *mesh = eik->mesh;

  size_t num_diff_edges = mesh3_get_diffractor_size(mesh, diff_index);
  size_t (*le)[2] = malloc(num_diff_edges*sizeof(size_t[2]));
  mesh3_get_diffractor(mesh, diff_index, le);

  for (size_t i = 0; i < num_diff_edges; ++i) {
    jet31t jet[2] = {eik_in->jet[le[i][0]], eik_in->jet[le[i][1]]};

    dbl3 x[2];
    mesh3_copy_vert(mesh, le[i][0], x[0]);
    mesh3_copy_vert(mesh, le[i][1], x[1]);

    bb31 T;
    bb31_init_from_jets(&T, jet, x);

    add_diff_bc_for_edge_from_bb31(eik, le[i], &T);
  }

  free(le);
}

static bool OK_edge_inds(uint2 const le) {
  return le[0] != (size_t)NO_INDEX && le[1] != (size_t)NO_INDEX;
}

static size_t add_next_edge_to_queue(eik3_s *eik, mesh1_s const *diff_mesh,
                                     size_t lhat, uint2 const l, size_t la,
                                     array_s *queue) {
  size_t le[2];
  mesh1_ve(diff_mesh, la, le);
  assert(le[0] != (size_t)NO_INDEX || le[1] != (size_t)NO_INDEX);

  uint2 l_sorted = {l[0], l[1]};
  SORT2(l_sorted[0], l_sorted[1]);

  bool added = false;

  for (size_t i = 0; i < 2; ++i) {
    if (le[i] == (size_t)NO_INDEX)
      continue;

    uint2 ev;
    mesh1_ev(diff_mesh, le[i], ev);
    SORT2(ev[0], ev[1]);

    if (l_sorted[0] == ev[0] && l_sorted[1] == ev[1])
      continue;

    if (array_contains(queue, &ev) ||
        utri_cache_contains_inds(eik->old_diff_utri, lhat, ev))
      continue;

    array_append(queue, &ev);
    added = true;
  }

  assert(added);

  return 1;
}

static void do_diffractor_updates(eik3_s *eik, mesh1_s const *diff_mesh,
                                  size_t l, uint2 const le_guess) {
  assert(OK_edge_inds(le_guess));

  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(uint2), ARRAY_DEFAULT_CAPACITY);

  uint2 le = {le_guess[0], le_guess[1]};
  SORT2(le[0], le[1]);
  array_append(queue, &le);

  while (isinf(eik->jet[l].f)) {
    assert(!array_is_empty(queue));
    array_pop_front(queue, &le);

    par3_s par;

    /* Do the current triangle update */
    assert(!utri_cache_contains_inds(eik->old_diff_utri, l, le));
    do_utri(eik, l, le[0], le[1], eik->old_diff_utri, &par);

    /* If we managed to update `l`, break early */
    if (isfinite(eik->jet[l].f))
      break;

    /* TODO: for now... */
    assert(!par3_is_empty(&par));

    /* Get the active indices */
    size_t la[3];
    size_t na = par3_get_active_inds(&par, la);
    assert(na == 1);

    size_t num_added = add_next_edge_to_queue(eik,diff_mesh,l,le,la[0],queue);
    assert(num_added == 1);
  }

  array_deinit(queue);
  array_dealloc(&queue);
}

update_inds_s get_diff_update_inds(mesh1_s const *diff_mesh,
                                   update_inds_s parent_update_inds, size_t l) {
  update_inds_s update_inds = parent_update_inds;
  update_inds.lhat = l;

  size_t num_active = 0;
  for (size_t i = 0; i < 3; ++i)
    num_active += update_inds.l[i] != (size_t)NO_INDEX;
  assert(num_active == 0 || num_active == 2);

  if (num_active == 0) {
    size_t le[2];
    mesh1_ve(diff_mesh, parent_update_inds.lhat, le);

    bool le_active[2] = {
      le[0] != (size_t)NO_INDEX,
      le[1] != (size_t)NO_INDEX
    };

    assert(le_active[0] || le_active[1]);

    mesh1_ev(diff_mesh, le_active[0] ? le[0] : le[1], update_inds.l);
  }

  return update_inds;
}

/* Add edge diffraction BCs in a tube surrounding diff_index */
void eik3_add_diff_bcs(eik3_s *eik, eik3_s const *eik_in, size_t diff_index, dbl rfac) {
  assert(eik->mesh == eik_in->mesh);
  mesh3_s const *mesh = eik3_get_mesh(eik);

  /* Add BCs for the diffractor */
  add_diff_bcs_for_diffractor(eik, eik_in, diff_index);

  /* Get the polygonal curve corresponding to the diffractor. We'll
   * use this to do a local search for updates to accelerate our
   * BFS. */
  mesh1_s *diff_mesh = mesh3_get_diff_mesh(mesh, diff_index);

  /* The queue driving the BFS */
  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(update_inds_s), ARRAY_DEFAULT_CAPACITY);

  /* Seed the BFS's queue with updates for each of the nodes on the
   * diffractor */
  for (size_t le = 0; le < mesh1_nedges(diff_mesh); ++le) {
    uint2 l;
    mesh1_ev(diff_mesh, le, l);
    for (size_t i = 0; i < 2; ++i) {
      update_inds_s update_inds = make_empty_update_inds(l[i]);
      if (!array_contains(queue, &update_inds))
        array_append(queue, &update_inds);
    }
  }

  /* Main BFS loop: */
  while (!array_is_empty(queue)) {
    update_inds_s update_inds;
    array_pop_front(queue, &update_inds);

    size_t l = update_inds.lhat, *le = &update_inds.l[0];

    /* If `l` isn't one of the initial points on the diffractor, set
     * the node to trial, insert it into the heap, and do updates from
     * the reflector, using `le` as a warm start. */
    if (OK_edge_inds(le) && eik3_is_far(eik, l)) {
      eik3_add_trial(eik, l, jet31t_make_empty());
      do_diffractor_updates(eik, diff_mesh, l, le);
    }

    /* Check if we're in the "factoring tube"... */
    if (eik->jet[l].f - get_par_eik_value(eik, l) > rfac)
      continue;

    /* ... if we are, then add this node's neighbors to the queue,
     * using the parent of `l` as a warm start for the reflector
     * updates for each added node. */

    size_t nvv = mesh3_nvv(mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(mesh, l, vv);

    for (size_t i = 0; i < nvv; ++i) {
      /* Skip this node if we've already updated it */
      if (isfinite(eik->jet[vv[i]].f))
        continue;

      /* Get the update indices for the child and append it to the
       * queue if we haven't done so already */
      update_inds_s child_update_inds =
        get_diff_update_inds(diff_mesh, update_inds, vv[i]);
      if (!array_contains(queue, &child_update_inds))
        array_append(queue, &child_update_inds);
    }

    free(vv);
  }

  freeze_bc_layer(eik);

  array_deinit(queue);
  array_dealloc(&queue);
}

static bool OK_face_inds(uint3 const lf) {
  return lf[0] != (size_t)NO_INDEX
    && lf[1] != (size_t)NO_INDEX && lf[2] != (size_t)NO_INDEX;
}

static bool add_face_to_queue(eik3_s const *eik, mesh2_s const *refl_mesh,
                              size_t lhat, size_t lf, array_s *queue) {
  uint3 l;
  mesh2_fv(refl_mesh, lf, l);
  SORT3(l[0], l[1], l[2]);

  if (array_contains(queue, &l) || did_utetra_already(eik, lhat, l))
    return false;

  array_append(queue, &l);
  return true;
}

static size_t add_vf_to_queue(eik3_s *eik, mesh2_s const *refl_mesh,
                              size_t lhat, size_t la, array_s *queue) {
  size_t nvf = mesh2_nvf(refl_mesh, la);
  size_t *lf = malloc(nvf*sizeof(uint3));
  mesh2_vf(refl_mesh, la, lf);

  size_t num_added = 0;
  for (size_t i = 0; i < nvf; ++i)
    if (add_face_to_queue(eik, refl_mesh, lhat, lf[i], queue))
      ++num_added;

  free(lf);

  return num_added;
}

static size_t add_opposite_face_to_queue(eik3_s *eik, mesh2_s const *refl_mesh,
                                         size_t lhat, uint3 const vf,
                                         size_t la[2], array_s *queue) {
  /* Find the inactive face vertex */
  size_t l = (size_t)NO_INDEX;
  for (size_t i = 0; i < 3; ++i)
    if (vf[i] != la[0] && vf[i] != la[1])
      l = vf[i];
  assert(l != (size_t)NO_INDEX);

  /* Get the current face index */
  size_t lf = mesh2_find_face(refl_mesh, vf);

  /* Find the face opposite from the inactive vertex */
  size_t lf_op = mesh2_fvf(refl_mesh, lf, l);
  if (lf_op == (size_t)NO_INDEX)
    return 0;

  /* And add it to the queue */
  return add_face_to_queue(eik, refl_mesh, lhat, lf_op, queue) ? 1 : 0;
}

static void do_utri_and_add_inc(eik3_s *eik, mesh2_s const *refl_mesh,
                                size_t lhat, uint3 const l, array_s *queue) {
  mesh3_s const *mesh = eik->mesh;
  assert(mesh3_is_diff_edge(mesh, l));

  par3_s par;
  do_utri(eik, lhat, l[0], l[1], eik->old_diff_utri, &par);
  assert(!par3_is_empty(&par));
  if (isfinite(eik->jet[lhat].f))
    return;

  size_t la[3];
  size_t na = par3_get_active_inds(&par, la);
  assert(na == 1);

  size_t nvf = mesh2_nvf(refl_mesh, la[0]);
  size_t *vf = malloc(nvf*sizeof(size_t));
  mesh2_vf(refl_mesh, la[0], vf);

  for (size_t i = 0; i < nvf; ++i) {
    uint3 fv;
    mesh2_fv(refl_mesh, vf[i], fv);
    SORT_UINT3(fv);
    if (!array_contains(queue, &fv) && !did_utetra_already(eik, lhat, fv))
      array_append(queue, &fv);
  }

  for (size_t i = 0; i < nvf; ++i) {
    uint3 fv;
    mesh2_fv(refl_mesh, vf[i], fv);
    for (size_t j = 0; j < 3; ++j) {
      uint3 ve_ = {fv[j], fv[(j + 1) % 3], (size_t)NO_INDEX};
      SORT_UINT3(ve_);
      if (ve_[0] != la[0] && ve_[1] != la[0])
        continue;
      if (mesh3_is_diff_edge(mesh, ve_)
          && !array_contains(queue, &ve_)
          && !utri_cache_contains_inds(eik->old_diff_utri, lhat, ve_))
        array_append(queue, &ve_);
    }
  }

  free(vf);
}

static void add_diff_utri_inc_on_utetra(eik3_s *eik, size_t lhat,
                                        uint3 const l, array_s *queue) {
  mesh3_s const *mesh = eik->mesh;
  for (size_t i = 0; i < 3; ++i) {
    size_t ne = mesh3_get_num_inc_diff_edges(mesh, l[i]);
    uint2 *le_diff = malloc(ne*sizeof(uint2));
    mesh3_get_inc_diff_edges(mesh, l[i], le_diff);
    for (size_t j = 0; j < ne; ++j) {
      uint3 le = {le_diff[j][0], le_diff[j][1], (size_t)NO_INDEX};
      SORT_UINT2(le);
      if (!array_contains(queue, &le)
          && !utri_cache_contains_inds(eik->old_diff_utri, lhat, le))
        array_append(queue, &le);
    }
    free(le_diff);
  }
}

static void do_utetra_and_add_inc(eik3_s *eik, mesh2_s const *refl_mesh,
                                  size_t lhat, uint3 const l, array_s *queue) {
  par3_s par;
  do_utetra(eik, lhat, l, &par);
  if (par3_is_empty(&par)) {
    add_diff_utri_inc_on_utetra(eik, lhat, l, queue);
    return;
  } else if (isfinite(eik->jet[lhat].f))
    return;

  size_t la[3];
  size_t na = par3_get_active_inds(&par, la);
  assert(na < 3);

  size_t num_added = 0;

  /* If there's only one active index, find the faces incident on
   * the active vertex */
  if (na == 1)
    num_added = add_vf_to_queue(eik, refl_mesh, lhat, la[0], queue);

  /* If there are two active indices, add the opposite face */
  if (na == 2)
    num_added = add_opposite_face_to_queue(eik, refl_mesh, lhat, l, la, queue);

  if (num_added == 0 && na == 2 && mesh3_is_diff_edge(eik->mesh, la)) {
    uint3 la_ = {la[0], la[1], (size_t)NO_INDEX};
    if (!array_contains(queue, &la_)
        && !utri_cache_contains_inds(eik->old_diff_utri, lhat, la_))
      array_append(queue, &la_);
  }

  /* This case should probably never happen... why?
   *
   * 1) this is the case where there's a single active index incident
   * on a diffracting *VERTEX*, and we haven't added any new faces

   * 2) but if this happens, then we should have already done the
   * incident tetrahedron updates
   *
   * 3) ... in which case, we should be able to accept a "common
   * vertex" update, but we failed to do so for some reason!
   *
   * TODO: if this happens, it's a bug, I think */
  if (num_added == 0
      && na == 1
      && mesh3_vert_is_terminal_diff_edge_vert(eik->mesh, la[0])
      && mesh3_get_num_inc_diff_edges(eik->mesh, la[0]) > 1)
    assert(false); // TODO: handle
}

static void do_reflector_updates(eik3_s *eik, mesh2_s const *refl_mesh,
                                 size_t lhat, uint3 const l_guess) {
  assert(OK_face_inds(l_guess));

  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(uint3), ARRAY_DEFAULT_CAPACITY);

  uint3 l = {l_guess[0], l_guess[1], l_guess[2]};
  SORT3(l[0], l[1], l[2]);
  array_append(queue, &l);

  while (isinf(eik->jet[lhat].f) && !array_is_empty(queue)) {
    array_pop_front(queue, &l);

    if (is_edge(l))
      do_utri_and_add_inc(eik, refl_mesh, lhat, l, queue);
    else if (is_face(l))
      do_utetra_and_add_inc(eik, refl_mesh, lhat, l, queue);
  }

  array_deinit(queue);
  array_dealloc(&queue);
}

static update_inds_s get_refl_update_inds(mesh2_s const *refl_mesh,
                                          update_inds_s parent_update_inds,
                                          size_t l) {
  update_inds_s update_inds = parent_update_inds;
  update_inds.lhat = l;

  size_t num_active = 0;
  for (size_t i = 0; i < 3; ++i)
    num_active += update_inds.l[i] != (size_t)NO_INDEX;
  assert(num_active == 0 || num_active == 3);

  if (num_active == 0) {
    size_t nvf = mesh2_nvf(refl_mesh, parent_update_inds.lhat);
    size_t *vf = malloc(nvf*sizeof(size_t));
    mesh2_vf(refl_mesh, parent_update_inds.lhat, vf);

    assert(nvf > 0);

    mesh2_fv(refl_mesh, vf[0], update_inds.l);

    free(vf);
  }

  return update_inds;
}

/* Setup routine common to both `eik3_add_refl_bcs` and
 * `eik3_add_refl_bcs_with_fac`. */
static mesh2_s *
init_refl_bcs(eik3_s *eik, eik3_s const *eik_in, size_t refl_index) {
  /* Both `eik` and `eik_in` should have the same domain. */
  assert(eik->mesh == eik_in->mesh);
  mesh3_s const *mesh = eik->mesh;

  /* Get the reflector mesh. We'll use this to do a local search for
   * updates to accelerate our BFS. */
  mesh2_s *refl_mesh = mesh3_get_refl_mesh(mesh, refl_index);

  /* Add reflection BCs for each vertex on the reflector. */
  for (size_t lf = 0; lf < mesh2_nfaces(refl_mesh); ++lf) {
    dbl33 R;
    mesh2_get_R_for_face(refl_mesh, lf, R);

    uint3 l;
    mesh2_fv(refl_mesh, lf, l);

    /* Add BCs and insert each node into the update queue.*/
    for (size_t i = 0; i < 3; ++i) {
      if (eik3_has_BCs(eik, l[i]))
        continue;

      /* ... add the BCs... */
      jet31t jet = eik3_get_jet(eik_in, l[i]);
      dbl33_dbl3_mul_inplace(R, jet.Df);
      eik3_add_bc(eik, l[i], jet);
    }
  }

  /* Add BCs for each diffracting edge incident on the reflector. */
  size_t num_diffractors = mesh3_get_num_diffs_inc_on_refl(mesh, refl_index);
  size_t *diff_index = malloc(num_diffractors*sizeof(size_t));
  mesh3_get_diffs_inc_on_refl(mesh, refl_index, diff_index);
  for (size_t i = 0; i < num_diffractors; ++i)
    add_diff_bcs_for_diffractor(eik, eik_in, diff_index[i]);
  free(diff_index);

  return refl_mesh;
}

void eik3_add_refl_bcs(eik3_s *eik, eik3_s const *eik_in, size_t refl_index) {
  /* Both `eik` and `eik_in` should have the same domain. */
  assert(eik->mesh == eik_in->mesh);
  mesh3_s const *mesh = eik->mesh;

  /* Get the reflector mesh. We'll use this to do a local search for
   * updates to accelerate our BFS. */
  mesh2_s *refl_mesh = mesh3_get_refl_mesh(mesh, refl_index);

  /* Add reflection BCs for each vertex on the reflector. */
  for (size_t lf = 0; lf < mesh2_nfaces(refl_mesh); ++lf) {
    dbl33 R;
    mesh2_get_R_for_face(refl_mesh, lf, R);

    uint3 l;
    mesh2_fv(refl_mesh, lf, l);

    /* Add BCs and insert each node into the update queue.*/
    for (size_t i = 0; i < 3; ++i) {
      if (!eik3_is_far(eik, l[i])) continue;
      jet31t jet = eik3_get_jet(eik_in, l[i]);
      dbl33_dbl3_mul_inplace(R, jet.Df);
      eik3_add_trial(eik, l[i], jet);
    }
  }
}

/* Add reflection BCs using a factoring radius. Each node reachable by
 * a ray of length at most `rfac` will be updated using a BFS-based
 * update algorithm over the reflector itself.
 *
 * CAUTION: I have not looked into how this copes with visibility
 * issues. If there are occluding obstacles inside the factoring
 * region, bad things could happen. */
void eik3_add_refl_bcs_with_fac(eik3_s *eik, eik3_s const *eik_in,
                                size_t refl_index, dbl rfac) {
  /* Add reflection BCs. */
  mesh2_s *refl_mesh = init_refl_bcs(eik, eik_in, refl_index);

  /* The queue driving the BFS. We store the node to update along with
   * a starting guess for a face to start with. */
  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(update_inds_s), ARRAY_DEFAULT_CAPACITY);

  /* Add the BC nodes into the queue backing the BFS. */
  for (size_t lf = 0; lf < mesh2_nfaces(refl_mesh); ++lf) {
    uint3 l;
    mesh2_fv(refl_mesh, lf, l);
    for (size_t i = 0; i < 3; ++i) {
      if (eik3_has_BCs(eik, l[i])) {
        update_inds_s update_inds = make_empty_update_inds(l[i]);
        if (!array_contains(queue, &update_inds))
          array_append(queue, &update_inds);
      }
    }
  }

  /* Main BFS loop: */
  while (!array_is_empty(queue)) {
    update_inds_s update_inds;
    array_pop_front(queue, &update_inds);

    size_t l = update_inds.lhat, *lf = &update_inds.l[0];

    /* If `l` is off the reflector, set the node to trial, insert it
     * into the heap, and do updates from the reflector, using `lf` as
     * a warm start. */
    if (OK_face_inds(lf) && eik3_is_far(eik, l)) {
      eik3_add_trial(eik, l, jet31t_make_empty());
      do_reflector_updates(eik, refl_mesh, l, lf);
    }

    /* Check if we're in the "factoring tube"... */
    if (eik->jet[l].f - get_par_eik_value(eik, l) > rfac)
      continue;

    /* Get `l`'s neighbors */
    size_t nvv = mesh3_nvv(eik->mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(eik->mesh, l, vv);

    /* Since we're in the "factoring tube" now, add this node's
     * neighbors to the queue, using the parent of `l` as a warm start
     * for the reflector updates for each added node. */
    for (size_t i = 0; i < nvv; ++i) {
      /* Skip this node if we've already updated it */
      if (isfinite(eik->jet[vv[i]].f))
        continue;

      update_inds_s child_update_inds =
        get_refl_update_inds(refl_mesh, update_inds, vv[i]);
      if (!array_contains(queue, &child_update_inds))
        array_append(queue, &child_update_inds);
    }

    free(vv);
  }

  freeze_bc_layer(eik);

  /* Clean up */

  array_deinit(queue);
  array_dealloc(&queue);

  mesh2_deinit(refl_mesh);
  mesh2_dealloc(&refl_mesh);
}

bool eik3_is_far(eik3_s const *eik, size_t l) {
  return eik->state[l] == FAR;
}

bool eik3_is_trial(eik3_s const *eik, size_t l) {
  return eik->state[l] == TRIAL;
}

bool eik3_is_valid(eik3_s const *eik, size_t l) {
  return eik->state[l] == VALID;
}

mesh3_s const *eik3_get_mesh(eik3_s const *eik) {
  return eik->mesh;
}

array_s const *eik3_get_trial_inds(eik3_s const *eik) {
  return eik->trial_inds;
}

array_s const *eik3_get_bc_inds(eik3_s const *eik) {
  return eik->bc_inds;
}

dbl eik3_get_T(eik3_s const *eik, size_t l) {
  return eik->jet[l].f;
}

jet31t eik3_get_jet(eik3_s const *eik, size_t l) {
  return eik->jet[l];
}

void eik3_set_jet(eik3_s *eik, size_t l, jet31t jet) {
  eik->jet[l] = jet;
}

jet31t *eik3_get_jet_ptr(eik3_s const *eik) {
  return eik->jet;
}

state_e *eik3_get_state_ptr(eik3_s const *eik) {
  return eik->state;
}

par3_s eik3_get_par(eik3_s const *eik, size_t l) {
  return eik->par[l];
}

void eik3_set_par(eik3_s *eik, size_t l, par3_s par) {
  eik->par[l] = par;
}

bool eik3_has_par(eik3_s const *eik, size_t l) {
  return !par3_is_empty(&eik->par[l]);
}

bool eik3_has_BCs(eik3_s const *eik, size_t l) {
  return array_contains(eik->bc_inds, &l);
}

size_t const *eik3_get_accepted_ptr(eik3_s const *eik) {
  return eik->accepted;
}

size_t eik3_num_bc(eik3_s const *eik) {
  return array_size(eik->bc_inds);
}

void eik3_get_edge_T(eik3_s const *eik, size_t const le[2], bb31 *T) {
  assert(mesh3_is_edge(eik->mesh, le));

  dbl3 x[2];
  mesh3_copy_vert(eik->mesh, le[0], x[0]);
  mesh3_copy_vert(eik->mesh, le[1], x[1]);

  jet31t jet[2] = {eik->jet[le[0]], eik->jet[le[1]]};

  bb31_init_from_jets(T, jet, x);
}

dbl eik3_get_max_T(eik3_s const *eik) {
  dbl T_max = -INFINITY;
  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l)
    T_max = fmax(T_max, eik->jet[l].f);
  return T_max;
}

static void init_org(mesh3_s const *mesh, dbl *org) {
  /* Initialize `org` to `NAN` */
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    org[l] = NAN;

  // /* Set the origin of all nodes on a diffracting edge to 0 */
  // for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
  //   if (isnan(org[l]) && mesh3_vert_incident_on_diff_edge(mesh, l))
  //     org[l] = 0;
}

void eik3_init_org_from_BCs(eik3_s const *eik, dbl *org) {
  init_org(eik->mesh, org);

  /* Set the origin of all BC nodes to 1 */
  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);
    if (isnan(org[l]))
      org[l] = 1;
  }
}

void eik3_init_org_for_refl(eik3_s const *eik, dbl *org, size_t refl_index,
                            dbl const *org_in) {
  mesh3_s const *mesh = eik->mesh;

  init_org(mesh, org);

  size_t nf = mesh3_get_reflector_size(mesh, refl_index);
  uint3 *lf = malloc(nf*sizeof(uint3));
  mesh3_get_reflector(mesh, refl_index, lf);

  for (size_t i = 0; i < nf; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (isnan(org[lf[i][j]]))
        org[lf[i][j]] = org_in[lf[i][j]];

  free(lf);
}

void eik3_prop_org(eik3_s const *eik, dbl *org) {
  // mesh3_s const *mesh = eik->mesh;

  // /* Now, transport the origins, skipping already set values */
  // eik3_transport_dbl(eik, org, true);

  // /* We set `org` to 0.5 for each vertex on a diffracting edge with
  //  * BCs. */
  // for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
  //   if (mesh3_vert_incident_on_diff_edge(mesh, l)
  //       && array_contains(eik->bc_inds, &l))
  //     org[l] = 0.5;

  // /* We set `org` to 0.5 for each vertex on a diffracting edge that
  //  * was updated from a node with an origin equal to 1. */
  // for (size_t l = 0; l < mesh3_nverts(mesh); ++l) {
  //   if (!mesh3_vert_incident_on_diff_edge(mesh, l))
  //     continue;

  //   par3_s const *par = &eik->par[l];

  //   size_t l_active[3];
  //   size_t num_active = par3_get_active_inds(par, l_active);

  //   for (size_t i = 0; i < num_active; ++i)
  //     if (org[l_active[i]] == 1)
  //       org[l] = 0.5;
  // }

  mesh3_s const *mesh = eik->mesh;
  size_t nverts = mesh3_nverts(mesh);

  bool *diffracting = calloc(nverts, sizeof(bool));

  // for (size_t l = 0; l < nverts; ++l) {
  //   if (!isnan(org[l])) continue;

  //   par3_s const *par = &eik->par[l];
  //   assert(!par3_is_empty(par));

  //   uint3 la;
  //   dbl3 b;
  //   size_t na = par3_get_active(par, la, b);

  //   if (na == 1 && mesh3_vert_incident_on_diff_edge(mesh, la[0])) {
  //     org[la[0]] = 0;
  //     diffracting[la[0]] = true;
  //   }

  //   if (na == 2 && mesh3_is_diff_edge(mesh, la)) {
  //     org[la[0]] = org[la[1]] = 0;
  //     diffracting[la[0]] = diffracting[la[1]] = true;
  //   }
  // }

  // for (size_t l = 0; l < nverts; ++l) {
  //   if (!mesh3_vert_incident_on_diff_edge(mesh, l))
  //     continue;

  //   dbl const *DT = &eik->jet[l].Df[0];

  //   if (!mesh3_local_ray_in_vertex_cone(mesh, DT, l))
  //     continue;

  //   diffracting[l] = true;
  //   org[l] = 0.0;
  // }

  for (size_t l = 0; l < nverts; ++l) {
    if (mesh3_vert_incident_on_diff_edge(mesh, l)) {
      diffracting[l] = true;
      org[l] = 0.0;
    }
  }

  eik3_transport_dbl(eik, org, true);

  for (size_t i = 0; i < nverts; ++i) {
    size_t l = eik->accepted[i];

    if (!diffracting[l]) continue;

    par3_s const *par = &eik->par[l];
    if (par3_is_empty(par))
      continue;

    uint3 la = {NO_INDEX, NO_INDEX, NO_INDEX};
    dbl3 b = {NAN, NAN, NAN};
    size_t na = par3_get_active(par, la, b);
    (void)na;

    for (size_t j = 0; j < na; ++j)
      if (diffracting[la[j]])
        la[j] = NO_INDEX;

    dbl3 orgpar = {NAN, NAN, NAN};
    dbl3_gather(org, la, orgpar);

    dbl nanmin = dbl3_nanmin(orgpar);
    if (nanmin > 0.5)
      org[l] = 0.5;
  }

  free(diffracting);
}

bool eik3_updated_from_diff_edge(eik3_s const *eik, size_t l) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  par3_s par = eik3_get_par(eik, l);

  size_t npar = 0;
  for (size_t i = 0; i < 3; ++i)
    npar += par.l[i] != NO_PARENT;

  if (npar == 0 || npar == 3)
    return false;
  else if (npar == 1)
    return mesh3_vert_incident_on_diff_edge(mesh, par.l[0]);
  else /* npar == 2 */
    return mesh3_is_diff_edge(mesh, par.l);
}

static bool
any_cell_vert_updated_from_diff_edge(eik3_s const *eik, size_t cv[4]) {
  for (size_t i = 0; i < 4; ++i)
    if (eik3_updated_from_diff_edge(eik, cv[i]))
      return true;
  return false;
}

/* Approximate the Hessian at each vertex, storing the result for
 * vertex `l` at `D2T[l]`. The user should have already allocated and
 * initialized `D2T`. Entries which are `NAN` which will be filled,
 * and those which are finite will be left alone and used to compute
 * other values. */
void eik3_get_D2T(eik3_s const *eik, dbl33 *D2T) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  jet31t const *jet = eik->jet;

  /* we also want to initialize D2T for points which are immediately
   * downwind of the diffracting edge */
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l) {
    par3_s par = eik3_get_par(eik, l);
    size_t la[3];
    size_t na = par3_get_active_inds(&par, la);
    if (na == 2 && mesh3_is_diff_edge(mesh, la)) {
      /* Get the update point */
      dbl3 xhat;
      mesh3_copy_vert(mesh, l, xhat);

      /* Find the point of diffraction */
      dbl3 xe = {0, 0, xhat[2]}; // manually project...

      /* Unit tangent vector for diffracting edge */
      dbl3 te = {0, 0, 1};

      /* Compute unit vector pointing from `xe` to `xhat` */
      dbl3 tf;
      dbl3_sub(xhat, xe, tf);
      dbl rho = dbl3_normalize(tf); // cylindrical radius for `xhat`

      /* Get eikonal jet at `xhat` */
      jet31t J = eik3_get_jet(eik, l);

      /* Get the ray direction */
      /* Compute unit vector orthogonal to `te` and `tf` (this vector
       * will be orthogonal to the ray direction) */
      dbl3 q1;
      dbl3_cross(te, tf, q1);
      assert(fabs(dbl3_dot(q1, J.Df)) < 1e-13);

      /* Compute unit vector orthogonal to `q1` and the ray
       * direction */
      dbl3 q2;
      dbl3_cross(J.Df, q1, q2);

      /* Compute first curvature outer product */
      dbl33 outer1;
      dbl3_outer(q1, q1, outer1);
      dbl33_dbl_div_inplace(outer1, rho);

      /* Compute second curvature outer product */
      dbl33 outer2;
      dbl3_outer(q2, q2, outer2);
      dbl33_dbl_div_inplace(outer2, J.f);

      /* Sum them up to get the Hessian */
      dbl33_add(outer1, outer2, D2T[l]);
    }
  }

  /** Propagate D2T: */

  dbl33 *D2T_cell = malloc(4*mesh3_ncells(mesh)*sizeof(dbl33));

  bool *has_init = malloc(mesh3_nverts(mesh)*sizeof(bool));
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    has_init[l] = dbl33_isfinite(D2T[l]);

  /* first, compute the Hessian at each cell vertex */
  for (size_t lc = 0, lv[4]; lc < mesh3_ncells(mesh); ++lc) {
    mesh3_cv(mesh, lc, lv);

    /* copy in initial values of D2T */
    for (size_t i = 0; i < 4; ++i)
      if (has_init[lv[i]])
        dbl33_copy(D2T[lv[i]], D2T_cell[4*lc + i]);

    /* get T and DT */
    jet31t J[4];
    for (size_t i = 0; i < 4; ++i) {
      J[i] = jet[lv[i]];
    }

    /* set up A */
    dbl4 A[3];
    for (size_t i = 0; i < 3; ++i) {
      dbl4_zero(A[i]);
      A[i][i] = 1;
      A[i][3] = -1;
    }

    /* get cell verts */
    dbl43 X;
    for (size_t i = 0; i < 4; ++i)
      mesh3_copy_vert(mesh, lv[i], X[i]);

    /* set up dX */
    dbl33 dX;
    for (size_t i = 0; i < 3; ++i)
      dbl3_sub(X[i], X[3], dX[i]);

    dbl33 dXinv, dXinvT;
    dbl33_copy(dX, dXinv);
    dbl33_invert(dXinv);
    dbl33_transposed(dXinv, dXinvT);

    /* set up bb33 */
    bb33 bb;
    bb33_init_from_jets(&bb, J, X);

    /* compute Hessian at each vertex */
    for (size_t i = 0; i < 4; ++i) {
      if (has_init[lv[i]])
        continue;

      dbl4 b;
      dbl4_e(b, i);

      /* compute Hessian in affine coordinates */
      dbl33 D2T_affine;
      for (size_t p = 0; p < 3; ++p) {
        for (size_t q = 0; q < 3; ++q) {
          dbl4 a[2];
          dbl4_copy(A[p], a[0]); // blech
          dbl4_copy(A[q], a[1]); // blech
          D2T_affine[p][q] = bb33_d2f(&bb, b, a);
        }
      }

      /* transform back to Cartesian and store with cell vertex */
      dbl33 tmp;
      dbl33_mul(dXinv, D2T_affine, tmp);
      dbl33_mul(tmp, dXinvT, D2T_cell[4*lc + i]);
    }
  }

  /* zero out D2T */
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    if (!has_init[l])
      dbl33_zero(D2T[l]);

  /* number of terms in weighted average */
  size_t *N = calloc(mesh3_nverts(mesh), sizeof(size_t));

  /* accumulate each D2T_cell entry into D2T */
  for (size_t lc = 0, cv[4]; lc < mesh3_ncells(mesh); ++lc) {
    mesh3_cv(mesh, lc, cv);

    /* skip this cell if its data is invalid */
    if (dbl33_isnan(D2T_cell[4*lc]))
      continue;

    for (size_t i = 0; i < 4; ++i) {
      if (has_init[cv[i]])
        continue;

      /* If this vertex was updated from a diff edge, don't use data
       * from a cell which is incident on a diff edge... */
      if (eik3_updated_from_diff_edge(eik, cv[i]) &&
          mesh3_cell_incident_on_diff_edge(mesh, lc))
        continue;

      /* If this vertex is incident on a diff edge, don't use data
       * from a cell which was updated from a diff edge */
      if (mesh3_vert_incident_on_diff_edge(mesh, cv[i]) &&
          any_cell_vert_updated_from_diff_edge(eik, cv))
        continue;

      dbl33_add_inplace(D2T[cv[i]], D2T_cell[4*lc + i]);
      ++N[cv[i]]; /* increment number of terms in average */
    }
  }

  /* normalize each entry by the number of incident cells */
  for (size_t lv = 0; lv < mesh3_nverts(mesh); ++lv) {
    if (has_init[lv])
      continue;
    size_t nvc = mesh3_nvc(mesh, lv);
    dbl33_dbl_div_inplace(D2T[lv], nvc);
  }

  free(N);
  free(has_init);
  free(D2T_cell);
}

void eik3_init_A_pt_src(eik3_s const *eik, dbl3 const xsrc, dbl *A) {
  assert(array_size(eik->bc_inds) == 1);

  /* Initialize with NANs */
  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l)
    A[l] = NAN;

  /* Set amplitude of nodes with point source BCs to 1/r */
  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);
    A[l] = 1/dbl3_dist(mesh3_get_vert_ptr(eik->mesh, l), xsrc);
  }

  /* We also need to go through and set the values of `A` for nodes
   * which don't have BCs but which nevertheless only have `lsrc` for
   * a parent to `1/r` */
  size_t lsrc = mesh3_get_vert_index(eik->mesh, xsrc);
  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l) {
    size_t la[3];
    size_t na = par3_get_active_inds(&eik->par[l], la);
    assert(na <= 1 || (la[0] != lsrc && la[1] != lsrc && la[2] != lsrc));
    if (na == 1 && la[0] == lsrc) {
      assert(!isfinite(A[l]));
      A[l] = 1/dbl3_dist(mesh3_get_vert_ptr(eik->mesh, l), xsrc);
    }
  }
}

void eik3_init_A_refl(eik3_s const *eik, dbl const *A_in, dbl *A) {
  /* Initialize with NANs */
  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l)
    A[l] = NAN;

  /* Copy amplitude from incident field */
  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);
    A[l] = A_in[l];
  }
}

void eik3_init_A_diff(eik3_s const *eik, dbl const *A_in, dbl *A) {
  /* Initialize with NANs */
  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l)
    A[l] = NAN;

  /* Copy amplitude from incident field */
  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);
    A[l] = A_in[l];
  }
}

void eik3_prop_A(eik3_s const *eik, dbl33 const *D2T, dbl *A) {
  mesh3_s const *mesh = eik->mesh;

  for (size_t i = 0, l; i < mesh3_nverts(mesh); ++i) {
    l = eik->accepted[i];

    if (!isnan(A[l]))
      continue;

    par3_s par = eik3_get_par(eik, l);
    if (par3_is_empty(&par))
      continue;

    dbl3 lam, abslam;
    size_t perm[3];
    dbl33_eigvals_sym(D2T[l], lam);
    dbl3_abs(lam, abslam);
    dbl3_argsort(abslam, perm);

    dbl kappa1 = lam[perm[2]], kappa2 = lam[perm[1]];

    uint3 la;
    dbl3 ba;
    size_t na = par3_get_active(&par, la, ba);
    assert(na > 0);

    dbl A_lam = 1;
    for (size_t j = 0; j < na; ++j) {
      assert(isfinite(A[la[j]]));
      A_lam *= pow(A[la[j]], ba[j]);
    }

    dbl3 xlam = {0, 0, 0};
    for (size_t j = 0; j < na; ++j) {
      dbl3 xj;
      mesh3_copy_vert(mesh, la[j], xj);
      for (size_t k = 0; k < 3; ++k)
        xlam[k] += ba[j]*xj[k];
    }

    dbl3 x;
    mesh3_copy_vert(mesh, l, x);

    dbl L = dbl3_dist(x, xlam);

    A[l] = A_lam*exp(-L*(kappa1 + kappa2)/2);
  }
}

void eik3_get_t_in(eik3_s const *eik, dbl3 *t_in) {
  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l)
    dbl3_nan(t_in[l]);

  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);
    dbl3_copy(eik->jet[l].Df, t_in[l]);
  }

  eik3_transport_unit_vector(eik, t_in, true);
}

void eik3_get_t_out(eik3_s const *eik, dbl3 *t_out) {
  mesh3_s const *mesh = eik->mesh;

  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    dbl3_nan(t_out[l]);

  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);

    /* The `t_out` vector is undefined on a BC node that lies on a
     * diffracting edge...
     *
     * TODO: usually... Will probably need to fix this later */
    if (mesh3_vert_incident_on_diff_edge(mesh, l))
      continue;

    dbl3_copy(eik->jet[l].Df, t_out[l]);

    if (!mesh3_bdv(mesh, l))
      continue;

    dbl33 R;
    mesh3_get_R_for_interior_reflector_vertex(mesh, l, R);

    dbl33_dbl3_mul_inplace(R, t_out[l]);
  }

  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    if (eik3_updated_from_diff_edge(eik, l))
      dbl3_copy(eik->jet[l].Df, t_out[l]);

  eik3_transport_unit_vector(eik, t_out, true);
}
