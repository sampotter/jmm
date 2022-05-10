#include "eik3.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include "alist.h"
#include "array.h"
#include "bb.h"
#include "edge.h"
#include "heap.h"
#include "log.h"
#include "macros.h"
#include "mat.h"
#include "mesh3.h"
#include "utetra.h"
#include "util.h"
#include "utri.h"
#include "vec.h"

// TODO:
//
// - There's a common pattern that gets reused here a lot:
//
//   1. Use a mesh iterator to enumerate the bases of a set of
//      `utetra` or `utri`
//   2. Do all these updates
//   3. Sort them
//   4. Commit an update based on what's going on with the first
//      couple of updates
//
//   I rewrote this from scratch each time I did it because I wasn't
//   sure how often I'd be doing it! But there you are. Ideally, we
//   should refactor these "updates sequences" into new module (e.g.,
//   utetras.h, utris.h), and call them from here. It will make the
//   code in this file much easier to understand and likely suss out a
//   few bugs.

/* A structure managing a jet marching method solving the eikonal
 * equation in 3D on an unstructured tetrahedron mesh.
 *
 * NOTE: this is just for s = 1 at the moment. Will extend this to
 * handle s != later. */
struct eik3 {
  mesh3_s *mesh;
  jet31t *jet;
  state_e *state;
  int *pos;
  par3_s *par;
  heap_s *heap;

  /* In some cases, we'll skip old updates that might be useful at a
   * later stage. We keep track of them here. */
  array_s *old_utetra;
  array_s *old_bd_utri; // old two-point boundary `utri`
  array_s *old_diff_utri; // old two-point updates from diff. edges

  /* If `has_bc[l] == true`, then node `l` had boundary conditions
   * specified. */
  bool *has_bc;

  /* Useful statistics for debugging */
  size_t num_accepted; /* number of nodes fixed by `eik3_step` */

  /* An array containing the order in which the individual nodes were
   * accepted. That is, `accepted[i] == l` means that `eik3_step()`
   * returned `l` when it was called for the `i`th time. */
  size_t *accepted;

  bool is_initialized;
};

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

void eik3_init(eik3_s *eik, mesh3_s *mesh) {
  eik->mesh = mesh;

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

  array_alloc(&eik->old_bd_utri);
  array_init(eik->old_bd_utri, sizeof(utri_s *), 16);

  array_alloc(&eik->old_diff_utri);
  array_init(eik->old_diff_utri, sizeof(utri_s *), 16);

  eik->has_bc = calloc(nverts, sizeof(bool));

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
    utri_s *utri;

    for (size_t i = 0; i < array_size(eik->old_utetra); ++i) {
      array_get(eik->old_utetra, i, &utetra);
      utetra_dealloc(&utetra);
    }
    array_deinit(eik->old_utetra);
    array_dealloc(&eik->old_utetra);

    for (size_t i = 0; i < array_size(eik->old_bd_utri); ++i) {
      array_get(eik->old_bd_utri, i, &utri);
      utri_dealloc(&utri);
    }
    array_deinit(eik->old_bd_utri);
    array_dealloc(&eik->old_bd_utri);

    for (size_t i = 0; i < array_size(eik->old_diff_utri); ++i) {
      array_get(eik->old_diff_utri, i, &utri);
      utri_dealloc(&utri);
    }
    array_deinit(eik->old_diff_utri);
    array_dealloc(&eik->old_diff_utri);
  }

  free(eik->has_bc);
  eik->has_bc = NULL;

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

void eik3_dump_has_bc(eik3_s const *eik, char const *path) {
  FILE *fp = fopen(path, "wb");
  fwrite(eik->has_bc, sizeof(eik->has_bc[0]), mesh3_nverts(eik->mesh), fp);
  fclose(fp);
}

static void adjust(eik3_s *eik, size_t l) {
  assert(eik->state[l] == TRIAL);
  assert(l < mesh3_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l]);
}

static void commit_utri(eik3_s *eik, size_t lhat, utri_s const *utri) {
  /* TODO: see comment about caustics in `commit_utetra` */
  assert(utri_get_value(utri) < eik->jet[lhat].f);

  utri_get_jet31t(utri, &eik->jet[lhat]);

  eik3_set_par(eik, lhat, utri_get_par(utri));
}

/** Functions for `do_bd_utri` and `do_diff_utri`: */

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

/* Look through the cache of old edge-diffracted two-point updates for
 * an update which matches `utri` (is incident on it and which shares
 * the same active index). If we find one, delete it from the cache
 * and return it. */
static utri_s *
find_and_delete_cached_utri(array_s *utri_cache, utri_s const *utri) {
  /* get indices of current `utri` */
  size_t l = utri_get_l(utri);
  size_t l_active = utri_get_active_ind(utri);
  size_t l_inactive = utri_get_inactive_ind(utri);

  utri_s *utri_other = NULL;

  /* iterate over the other `utri` in the cache... */
  for (size_t i = 0; i < array_size(utri_cache); ++i) {
    array_get(utri_cache, i, &utri_other);

    /* if this is a distinct `utri`, with the same target index and
     * the same active index (so, the inactive index must be
     * different!) ... */
    if (l == utri_get_l(utri_other) &&
        l_active == utri_get_active_ind(utri_other) &&
        l_inactive != utri_get_inactive_ind(utri_other)) {
      /* ... then delete it and break from the loop */
      array_delete(utri_cache, i);
      break;
    }

    /* make sure to set `utri_other` back to `NULL` here---if we don't
     * end up breaking, on the last loop, we will return `NULL` to
     * signal that we didn't find a matching `utri` to delete */
    utri_other = NULL;
  }

  return utri_other;
}

/* Check whether `utri` has been stored in the cache for
 * edge-diffracted updates already. */
static bool utri_is_cached(array_s const *utri_cache, utri_s const *utri) {
  utri_s const *utri_other = NULL;
  for (size_t i = 0; i < array_size(utri_cache); ++i) {
    array_get(utri_cache, i, &utri_other);
    if (utris_have_same_inds(utri, utri_other))
      return true;
  }
  return false;
}

static void
do_utri(eik3_s *eik, size_t l, size_t l0, size_t l1, array_s *utri_cache) {
  utri_spec_s spec = utri_spec_from_eik(eik, l, l0, l1);

  utri_s *utri;
  utri_alloc(&utri);
  utri_init(utri, &spec);

  for (size_t j = 0; j < array_size(utri_cache); ++j) {
    utri_s *utri_other;
    array_get(utri_cache, j, &utri_other);
    if (utris_have_same_inds(utri, utri_other))
      goto cleanup;
  }

  if (utri_is_degenerate(utri))
    goto cleanup;

  utri_solve(utri);

  if (utri_get_value(utri) >= eik->jet[l].f)
    goto cleanup;

  if (utri_ray_is_occluded(utri, eik))
    goto cleanup;

  if (utri_has_interior_point_solution(utri)) {
    commit_utri(eik, l, utri);
    adjust(eik, l);
    goto cleanup;
  }

  /* see if we can find an old triangle update which matches the
   * current one, and commit the update if we can
   *
   * (if successful, delete the old one!) */
  utri_s *utri_other = find_and_delete_cached_utri(utri_cache, utri);
  if (utri_other) {
    commit_utri(eik, l, utri);
    adjust(eik, l);
    utri_dealloc(&utri_other);
    goto cleanup;
  }

  /* if we failed, we cache this update for later (... if we
   * haven't already) */
  if (!utri_is_cached(utri_cache, utri)) {
    array_append(utri_cache, &utri);
    return; /* don't dealloc in this case! */
  }

cleanup:
  utri_dealloc(&utri);
}

void eik3_do_diff_utri(eik3_s *eik, size_t l, size_t l0, size_t l1) {
  do_utri(eik, l, l0, l1, eik->old_diff_utri);
}

static void do_utris_if(eik3_s *eik, size_t l, size_t l0, array_s *utri_cache,
                        bool (*pred)(eik3_s const *, size_t const[2])) {
  assert(l != l0);

  /* find the diffracting edges incident on l0 with VALID indices */
  array_s *l1_arr;
  array_alloc(&l1_arr);
  array_init(l1_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);
  get_valid_inc_edges(eik, l0, l1_arr, pred);

  for (size_t i = 0, l1; i < array_size(l1_arr); ++i) {
    array_get(l1_arr, i, &l1);
    if (l == l1) continue;
    do_utri(eik, l, l0, l1, utri_cache);
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
utetras_bracket_ray(utetra_s const *utetra, size_t l0,
                    utetra_s **utetra_other, size_t num_utetra_other) {
  if (num_utetra_other == 0)
    return false;

  size_t num_interior = utetra_get_num_interior_coefs(utetra);
  assert(num_interior == 0 || num_interior == 2);

  bool l0_is_active = utetra_index_is_active(utetra, l0);

  if (num_interior == 2)
    return l0_is_active;

  if (!l0_is_active)
    return false;

  /* Get the first utetra... if any of the vertices are boundary
   * vertices, we should do something a little different here. Namely,
   * we should make sure that l_start is on the boundary! But I'm not
   * sure whether this is actually that important... see comment
   * below. */

  size_t i_current = 0;

  size_t l_other[2] = {NO_INDEX, NO_INDEX};
  utetra_get_other_inds(utetra_other[i_current], l0, l_other);

  size_t l_start = l_other[0], l_current = l_other[1];

  /* Next, we walk around the ring of utetra until we've visited all
   * of them... */
  size_t num_visited = 1;
  while (num_visited++ < num_utetra_other) {
    /* Search through the utetra for the one we should step to next */
    for (size_t i = 0; i < num_utetra_other; ++i) {
      if (i == i_current)
        continue;

      utetra_get_other_inds(utetra_other[i], l0, l_other);
      if (l_current != l_other[0] && l_current != l_other[1])
        continue;

      /* Update l_current and i_current and break */
      l_current = l_current == l_other[0] ? l_other[1] : l_other[0];
      i_current = i;
      break;
    }
  }

  /* We're done walking... check if we're back where we started */
  return l_start == l_current;
}


static utetra_s **
find_and_delete_cached_utetra(eik3_s *eik, size_t l0,
                              utetra_s const *utetra, size_t *num_utetra) {
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

    if (!utetras_have_same_minimizer(utetra, utetra_other))
      continue;

    /* Initially just append here so we can check whether we actually
     * want to return these utetras! */
    array_append(i_arr, &i);
  }

  /* Get the pointers to the utetra we found */
  *num_utetra = array_size(i_arr);
  utetra_s **utetra_found = malloc(*num_utetra*sizeof(utetra_s *));
  for (size_t j = 0, i; j < array_size(i_arr); ++j) {
    array_get(i_arr, j, &i);
    array_get(eik->old_utetra, i, &utetra_found[j]);
    assert(!utetras_have_same_inds(utetra, utetra_found[j]));
  }

  /* We're good if we've found at least one utetra, and either: 1)
   * we're dealing with an optimum on the interior of an edge, or 2)
   * we have an optimum on a corner and we can "walk" the ring of
   * surrounding utetra */
  bool found = *num_utetra > 0 &&
    utetras_bracket_ray(utetra, l0, utetra_found, *num_utetra);
  if (found) {
    array_delete_all(eik->old_utetra, i_arr);
  } else {
    free(utetra_found);
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

void eik3_do_utetra(eik3_s *eik, size_t l, size_t l0, size_t l1, size_t l2) {
  utetra_spec_s spec = utetra_spec_from_eik_and_inds(eik, l, l0, l1, l2);

  utetra_s *utetra;
  utetra_alloc(&utetra);
  utetra_init(utetra, &spec);

  if (utetra_is_backwards(utetra, eik))
    goto cleanup;

  for (size_t j = 0; j < array_size(eik->old_utetra); ++j) {
    utetra_s *utetra_other;
    array_get(eik->old_utetra, j, &utetra_other);
    if (utetras_have_same_inds(utetra, utetra_other))
      goto cleanup;
  }

  if (utetra_is_degenerate(utetra))
    goto cleanup;

  utetra_solve(utetra, /* warm start: */ NULL);

  if (utetra_get_value(utetra) >= eik->jet[l].f)
    goto cleanup;

  if (utetra_ray_is_occluded(utetra, eik))
    goto cleanup;

  if (utetra_has_interior_point_solution(utetra)) {
    commit_utetra(eik, l, utetra);
    adjust(eik, l);
    goto cleanup;
  }

  size_t num_interior = utetra_get_num_interior_coefs(utetra);
  assert(num_interior == 0 || num_interior == 2);

  size_t num_utetra_other;
  utetra_s **utetra_other =
    find_and_delete_cached_utetra(eik, l0, utetra, &num_utetra_other);
  if (utetra_other) {
    commit_utetra(eik, l, utetra);
    adjust(eik, l);
    for (size_t j = 0; j < num_utetra_other; ++j)
      utetra_dealloc(&utetra_other[j]);
    free(utetra_other);
    goto cleanup;
  }

  if (!utetra_cached_already(eik, utetra)) {
    array_append(eik->old_utetra, &utetra);
    return; /* don't dealloc! */
  }

cleanup:
  utetra_dealloc(&utetra);
}

static void do_utetra_fan(eik3_s *eik, size_t l, size_t l0) {
  /* Get the fan of `VALID` triangles incident on `l0`. */
  array_s *le_arr;
  array_alloc(&le_arr);
  array_init(le_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);
  get_update_fan(eik, l0, le_arr);

  /* Do them all */
  for (size_t i = 0, le[2]; i < array_size(le_arr); ++i) {
    array_get(le_arr, i, &le);
    eik3_do_utetra(eik, l, l0, le[0], le[1]);
  }

  /* Clean up */
  array_deinit(le_arr);
  array_dealloc(&le_arr);
}

static void update(eik3_s *eik, size_t l, size_t l0) {
  bool l0_is_on_diff_edge = mesh3_vert_incident_on_diff_edge(eik->mesh, l0);

  /* If `l0` is incident on a diffracting edge, look for corresponding
   * two-point updates to do. Do not do any other types of updates! */
  if (l0_is_on_diff_edge) {
    do_utris_if(eik, l, l0, eik->old_diff_utri, is_diff_edge);
    return;
  }

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

size_t eik3_peek(eik3_s const *eik) {
  return heap_front(eik->heap);
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

/* Remove old tetrahedron and two-point boundary updates targeting
 * `l0` from `eik`. */
static void purge_old_utetra(eik3_s *eik, size_t l0) {
  utetra_s *old_utetra;
  for (size_t i = array_size(eik->old_utetra); i > 0; --i) {
    array_get(eik->old_utetra, i - 1, &old_utetra);
    if (utetra_get_l(old_utetra) == l0) {
      utetra_dealloc(&old_utetra);
      array_delete(eik->old_utetra, i - 1);
    }
  }

  utri_s *old_utri;
  for (size_t i = array_size(eik->old_bd_utri); i > 0; --i) {
    array_get(eik->old_bd_utri, i - 1, &old_utri);
    if (utri_get_l(old_utri) == l0) {
      utri_dealloc(&old_utri);
      array_delete(eik->old_bd_utri, i - 1);
    }
  }
}

jmm_error_e eik3_step(eik3_s *eik, size_t *l0) {
  /* Pop the first node from the front of the heap and set it to
   * VALID. Do some basic sanity checks. */
  *l0 = heap_front(eik->heap);
  assert(isfinite(eik->jet[*l0].f));
  assert(eik->state[*l0] == TRIAL);
  heap_pop(eik->heap);
  eik->state[*l0] = VALID;

  /* If the eikonal of the newly VALID node isn't finite, something
   * bad happened. */
  if (!isfinite(eik->jet[*l0].f))
    return JMM_ERROR_RUNTIME_ERROR;

  purge_old_utetra(eik, *l0);
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

bool eik3_is_solved(eik3_s const *eik) {
  return eik->num_accepted == mesh3_nverts(eik->mesh);
}

void eik3_add_trial(eik3_s *eik, size_t l, jet31t jet) {
  if (eik->has_bc[l]) {
    log_warn("tried to add BCs for node %lu more than once\n", l);
    return;
  }

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

  eik->has_bc[l] = true;
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

mesh3_s *eik3_get_mesh(eik3_s const *eik) {
  return eik->mesh;
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
  return eik->has_bc[l];
}

size_t const *eik3_get_accepted_ptr(eik3_s const *eik) {
  return eik->accepted;
}
