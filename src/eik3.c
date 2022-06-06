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
#include "eik3_transport.h"
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
  mesh3_s const *mesh;

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

void eik3_init(eik3_s *eik, mesh3_s const *mesh) {
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

  array_alloc(&eik->bc_inds);
  array_init(eik->bc_inds, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  array_alloc(&eik->trial_inds);
  array_init(eik->trial_inds, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  eik->has_bc = calloc(nverts, sizeof(bool));

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

  array_deinit(eik->bc_inds);
  array_dealloc(&eik->bc_inds);

  array_deinit(eik->trial_inds);
  array_dealloc(&eik->trial_inds);

  free(eik->has_bc);
  eik->has_bc = NULL;

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

/* TODO: delete me!! */
static void do_utetra(eik3_s *eik, size_t lhat, size_t l[3]);
static void do_diff_utri(eik3_s *eik, size_t lhat, size_t l[2]);
static void do_bc_updates(eik3_s *eik, size_t l, size_t l0) {
  assert(eik->has_bc[l0]);
  assert(eik->state[l] == TRIAL);
  assert(eik->state[l0] == VALID);

  mesh3_s const *mesh = eik->mesh;
  array_s const *bc_inds = eik->bc_inds;

  /* Do boundary face BC updates... */

  size_t nf = mesh3_get_num_inc_bdf(mesh, l0);
  size_t (*lf)[3] = malloc(nf*sizeof(size_t[3]));
  mesh3_get_inc_bdf(mesh, l0, lf);

  for (size_t i = 0; i < nf; ++i)
    if (array_contains(bc_inds, &lf[i][0]) &&
        array_contains(bc_inds, &lf[i][1]) &&
        array_contains(bc_inds, &lf[i][2]))
      do_utetra(eik, l, lf[i]);

  /* Do diffracting edge BC updates... */

  size_t ne = mesh3_get_num_inc_diff_edges(mesh, l0);
  size_t (*le)[2] = malloc(ne*sizeof(size_t[2]));
  mesh3_get_inc_diff_edges(mesh, l0, le);

  for (size_t i = 0; i < ne; ++i)
    if (array_contains(bc_inds, &le[i][0]) &&
        array_contains(bc_inds, &le[i][1]))
      do_diff_utri(eik, l, le[i]);

  free(le);
  free(lf);
}

/** Functions for `do_bd_utri` and `do_diff_utri`: */

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

static bool did_utri_already(array_s const *utri_cache, utri_s const *utri) {
  for (size_t j = 0; j < array_size(utri_cache); ++j) {
    utri_s *utri_other;
    array_get(utri_cache, j, &utri_other);
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

  /* TODO: don't need to init utri to check this! */
  if (did_utri_already(utri_cache, utri))
    goto cleanup;

  if (utri_is_backwards(utri, eik))
    goto cleanup;

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

static void do_diff_utri(eik3_s *eik, size_t lhat, size_t l[2]) {
  do_utri(eik, lhat, l[0], l[1], eik->old_diff_utri);
}

static void do_utris_if(eik3_s *eik, size_t l, size_t l0, array_s *utri_cache,
                        bool (*pred)(eik3_s const *, size_t const[2])) {
  assert(l != l0);
  assert(!eik->has_bc[l0]);

  /* find the diffracting edges incident on l0 with VALID indices */
  array_s *l1_arr;
  array_alloc(&l1_arr);
  array_init(l1_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);
  get_valid_inc_edges(eik, l0, l1_arr, pred);

  for (size_t i = 0, l1; i < array_size(l1_arr); ++i) {
    array_get(l1_arr, i, &l1);
    if (l == l1)
      continue;

    /* If `l1` has boundary, we do the incident boundary updates
     * instead and skip this update. */
    if (eik->has_bc[l1]) {
      do_bc_updates(eik, l, l1);
      continue;
    }

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

  // TODO: not 100% sure what the right way to set this is...
  dbl const tol = 2*pow(mesh3_get_mean_edge_length(eik->mesh), 3);

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

    if (!utetras_have_same_minimizer(utetra, utetra_other, tol))
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

static bool did_utetra_already(eik3_s const *eik, utetra_s const *utetra) {
  for (size_t j = 0; j < array_size(eik->old_utetra); ++j) {
    utetra_s *utetra_other;
    array_get(eik->old_utetra, j, &utetra_other);
    if (utetras_have_same_inds(utetra, utetra_other))
      return true;
  }
  return false;
}

void do_utetra(eik3_s *eik, size_t lhat, size_t l[3]) {
  utetra_spec_s spec = utetra_spec_from_eik_and_inds(eik, lhat, l[0], l[1], l[2]);

  utetra_s *utetra;
  utetra_alloc(&utetra);
  utetra_init(utetra, &spec);

  /* TODO: we don't need to init the utetra in order to check this! */
  if (did_utetra_already(eik, utetra))
    goto cleanup;

  if (utetra_is_backwards(utetra, eik))
    goto cleanup;

  if (utetra_is_degenerate(utetra))
    goto cleanup;

  utetra_solve(utetra, /* warm start: */ NULL);

  if (utetra_get_value(utetra) >= eik->jet[lhat].f)
    goto cleanup;

  if (utetra_ray_is_occluded(utetra, eik))
    goto cleanup;

  if (utetra_has_interior_point_solution(utetra)) {
    commit_utetra(eik, lhat, utetra);
    adjust(eik, lhat);
    goto cleanup;
  }

  size_t num_interior = utetra_get_num_interior_coefs(utetra);
  assert(num_interior == 0 || num_interior == 2);

  size_t num_utetra_other;
  utetra_s **utetra_other =
    find_and_delete_cached_utetra(eik, l[0], utetra, &num_utetra_other);
  if (utetra_other) {
    commit_utetra(eik, lhat, utetra);
    adjust(eik, lhat);
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
  assert(!eik->has_bc[l0]);

  /* Get the fan of `VALID` triangles incident on `l0`. */
  array_s *le_arr;
  array_alloc(&le_arr);
  array_init(le_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);
  get_update_fan(eik, l0, le_arr);

  /* Do them all */
  for (size_t i = 0, le[2]; i < array_size(le_arr); ++i) {
    array_get(le_arr, i, &le);

    /* If the other valid node has a boundary data, then we should
     * update from the associated boundary updates and skip this
     * update. */
    if (eik->has_bc[le[0]] || eik->has_bc[le[1]]) {
      for (size_t j = 0; j < 2; ++j)
        if (eik->has_bc[le[j]])
          do_bc_updates(eik, l, le[j]);
      continue;
    }

    do_utetra(eik, l, (size_t[3]) {l0, le[0], le[1]});
  }

  /* Clean up */
  array_deinit(le_arr);
  array_dealloc(&le_arr);
}

static void update(eik3_s *eik, size_t l, size_t l0) {
  if (eik->has_bc[l0]) {
    do_bc_updates(eik, l, l0);
    return;
  }

  bool l0_is_on_diff_edge = mesh3_vert_incident_on_diff_edge(eik->mesh, l0);
  bool l_is_on_diff_edge = mesh3_vert_incident_on_diff_edge(eik->mesh, l);

  /* If `l0` is incident on a diffracting edge, look for corresponding
   * two-point updates to do. Do not do any other types of updates! */
  if (l0_is_on_diff_edge && !l_is_on_diff_edge) {
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
  assert(eik->has_bc[*l0] || !par3_is_empty(&eik->par[*l0])
         || array_contains(eik->trial_inds, l0));
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

  array_append(eik->trial_inds, &l);
}

void eik3_add_bc(eik3_s *eik, size_t l, jet31t jet) {
  assert(!array_contains(eik->bc_inds, &l));

  eik->jet[l] = jet;
  eik->state[l] = VALID;
  eik->has_bc[l] = true;
  eik->accepted[eik->num_accepted++] = l;

  array_append(eik->bc_inds, &l);
}

void do_1pt_update(eik3_s *eik, size_t l, size_t l0) {
  mesh3_s const *mesh = eik->mesh;

  dbl3 x, x0;
  mesh3_copy_vert(mesh, l, x);
  mesh3_copy_vert(mesh, l0, x0);

  jet31t jet;
  dbl3_sub(x, x0, jet.Df);
  jet.f = dbl3_normalize(jet.Df);

  if (jet.f >= eik->jet[l].f)
    return;

  eik->jet[l] = jet;

  eik->par[l].l[0] = l0;
  eik->par[l].b[0] = 1;

  adjust(eik, l);
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
    if (!eik3_has_BCs(eik, le[i]))
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

/* Set up BCs for each diffracting edge for this diffractor */
static void add_diff_bcs_for_diffractor(eik3_s *eik, eik3_s const *eik_in,
                                        size_t diff_index, array_s *le_arr) {
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

  for (size_t i = 0; i < num_diff_edges; ++i)
    if (!array_contains(le_arr, &le[i]))
      array_append(le_arr, &le[i]);

  free(le);
}

/* Add edge diffraction BCs in a tube surrounding diff_index */
void eik3_add_diff_bcs(eik3_s *eik, eik3_s const *eik_in, size_t diff_index, dbl rfac) {
  assert(eik->mesh == eik_in->mesh);
  mesh3_s const *mesh = eik3_get_mesh(eik);

  /* Array to keep track of the unique vertices on the diffractor */
  array_s *le_arr;
  array_alloc(&le_arr);
  array_init(le_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);

  /* Add BCs for the diffractor */
  add_diff_bcs_for_diffractor(eik, eik_in, diff_index, le_arr);

  /* Array containing the unique diffractor vertices */
  array_s *l_uniq_arr;
  array_alloc(&l_uniq_arr);
  array_init(l_uniq_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* ... fill it */
  for (size_t i = 0, le[2]; i < array_size(le_arr); ++i) {
    array_get(le_arr, i, &le);
    for (size_t j = 0; j < 2; ++j)
      if (!array_contains(l_uniq_arr, &le[j]))
        array_append(l_uniq_arr, &le[j]);
  }

  /* The queue driving the BFS */
  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* ... initialize it with the unique vertices on the diffractor. */
  for (size_t i = 0; i < array_size(l_uniq_arr); ++i)
    array_append(queue, array_get_ptr(l_uniq_arr, i));

  /* Main BFS loop: */
  while (!array_is_empty(queue)) {
    size_t l;
    array_pop_front(queue, &l);

    /* Set the node to trial and insert it into the heap */
    if (eik3_is_far(eik, l))
      eik3_add_trial(eik, l, jet31t_make_empty());

    /* If the current node isn't on the diffractor, then update it
     * from the diffractor
     *
     * TODO: should be optimized, otherwise it's O(N^4/3) time */
    if (!array_contains(l_uniq_arr, &l))
      for (size_t i = 0; i < array_size(le_arr); ++i)
        do_diff_utri(eik, l, array_get_ptr(le_arr, i));
    assert(isfinite(eik->jet[l].f));

    /* Check if we're in the "factoring tube"... */
    if (eik->jet[l].f - get_par_eik_value(eik, l) > rfac)
      continue;

    /* ... if we are, then add this node's neighbors to the queue */
    size_t nvv = mesh3_nvv(mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(mesh, l, vv);
    for (size_t i = 0; i < nvv; ++i)
      if (isinf(eik->jet[vv[i]].f) && !array_contains(queue, &vv[i]))
        array_append(queue, &vv[i]);
    free(vv);
  }

  array_deinit(queue);
  array_dealloc(&queue);

  array_deinit(l_uniq_arr);
  array_dealloc(&l_uniq_arr);

  array_deinit(le_arr);
  array_dealloc(&le_arr);
}

void eik3_add_refl_bcs(eik3_s *eik, eik3_s const *eik_in, size_t refl_index, dbl rfac) {
  assert(eik->mesh == eik_in->mesh);
  mesh3_s const *mesh = eik->mesh;

  /** Set up the reflection BCs: */

  size_t nf = mesh3_get_reflector_size(mesh, refl_index);
  size_t (*lf)[3] = malloc(nf*sizeof(size_t[3]));
  mesh3_get_reflector(mesh, refl_index, lf);

  /* Add reflection BCs for each vertex on the reflector */
  for (size_t i = 0; i < nf; ++i) {
    dbl33 R;
    mesh3_get_R_for_face(mesh, lf[i], R);

    for (size_t j = 0, l; j < 3; ++j) {
      l = lf[i][j];
      if (eik3_has_BCs(eik, l))
        continue;

      jet31t jet = eik3_get_jet(eik_in, l);
      dbl33_dbl3_mul_inplace(R, jet.Df);

      eik3_add_bc(eik, l, jet);
    }
  }

  array_s *l_uniq_arr;
  array_alloc(&l_uniq_arr);
  array_init(l_uniq_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < nf; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (!array_contains(l_uniq_arr, &lf[i][j]))
        array_append(l_uniq_arr, &lf[i][j]);

  /** Set up diff BCs for each incident diffractor: */

  /* Array used to keep track of all diffracting edges incident on the
   * reflector which we need to do updates from */
  array_s *le_arr;
  array_alloc(&le_arr);
  array_init(le_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);

  /* Add BCs for each diffracting edge and accumulate the unique edges
   * into `le_arr` */
  size_t num_diffractors = mesh3_get_num_diffs_inc_on_refl(mesh, refl_index);
  size_t *diff_index = malloc(num_diffractors*sizeof(size_t));
  mesh3_get_diffs_inc_on_refl(mesh, refl_index, diff_index);
  for (size_t i = 0; i < num_diffractors; ++i)
    add_diff_bcs_for_diffractor(eik, eik_in, diff_index[i], le_arr);
  free(diff_index);

  /** Using BFS, do boundary updates from the reflector's facets: */

  /* The queue driving the BFS */
  array_s *queue;
  array_alloc(&queue);
  array_init(queue, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* ... initialize it with the unique reflector vertices. */
  for (size_t i = 0; i < array_size(l_uniq_arr); ++i)
    array_append(queue, array_get_ptr(l_uniq_arr, i));

  /* Main BFS loop: */
  while (!array_is_empty(queue)) {
    size_t l;
    array_pop_front(queue, &l);

    /* Set the node to trial and insert it into the heap */
    if (eik3_is_far(eik, l))
      eik3_add_trial(eik, l, jet31t_make_empty());

    /* Do updates */
    if (!array_contains(l_uniq_arr, &l)) {
      for (size_t i = 0; i < array_size(le_arr); ++i)
        do_diff_utri(eik, l, array_get_ptr(le_arr, i));
      for (size_t i = 0; i < nf; ++i)
        do_utetra(eik, l, lf[i]);
    }
    assert(isfinite(eik->jet[l].f));

    /* Check if we're in the "factoring tube"... */
    if (eik->jet[l].f - get_par_eik_value(eik, l) > rfac)
      continue;

    /* ... if we are, then add this node's neighbors to the queue */
    size_t nvv = mesh3_nvv(mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(mesh, l, vv);
    for (size_t i = 0; i < nvv; ++i)
      if (isinf(eik->jet[vv[i]].f) && !array_contains(queue, &vv[i]))
        array_append(queue, &vv[i]);
    free(vv);
  }

  /* Clean up */

  array_deinit(queue);
  array_dealloc(&queue);

  array_deinit(le_arr);
  array_dealloc(&le_arr);

  free(lf);

  array_deinit(l_uniq_arr);
  array_dealloc(&l_uniq_arr);
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
  return eik->has_bc[l];
}

size_t const *eik3_get_accepted_ptr(eik3_s const *eik) {
  return eik->accepted;
}

size_t eik3_num_bc(eik3_s const *eik) {
  size_t num_bc = 0;
  for (size_t l = 0; l < mesh3_nverts(eik->mesh); ++l)
    num_bc += eik->has_bc[l];
  return num_bc;
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

void eik3_get_origins(eik3_s const *eik, dbl *origin) {
  mesh3_s const *mesh = eik->mesh;

  /* Initialize `origin` to `NAN` */
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    origin[l] = NAN;

  /* Set the origin of all BC nodes to 1 */
  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);
    origin[l] = 1;
  }

  /* Set the origin of all the initial `TRIAL` nodes to 1 */
  for (size_t i = 0, l; i < array_size(eik->trial_inds); ++i) {
    array_get(eik->trial_inds, i, &l);
    origin[l] = 1;
  }

  /* Set the origin of all nodes on a diffracting edge to 0 */
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    if (isnan(origin[l]) && mesh3_vert_incident_on_diff_edge(mesh, l))
      origin[l] = 0;

  /* Now, transport the origins, skipping already set values */
  eik3_transport_dbl(eik, origin, true);

  /* Finally, for each vertex on a diffracting edge, if it was updated
   * from a node with an origin equal to 1, we reset its origin value
   * to 0.5 (this corrects the level set approximating the shadow
   * boundary at these points) */
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l) {
    if (!mesh3_vert_incident_on_diff_edge(mesh, l))
      continue;

    par3_s const *par = &eik->par[l];

    size_t l_active[3];
    size_t num_active = par3_get_active_inds(par, l_active);

    for (size_t i = 0; i < num_active; ++i)
      if (origin[l_active[i]] == 1)
        origin[l] = 0.5;
  }
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

  size_t lsrc = mesh3_get_vert_index(eik->mesh, xsrc);
  A[lsrc] = INFINITY;

  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);
    A[l] = 1/dbl3_dist(mesh3_get_vert_ptr(eik->mesh, l), xsrc);
  }
}

void eik3_init_A_refl(eik3_s const *eik, dbl const *A_in, dbl *A) {
  for (size_t i = 0, l; i < array_size(eik->bc_inds); ++i) {
    array_get(eik->bc_inds, i, &l);
    A[l] = A_in[l];
  }
}

void eik3_init_A_diff(eik3_s const *eik, dbl const *A_in, dbl *A) {
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

    dbl A_lam = 1;
    assert(isfinite(par.b[0]));
    for (size_t j = 0; j < 3; ++j) {
      if (isfinite(par.b[j])) {
        assert(isfinite(A[par.l[j]]));
        A_lam *= pow(A[par.l[j]], par.b[j]);
      }
    }

    dbl3 xlam = {0, 0, 0};
    for (size_t j = 0; j < 3; ++j) {
      if (isfinite(par.b[j])) {
        dbl3 x_;
        mesh3_copy_vert(mesh, par.l[j], x_);
        for (size_t k = 0; k < 3; ++k)
          xlam[k] += par.b[j]*x_[k];
      }
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
