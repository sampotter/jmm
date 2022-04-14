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
#include "slerp.h"
#include "utetra.h"
#include "util.h"
#include "utri.h"
#include "vec.h"

typedef struct bde_bc {
  size_t le[2];
  bb31 bb;
} bde_bc_s;

bde_bc_s make_bde_bc(size_t const le[2], bb31 const *bb) {
  assert(le[0] != le[1]);
  bde_bc_s bc = {.le = {le[0], le[1]}, .bb = *bb};
  if (le[0] > le[1]) {
    SWAP(bc.le[0], bc.le[1]);
    bb31_reverse(&bc.bb);
  }
  return bc;
}

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
  jet32t *jet;
  state_e *state;
  int *pos;
  par3_s *par;
  heap_s *heap;

  /* In some cases, we'll skip old updates that might be useful at a
   * later stage. We keep track of them here. */
  array_s *old_updates;
  array_s *old_bd_utri; // old two-point boundary `utri`

  /* If `has_bc[l] == true`, then node `l` had boundary conditions
   * specified. */
  bool *has_bc;

  /* Boundary conditions for a (diffracting) boundary edge, which are
   * just cubic polynomials defined over the edge. These are used to
   * perform two-point updates from diffracting edges initially when
   * solving edge diffraction problems. */
  array_s *bde_bc;

  /* Convergence tolerances. The parameter `h` is an estimate of the
   * fineness of the mesh. The variables `h2` and `h3` are convenience
   * variables containing `h^2` and `h^3`. */
  dbl h, h2, h3;

  dbl slerp_tol;

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

  eik->jet = malloc(nverts*sizeof(jet32t));
  for (size_t l = 0; l < nverts; ++l)
    eik->jet[l] = jet32t_make_empty();

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

  array_alloc(&eik->old_updates);
  array_init(eik->old_updates, sizeof(utetra_s *), 16);

  array_alloc(&eik->old_bd_utri);
  array_init(eik->old_bd_utri, sizeof(utri_s *), 16);

  eik->has_bc = calloc(nverts, sizeof(bool));

  array_alloc(&eik->bde_bc);
  array_init(eik->bde_bc, sizeof(bde_bc_s), ARRAY_DEFAULT_CAPACITY);

  eik->h = mesh3_get_min_edge_length(mesh);
  eik->h2 = eik->h*eik->h;
  eik->h3 = eik->h*eik->h*eik->h;

  eik->slerp_tol = eik->h3;

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

  utetra_s *utetra;
  for (size_t i = 0; i < array_size(eik->old_updates); ++i) {
    array_get(eik->old_updates, i, &utetra);
    utetra_deinit(utetra);
    utetra_dealloc(&utetra);
  }
  array_deinit(eik->old_updates);
  array_dealloc(&eik->old_updates);

  utri_s *utri;
  for (size_t i = 0; i < array_size(eik->old_bd_utri); ++i) {
    array_get(eik->old_bd_utri, i, &utri);
    utri_dealloc(&utri);
  }
  array_deinit(eik->old_bd_utri);
  array_dealloc(&eik->old_bd_utri);

  free(eik->has_bc);
  eik->has_bc = NULL;

  array_deinit(eik->bde_bc);
  array_dealloc(&eik->bde_bc);

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

static bool commit_utri(eik3_s *eik, size_t lhat, utri_s const *utri) {
  if (utri_get_value(utri) >= eik->jet[lhat].f)
    return false;

  // TODO: set T
  // TODO: set DT
  // TODO: set D2T

  return true;
}

static void do_2pt_bd_updates(eik3_s *eik, size_t l, size_t l0) {
  size_t nve = mesh3_nve(eik->mesh, l);
  size_t (*ve)[2] = malloc(nve*sizeof(size_t[2]));
  mesh3_ve(eik->mesh, l, ve);

  utri_s **utri = calloc(nve, sizeof(size_t *));

  size_t l1;
  size_t nup = 0;
  for (size_t i = 0; i < nve; ++i) {
    utri[i] = NULL;
    if (ve[i][0] == l0 || ve[i][1] == l0) {
      l1 = ve[i][0] == l0 ? ve[i][1] : ve[i][0];
      size_t lf[3] = {l, l0, l1};
      if (eik->state[l1] == VALID && mesh3_is_bdf(eik->mesh, lf)) {
        utri_alloc(&utri[i]);
        utri_spec_s spec = utri_spec_from_eik(eik, l, l0, l1);
        ++nup;
        utri_init(utri[i], &spec);
        utri_solve(utri[i]);
      }
    }
  }

  qsort(utri, nve, sizeof(utri_s *), (compar_t)utri_cmp);

  for (size_t i = 0; i < nup; ++i)
    assert(utri[i] != NULL);
  for (size_t i = nup; i < nve; ++i)
    assert(utri[i] == NULL);

  /* Go through old boundary triangle updates and append any that have
   * the same update index and share an edge with the current updates
   * in `utri`. We will have solved these already, so don't need to
   * resolve them. */
  size_t copied_utri = 0;
  for (size_t i = 0; i < nve; ++i) {
    if (utri[i] == NULL)
      continue;
    utri_s *u_old;
    for (size_t j = array_size(eik->old_bd_utri); j > 0; --j) {
      array_get(eik->old_bd_utri, j - 1, &u_old);
      if (utri_get_l(u_old) != l ||
          !utri_opt_inc_on_other_utri(u_old, utri[i]))
        continue;
      array_delete(eik->old_bd_utri, j - 1);
      if (nup + copied_utri + 1 > nve)
        utri = realloc(utri, (nup + copied_utri + 1)*sizeof(utri_s *));
      utri[nup + copied_utri++] = u_old;
    }
  }
  nup += copied_utri;

  qsort(utri, nup, sizeof(utri_s *), (compar_t)utri_cmp);

  // Keep track of which utri to free. If we copy any over to
  // `old_bd_utri`, we want to make sure we don't accidentally free
  // them.
  bool *free_utri = malloc(nup*sizeof(bool));
  for (size_t i = 0; i < nup; ++i)
    free_utri[i] = true;

  // Try to commit a triangle update
  //
  // TODO: as always, this is a bit complicated. Would be nice to
  // simplify this or at least factor it out somewhere else.
  for (size_t i = 0; i < nup; ++i) {
    if (!utri_is_finite(utri[i]))
      break;
    if (utri_has_interior_point_solution(utri[i])) {
      if (utri_update_ray_is_physical(utri[i], eik) &&
          commit_utri(eik, l, utri[i])) {
        break;
      }
    } else if (i + 1 < nup &&
               utris_yield_same_update(utri[i], utri[i + 1]) &&
               utri_update_ray_is_physical(utri[i], eik) &&
               commit_utri(eik, l, utri[i])) {
      break;
    } else {
      array_append(eik->old_bd_utri, &utri[i]);
      free_utri[i] = false;
    }
  }

  for (size_t i = 0; i < nup; ++i) {
    if (utri[i] == NULL || !free_utri[i]) continue;
    utri_dealloc(&utri[i]);
  }

  free(free_utri);
  free(utri);
  free(ve);
}

static
size_t get_update_fan(eik3_s const *eik, size_t l0, size_t **l1, size_t **l2) {
  size_t nvc = mesh3_nvc(eik->mesh, l0);
  size_t *vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(eik->mesh, l0, vc);

  array_s *le_arr;
  array_alloc(&le_arr);
  array_init(le_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);

  bool is_valid[4];
  size_t cv[4], le[2], num_valid;
  for (size_t i = 0; i < nvc; ++i) {
    mesh3_cv(eik->mesh, vc[i], cv);

    num_valid = 0;
    for (size_t j = 0; j < 4; ++j)
      num_valid += is_valid[j] = eik->state[cv[j]] == VALID;

    if (num_valid != 3)
      continue;

    size_t k = 0;
    for (size_t j = 0; j < 4; ++j)
      if (cv[j] != l0 && is_valid[j])
        le[k++] = cv[j];
    assert(k == 2);

    SORT2(le[0], le[1]);

    if (eik->state[le[0]] != VALID || eik->state[le[1]] != VALID)
      continue;

    if (array_contains(le_arr, &le))
      continue;

    array_append(le_arr, &le);
  }

  size_t num_updates = array_size(le_arr);

  *l1 = malloc(num_updates*sizeof(size_t));
  *l2 = malloc(num_updates*sizeof(size_t));

  for (size_t i = 0; i < array_size(le_arr); ++i) {
    array_get(le_arr, i, &le);
    (*l1)[i] = le[0];
    (*l2)[i] = le[1];
  }

  array_deinit(le_arr);
  array_dealloc(&le_arr);

  free(vc);

  return num_updates;
}

static
void get_valid_incident_diff_edges(eik3_s const *eik, size_t l0,
                                   array_s *l1) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  int nvv = mesh3_nvv(mesh, l0);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l0, vv);

  size_t le[2] = {[0] = l0};
  for (int i = 0; i < nvv; ++i) {
    le[1] = vv[i];

    if (!mesh3_is_diff_edge(mesh, le))
      continue;

    if (!eik3_is_valid(eik, le[1]))
      continue;

    array_append(l1, &le[1]);
  }

  free(vv);
}

static void commit_utetra(eik3_s *eik, size_t lhat, utetra_s const *utetra) {
  /* If the value isn't a strict improvement, return early and do
   * nothing.
   *
   * TODO: at some point, we may want to look into dealing with the
   * case where the new value is approximately equal to the current
   * value. This is a sign that a caustic has formed. */
  if (utetra_get_value(utetra) >= eik->jet[lhat].f)
    return;

  utetra_get_jet32t(utetra, eik, &eik->jet[lhat]);
}

static void update(eik3_s *eik, size_t l, size_t l0) {
  // Next, if `l` is a boundary point, we want to do any two-point
  // updates that are immersed in the boundary. (These are "creeping
  // rays", which can be physical.)
  if (mesh3_bdv(eik->mesh, l) && mesh3_bdv(eik->mesh, l0))
    do_2pt_bd_updates(eik, l, l0);

  /**
   * First, find the "update triangle fan"
   */
  // TODO: would be better to use array_s for l1 and l2 and allocate
  // them here, so it's clear that it's this function's responsibility
  // to free them
  size_t *l1, *l2;
  size_t num_utetra = get_update_fan(eik, l0, &l1, &l2);
  if (num_utetra == 0) {
    free(l1);
    free(l2);
    return;
  }

  /**
   * Before doing tetrahedron updates, we want to check if there are
   * any diffracting edges updates that aren't adjacent to `l0`. These
   * won't be covered by `do_all_diff_edge_updates_and_adjust` in
   * `update_neighbors`.
   */

  // First, check which of the l1's and l2's are adjacent to l

  // TODO: the way we're checking for adjacent updates here is pretty
  // inefficient, but not sure if we can do better...

  bool *l_l1_adj = malloc(num_utetra*sizeof(bool));
  for (size_t i = 0; i < num_utetra; ++i)
    l_l1_adj[i] = mesh3_is_edge(eik->mesh, (size_t[2]) {l, l1[i]});

  bool *l_l2_adj = malloc(num_utetra*sizeof(bool));
  for (size_t i = 0; i < num_utetra; ++i)
    l_l2_adj[i] = mesh3_is_edge(eik->mesh, (size_t[2]) {l, l2[i]});

  // Count and mark the non-adjacent edges are diffracting edges

  size_t num_diff_edges = 0;
  bool *is_diff_edge = calloc(num_utetra, sizeof(bool));
  for (size_t i = 0; i < num_utetra; ++i) {
    if (l_l1_adj[i] || l_l2_adj[i]) {
      is_diff_edge[i] = false;
      continue;
    } else {
      size_t e[2] = {l1[i], l2[i]};
      is_diff_edge[i] = mesh3_is_diff_edge(eik->mesh, e);
      num_diff_edges += is_diff_edge[i];
    }
  }

  array_s *diff_utri;
  array_alloc(&diff_utri);
  array_init(diff_utri, sizeof(utri_s *), 2*num_diff_edges);

  utri_s *u;
  utri_spec_s spec;
  for (size_t i = 0; i < num_utetra; ++i) {
    if (!is_diff_edge[i]) continue;

    utri_alloc(&u);

    array_append(diff_utri, &u);

    spec = utri_spec_from_eik(eik, l, l1[i], l2[i]);
    spec.orig_index = i;
    utri_init(u, &spec);
    utri_solve(u);
    if (utri_has_interior_point_solution(u))
      continue;

    // TODO: <wtf?>

    size_t l_active = utri_get_active_ind(u);
    assert(l_active != (size_t)NO_INDEX);

    size_t l_inactive = utri_get_inactive_ind(u);
    assert(l_inactive != (size_t)NO_INDEX);

    size_t num_inc = mesh3_get_num_inc_diff_edges(eik->mesh, l_active);
    size_t (*le)[2] = malloc(num_inc*sizeof(size_t[2]));
    mesh3_get_inc_diff_edges(eik->mesh, l_active, le);

    for (size_t k = 0; k < num_inc; ++k) {
      /* Skip the current edge */
      if (((le[k][0] == l_active) ^ (le[k][0] == l_inactive)) &&
          ((le[k][1] == l_active) ^ (le[k][1] == l_inactive)))
        continue;

      /* Check if this edge is already in `diff_utri`. */
      bool already_found = false;
      for (size_t m = 0; m < array_size(diff_utri); ++m) {
        utri_s *u_;
        array_get(diff_utri, m, &u_);
        if (utri_contains_update_ind(u_, le[k][0]) &&
            utri_contains_update_ind(u_, le[k][1])) {
          already_found = true;
          break;
        }
      }
      if (already_found)
        continue;

      /* Add the edge! */
      utri_alloc(&u);
      array_append(diff_utri, &u);
      spec = utri_spec_from_eik(eik, l, le[k][0], le[k][1]);
      utri_init(u, &spec);
      utri_solve(u);
    }

    free(le);

    // TODO: </wtf?>
  }

  array_sort(diff_utri, (compar_t)utri_cmp);

  bool *updated_from_diff_edge = calloc(num_utetra, sizeof(bool));

  for (size_t i = 0, j; i < array_size(diff_utri); ++i) {
    utri_s *u, *u_;
    array_get(diff_utri, i, &u);
    if (!utri_is_finite(u))
      continue;
    if (utri_has_interior_point_solution(u)) {
      if (utri_update_ray_is_physical(u, eik) &&
          commit_utri(eik, l, u)) {
        if (utri_has_orig_index(u)) {
          j = utri_get_orig_index(u);
          updated_from_diff_edge[j] = true;
        }
        break;
      }
    } else if (i + 1 < array_size(diff_utri)) {
      array_get(diff_utri, i + 1, &u_);
      if (utris_yield_same_update(u, u_) &&
          utri_update_ray_is_physical(u, eik) &&
          commit_utri(eik, l, u)) {
        if (utri_has_orig_index(u)) {
          j = utri_get_orig_index(u);
          updated_from_diff_edge[j] = true;
        }
        if (utri_has_orig_index(u_)) {
          j = utri_get_orig_index(u_);
          updated_from_diff_edge[j] = true;
        }
      }
    }
  }

  for (size_t i = 0; i < array_size(diff_utri); ++i) {
    utri_s *u;
    array_get(diff_utri, i, &u);
    utri_dealloc(&u);
  }

  array_deinit(diff_utri);
  array_dealloc(&diff_utri);

  free(is_diff_edge);

  free(l_l1_adj);
  free(l_l2_adj);

  /**
   * Now we move on to doing tetrahedron updates
   */

  // Allocate tetrahedron updates
  utetra_s **utetra = malloc(num_utetra*sizeof(utetra_s *));

  // Do each tetrahedron update and sort
  for (size_t i = 0; i < num_utetra; ++i) {
    utetra_alloc(&utetra[i]);

    utetra_spec_s spec = utetra_spec_from_eik_and_inds(eik, l, l0, l1[i], l2[i]);
    if (!utetra_init(utetra[i], &spec))
      continue;

    // This is a gross hack. What we do here is prioritize a
    // diffracting edge that's incident on this tetrahedron update. It
    // might yield a somewhat higher value, but when we're close to a
    // diffracting edge, it's important to correct the ray to ensure
    // that it emits from the diffracting edge. So, we skip the
    // tetrahedron update here.
    if (updated_from_diff_edge[i])
      continue;

    // TODO: move into utetra_init?
    if (utetra_is_degenerate(utetra[i]))
      continue;

    utetra_solve(utetra[i], NULL);
  }

  free(updated_from_diff_edge);

  // Go through old updates and append any that have the same update
  // index (`l`) and share an edge with the updates currently in
  // `utetra`. Note: we will have already solved these updates! No
  // need to redo them.
  //
  // TODO: this is a bit of a mess :-(
  size_t copied_utetra = 0;
  for (size_t i = 0; i < num_utetra; ++i) {
    utetra_s *old_utetra;
    for (size_t j = array_size(eik->old_updates); j > 0; --j) {
      array_get(eik->old_updates, j - 1, &old_utetra);
      // TODO: should filter out the entries of `old_updates` that
      // have `lhat == l` beforehand so that we don't have to do this
      // check for each element of `utetra[]`... wasteful
      if (utetra_get_l(old_utetra) != l)
        continue;
      if (!utetra_opt_inc_on_other_utetra(old_utetra, utetra[i]))
        continue;
      // TODO: ACTUALLY, need to check if `old_utetra`'s optimum is
      // incident on `utetra[i]`'s base!
      array_delete(eik->old_updates, j - 1);
      utetra = realloc(
        utetra, (num_utetra + copied_utetra + 1)*sizeof(utetra_s *));
      utetra[num_utetra + copied_utetra++] = old_utetra;
    }
  }
  num_utetra += copied_utetra;

  // Sort the updates by their eikonal value
  qsort(utetra, num_utetra, sizeof(utetra_s *), (compar_t)utetra_cmp);

  // Keep track of which updates to free. If we copy any updates over
  // to `old_updates`, we want to make sure we don't accidentally free
  // them.
  bool *should_free_utetra = malloc(num_utetra*sizeof(bool));
  for (size_t i = 0; i < num_utetra; ++i)
    should_free_utetra[i] = true;

  // See if we can commit a tetrahedron update
  for (size_t i = 0; i < num_utetra; ++i) {
    if (!isfinite(utetra_get_value(utetra[i])))
      break;
    if (utetra_has_interior_point_solution(utetra[i]) &&
        utetra_update_ray_is_physical(utetra[i], eik)) {
      commit_utetra(eik, l, utetra[i]);
      break;
    } else {
      size_t num_int = utetra_get_num_interior_coefs(utetra[i]);
      assert(num_int == 1 || num_int == 2);
      size_t num_adj = 4 - num_int;
      if (i + num_adj <= num_utetra &&
          utetras_yield_same_update((utetra_s const **)&utetra[i], num_adj) &&
          utetra_update_ray_is_physical(utetra[i], eik)) {
        commit_utetra(eik, l, utetra[i]);
        break;
      } else {
        array_append(eik->old_updates, &utetra[i]);
        should_free_utetra[i] = false;
      }
    }
  }

  for (size_t i = 0; i < num_utetra; ++i)
    if (should_free_utetra[i]) {
      utetra_deinit(utetra[i]);
      utetra_dealloc(&utetra[i]);
    }
  free(utetra);

  free(should_free_utetra);

  free(l1);
  free(l2);
}

static void adjust(eik3_s *eik, size_t l) {
  assert(eik->state[l] == TRIAL);
  assert(l < mesh3_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l]);
}

size_t eik3_peek(eik3_s const *eik) {
  return heap_front(eik->heap);
}

void do_diff_edge_updates_and_adjust(eik3_s *eik, size_t l0, size_t l1,
                                     size_t *l0_nb, int l0_nnb) {
  int l1_nnb = mesh3_nvv(eik->mesh, l1);
  size_t *l1_nb = malloc(l1_nnb*sizeof(size_t));
  mesh3_vv(eik->mesh, l1, l1_nb);

  array_s *nb;
  array_alloc(&nb);
  array_init(nb, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  // TODO: should really do this work in the calling function, passing
  // in a copy of the array each time...
  for (int i = 0; i < l0_nnb; ++i)
    if (eik3_is_trial(eik, l0_nb[i]))
      array_append(nb, &l0_nb[i]);

  for (int i = 0; i < l1_nnb; ++i)
    if (eik3_is_trial(eik, l1_nb[i]) && !array_contains(nb, &l1_nb[i]))
      array_append(nb, &l1_nb[i]);

  free(l1_nb);

  int nnb = array_size(nb);

  utri_s *utri;
  utri_alloc(&utri);

  size_t l;

  // TODO: a problem with what I'm doing here: may need to do adjacent
  // diffracting edge updates. This could be a problem if I have
  // curved obstacle, or if a ray goes around the corner of an
  // obstacle. We'll see how far we can get with this for now...

  /* Do a triangle update for each neighbor of the diff. edge */
  for (int i = 0; i < nnb; ++i) {
    array_get(nb, i, &l);

    if (l == l0 || l == l1)
      continue;

    utri_spec_s spec = utri_spec_from_eik(eik, l, l0, l1);

    utri_init(utri, &spec);

    if (utri_is_degenerate(utri))
      continue;

    utri_solve(utri);

    if ((utri_has_interior_point_solution(utri) ||
         utri_emits_terminal_ray(utri, eik)) &&
        utri_update_ray_is_physical(utri, eik) &&
        commit_utri(eik, l, utri))
      adjust(eik, l);
  }

  utri_dealloc(&utri);

  array_deinit(nb);
  array_dealloc(&nb);
}

void do_all_diff_edge_updates_and_adjust(eik3_s *eik, size_t l0, size_t *l0_nb,
                                         int l0_nnb) {
  array_s *l1;
  array_alloc(&l1);
  array_init(l1, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  get_valid_incident_diff_edges(eik, l0, l1);

  for (size_t i = 0; i < array_size(l1); ++i)
    do_diff_edge_updates_and_adjust(
      eik, l0, *(size_t *)array_get_ptr(l1, i), l0_nb, l0_nnb);

  array_deinit(l1);
  array_dealloc(&l1);
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

  // Find newly VALID diffracting edges (l0, l1) and update all
  // neighboring nodes---these are nodes neighboring *both* l0 and
  // l1. If we don't look at all of these neighbors, we could miss
  // important edge updates.
  do_all_diff_edge_updates_and_adjust(eik, l0, nb, nnb);

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
static void purge_old_updates(eik3_s *eik, size_t l0) {
  utetra_s *old_utetra;
  for (size_t i = array_size(eik->old_updates); i > 0; --i) {
    array_get(eik->old_updates, i - 1, &old_utetra);
    if (utetra_get_l(old_utetra) == l0) {
      utetra_deinit(old_utetra);
      utetra_dealloc(&old_utetra);
      array_delete(eik->old_updates, i - 1);
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
  assert(eik->state[*l0] == TRIAL);
  heap_pop(eik->heap);
  eik->state[*l0] = VALID;

  /* If the eikonal of the newly VALID node isn't finite, something
   * bad happened. */
  if (!isfinite(eik->jet[*l0].f))
    return JMM_ERROR_RUNTIME_ERROR;

  purge_old_updates(eik, *l0);
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

void eik3_add_trial(eik3_s *eik, size_t l, jet32t jet) {
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

jet32t eik3_get_jet32t(eik3_s const *eik, size_t l) {
  return eik->jet[l];
}

void eik3_set_jet(eik3_s *eik, size_t l, jet32t jet) {
  eik->jet[l] = jet;
}

jet32t *eik3_get_jet_ptr(eik3_s const *eik) {
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

dbl eik3_get_slerp_tol(eik3_s const *eik) {
  return eik->slerp_tol;
}

bool eik3_has_BCs(eik3_s const *eik, size_t l) {
  return eik->has_bc[l];
}

dbl eik3_get_h(eik3_s const *eik) {
  return eik->h;
}

bool eik3_get_refl_bdf_inc_on_diff_edge(eik3_s const *eik, size_t const le[2],
                                        size_t lf[3]) {
  size_t nbdf = mesh3_get_num_bdf_inc_on_edge(eik->mesh, le);
  size_t (*bdf)[3] = malloc(nbdf*sizeof(size_t[3]));
  mesh3_get_bdf_inc_on_edge(eik->mesh, le, bdf);

  size_t i;

  for (i = 0; i < nbdf; ++i) {
    if (eik3_has_BCs(eik, bdf[i][0]) &&
        eik3_has_BCs(eik, bdf[i][1]) &&
        eik3_has_BCs(eik, bdf[i][2])) {
      memcpy(lf, bdf[i], sizeof(size_t[3]));
      break;
    }
  }

  bool found = i < nbdf;

#if JMM_DEBUG
  for (i = i + 1; i < nbdf; ++i) {
    if (eik3_has_BCs(eik, bdf[i][0]) &&
        eik3_has_BCs(eik, bdf[i][1]) &&
        eik3_has_BCs(eik, bdf[i][2])) {
      assert(false);
    }
  }
#endif

  free(bdf);

  return found;
}

size_t const *eik3_get_accepted_ptr(eik3_s const *eik) {
  return eik->accepted;
}
