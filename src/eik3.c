#include "eik3.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

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
  jet3 *jet;
  state_e *state;
  int *pos;
  par3_s *par;
  heap_s *heap;

  /* In some cases, we'll skip old updates that might be useful at a
   * later stage. We keep track of them here. */
  array_s *old_updates;
  array_s *old_bd_utri; // old two-point boundary `utri`

  /* Vertices incident on faces with "reflection" boundary data. */
  bool *refl_bdv;

  /* Statistics... */
  int num_accepted; // number of nodes fixed using `eik3_step`
};

void eik3_alloc(eik3_s **eik) {
  *eik = malloc(sizeof(eik3_s));
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

  eik->jet = malloc(nverts*sizeof(jet3));
  for (size_t l = 0; l < nverts; ++l) {
    eik->jet[l] = (jet3) {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN};
  }

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

  array_alloc(&eik->old_updates);
  array_init(eik->old_updates, sizeof(utetra_s *), 16);

  array_alloc(&eik->old_bd_utri);
  array_init(eik->old_bd_utri, sizeof(utri_s *), 16);

  eik->refl_bdv = calloc(mesh3_nverts(mesh), sizeof(bool));
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
    utri_deinit(utri);
    utri_dealloc(&utri);
  }
  array_deinit(eik->old_bd_utri);
  array_dealloc(&eik->old_bd_utri);

  free(eik->refl_bdv);
  eik->refl_bdv = NULL;
}

static bool can_update_from_point(eik3_s const *eik, size_t l) {
  return eik->state[l] == VALID && !eik3_is_point_source(eik, l);
}

static void do_1pt_update(eik3_s *eik, size_t l, size_t l0) {
  // Compute new jet for one point update
  jet3 jet;
  dbl const *x = mesh3_get_vert_ptr(eik->mesh, l);
  dbl const *x0 = mesh3_get_vert_ptr(eik->mesh, l0);
  dbl3_sub(x, x0, &jet.fx);
  jet.f = dbl3_normalize(&jet.fx);

  // TODO: we should try to check if this update leaves the domain...

  // Commit the update
  //
  // TODO: we assume here that if we're actually doing a one-point
  // update "for keeps", we know that it's going to improve the
  // solution. This may not be a great approach, but it's what we're
  // doing for now...
  assert(jet.f <= eik->jet[l].f);

  eik3_set_jet(eik, l, jet);

  par3_s par = {
    .l = {l0, NO_PARENT, NO_PARENT},
    .b = {1, NAN, NAN}
  };

  eik3_set_par(eik, l, par);
}

static bool commit_tri_update(eik3_s *eik, size_t lhat, utri_s const *utri) {
  if (utri_get_value(utri) >= eik->jet[lhat].f)
    return false;

  jet3 jet;
  utri_get_jet(utri, &jet);
  eik3_set_jet(eik, lhat, jet);

  eik3_set_par(eik, lhat, utri_get_par(utri));

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
      size_t f[3] = {l, l0, l1};
      if (can_update_from_point(eik, l1) && !eik3_is_point_source(eik, l1) &&
          mesh3_bdf(eik->mesh, f)) {
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
          commit_tri_update(eik, l, utri[i])) {
        break;
      }
    } else if (i + 1 < nup &&
               utris_yield_same_update(utri[i], utri[i + 1]) &&
               utri_update_ray_is_physical(utri[i], eik) &&
               commit_tri_update(eik, l, utri[i])) {
      break;
    } else {
      array_append(eik->old_bd_utri, &utri[i]);
      free_utri[i] = false;
    }
  }

  for (size_t i = 0; i < nup; ++i) {
    if (utri[i] == NULL || !free_utri[i]) continue;
    utri_deinit(utri[i]);
    utri_dealloc(&utri[i]);
  }

  free(free_utri);
  free(utri);
  free(ve);
}

static bool can_update_from_face(eik3_s const *eik,
                                 size_t l0, size_t l1, size_t l2) {
  return can_update_from_point(eik, l0) &&
    can_update_from_point(eik, l1) &&
    can_update_from_point(eik, l2);
}

static
size_t get_update_fan(eik3_s const *eik, size_t l0, size_t **l1, size_t **l2) {
  // 1. get faces surrounding l0

  int nve = mesh3_nve(eik->mesh, l0);
  size_t (*ve)[2] = malloc(nve*sizeof(size_t[2]));
  mesh3_ve(eik->mesh, l0, ve);

  // 2. count number of cells with exactly three valid or shadow states

  int ntri = 0;
  bool *updateable_tri = malloc(nve*sizeof(bool));
  for (int i = 0; i < nve; ++i)
    ntri += updateable_tri[i] =
      can_update_from_face(eik, l0, ve[i][0], ve[i][1]);

  // 3. return if we didn't find any triangles we can update from

  if (ntri == 0)
    goto cleanup;

  // 4. allocate some space for the l1 and l2 indices that we're going
  // to pull out of these update triangles

  *l1 = malloc(ntri*sizeof(size_t));
  *l2 = malloc(ntri*sizeof(size_t));

  // 5. walk the update triangles and grab the vertices that are
  // distinct from l0, storing them in l1 and l2

  int k = 0;
  for (int i = 0; i < nve; ++i) {
    if (!updateable_tri[i])
      continue;
    (*l1)[k] = ve[i][0];
    (*l2)[k++] = ve[i][1];
  }
  assert(k == ntri);

  // 6. clean up and return

cleanup:
  free(updateable_tri);
  free(ve);

  return ntri;
}

static
void get_valid_incident_diff_edges(eik3_s const *eik, size_t l0, array_s *l1) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  int nvv = mesh3_nvv(mesh, l0);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l0, vv);

  edge_s edge = {.l = {[0] = l0}};
  for (int i = 0; i < nvv; ++i) {
    edge.l[1] = vv[i];
    if (can_update_from_point(eik, edge.l[1]) &&
        mesh3_is_diff_edge(mesh, edge.l))
      array_append(l1, &edge.l[1]);
  }

  free(vv);
}

static bool commit_tetra_update(eik3_s *eik, size_t lhat, utetra_s const *utetra) {
  if (utetra_get_value(utetra) >= eik->jet[lhat].f)
    return false;

  jet3 jet;
  utetra_get_jet(utetra, &jet);
  eik3_set_jet(eik, lhat, jet);

  eik3_set_par(eik, lhat, utetra_get_parent(utetra));

  return true;
}

static void update(eik3_s *eik, size_t l, size_t l0) {
  /**
   * The first thing we do is check if we're trying to update from a
   * point source. In this case, we just solve the two-point BVP to
   * high accuracy and return immediately.
   */
  if (eik3_is_point_source(eik, l0)) {
    do_1pt_update(eik, l, l0);
    return;
  }

  // Next, if `l` is a boundary point, we want to do any two-point
  // updates that are immersed in the boundary. (These are "creeping
  // rays", which can be physical.)
  if (mesh3_bdv(eik->mesh, l) && mesh3_bdv(eik->mesh, l0))
    do_2pt_bd_updates(eik, l, l0);

  /**
   * First, find the "update triangle fan"
   */
  size_t *l1, *l2;
  size_t num_utetra = get_update_fan(eik, l0, &l1, &l2);
  if (num_utetra == 0)
    return;

  /**
   * Before doing tetrahedron updates, we want to check if there are
   * any diffracting edges updates that aren't adjacent to `l0`. These
   * won't be covered by `do_all_diff_edge_updates_and_adjust` in
   * `update_neighbors`.
   */

  // First, check which of the l1's and l2's are adjacent to l0

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
  bool *is_diff_edge = malloc(num_utetra*sizeof(bool));
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

    // TODO: the following section cries out for refactoring...

    size_t l_active = utri_get_active_ind(u);
    assert(l_active != (size_t)NO_INDEX);

    size_t l_inactive = utri_get_inactive_ind(u);
    assert(l_inactive != (size_t)NO_INDEX);

    size_t num_inc = mesh3_get_num_inc_diff_edges(eik->mesh, l_active);
    size_t (*le)[2] = malloc(num_inc*sizeof(size_t[2]));
    mesh3_get_inc_diff_edges(eik->mesh, l_active, le);

    for (size_t k = 0; k < num_inc; ++k) {
      /* Skip the current edge */
      if ((le[k][0] == l_active ^ le[k][0] == l_inactive) &&
          (le[k][1] == l_active ^ le[k][1] == l_inactive))
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
          commit_tri_update(eik, l, u)) {
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
          commit_tri_update(eik, l, u)) {
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
    utri_deinit(u);
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
    if (utetra_has_interior_point_solution(utetra[i])) {
      if (utetra_update_ray_is_physical(utetra[i], eik) &&
          commit_tetra_update(eik, l, utetra[i]))
        break;
    } else {
      size_t num_int = utetra_get_num_interior_coefs(utetra[i]);
      assert(num_int == 1 || num_int == 2);
      size_t num_adj = 4 - num_int;
      if (i + num_adj <= num_utetra &&
          utetras_yield_same_update((utetra_s const **)&utetra[i], num_adj) &&
          utetra_update_ray_is_physical(utetra[i], eik) &&
          commit_tetra_update(eik, l, utetra[i])) {
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
  assert(!eik3_is_point_source(eik, l0));
  assert(!eik3_is_point_source(eik, l1));

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

  // TODO: don't need to allocate a bunch of separate `utri`s...
  utri_s **utri = malloc(nnb*sizeof(utri_s *));

  size_t l;

  // TODO: a problem with what I'm doing here: may need to do adjacent
  // diffracting edge updates. This could be a problem if I have
  // curved obstacle, or if a ray goes around the corner of an
  // obstacle. We'll see how far we can get with this for now...

  for (int i = 0; i < nnb; ++i) {
    array_get(nb, i, &l);

    // Do a triangle update for each neighbor l of the diffracting
    // edge (l0, l1)...
    utri_alloc(&utri[i]);
    utri_spec_s spec = utri_spec_from_eik(eik, l, l0, l1);
    utri_init(utri[i], &spec);
    utri_solve(utri[i]);

    // ... and attempt to commit it.
    if (utri_has_interior_point_solution(utri[i]) &&
        utri_update_ray_is_physical(utri[i], eik) &&
        commit_tri_update(eik, l, utri[i]))
      adjust(eik, l);
  }

  // Free triangle updates
  for (int i = 0; i < nnb; ++i) {
    if (utri[i] == NULL) continue;
    utri_deinit(utri[i]);
    utri_dealloc(&utri[i]);
  }
  free(utri);

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
  state_e l0_state = eik->state[l0];
  assert(l0_state == VALID);

  size_t l; // Node l is a neighbor of node l0

  // Get i0's neighboring nodes.
  int nnb = mesh3_nvv(eik->mesh, l0);
  size_t *nb = malloc(nnb*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, nb);

  if (l0_state == VALID) {
    // Set FAR nodes to TRIAL and insert them into the heap.
    for (int i = 0; i < nnb; ++i) {
      if (eik->state[l = nb[i]] == FAR) {
        eik->state[l] = TRIAL;
        heap_insert(eik->heap, l);
      }
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
      utri_deinit(old_utri);
      utri_dealloc(&old_utri);
      array_delete(eik->old_bd_utri, i - 1);
    }
  }
}

size_t eik3_step(eik3_s *eik) {
  size_t l0 = heap_front(eik->heap);

  assert(eik->state[l0] == TRIAL);
  assert(isfinite(eik->jet[l0].f));

  heap_pop(eik->heap);

  eik->state[l0] = VALID;

  purge_old_updates(eik, l0);
  update_neighbors(eik, l0);

  ++eik->num_accepted;

  return l0;
}

void eik3_solve(eik3_s *eik) {
  while (heap_size(eik->heap) > 0) {
    (void)eik3_step(eik);
  }
}

void eik3_add_trial(eik3_s *eik, size_t l, jet3 jet) {
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

  /* Update whether this is a vertex with a reflection BC */
  eik->refl_bdv[l] = mesh3_bdv(eik->mesh, l);
}

bool eik3_is_point_source(eik3_s const *eik, size_t l) {
  return isfinite(eik->jet[l].f) && jet3_is_nan(&eik->jet[l]);
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

jet3 eik3_get_jet(eik3_s const *eik, size_t l) {
  return eik->jet[l];
}

void eik3_set_jet(eik3_s *eik, size_t l, jet3 jet) {
  eik->jet[l] = jet;
}

jet3 *eik3_get_jet_ptr(eik3_s const *eik) {
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

void eik3_get_DT(eik3_s const *eik, size_t l, dbl DT[3]) {
  memcpy(DT, &eik->jet[l].fx, 3*sizeof(dbl));
}

bool eik3_is_refl_bdf(eik3_s const *eik, size_t const l[3]) {
  return eik->refl_bdv[l[0]] || eik->refl_bdv[l[1]] ||
    eik->refl_bdv[l[2]];
}
