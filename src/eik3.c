#include "eik3.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include "array.h"
#include "bb.h"
#include "edgemap.h"
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

/**
 * A structure managing a jet marching method solving the eikonal
 * equation in 3D on an unstructured tetrahedron mesh.
 *
 * NOTE: this is just for s = 1 at the moment. Will extend this to
 * handle s != later.
 */
struct eik3 {
  mesh3_s *mesh;
  jet3 *jet;
  state_e *state;
  int *pos;
  par3_s *par;
  heap_s *heap;
  edgemap_s *cutset;

  /* A margin for the cutedge parameter `t`. When being computed, and
   * before clamping to [0, 1], we assert that t falls in the range
   * [-`tlim`, 1 + `tlim`], where `tlim` is set to the square of
   * `mesh3_get_min_tetra_alt(eik->mesh)`. */
  dbl tlim;

  /**
   * In some cases, we'll skip old updates that might be useful at a
   * later stage. We keep track of them here.
   */
  array_s *old_updates;

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

  edgemap_alloc(&eik->cutset);
  edgemap_init(eik->cutset, sizeof(cutedge_s));

  eik->tlim = pow(mesh3_get_min_tetra_alt(eik->mesh), 2);

  array_alloc(&eik->old_updates);
  array_init(eik->old_updates, sizeof(utetra_s *), 16);
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

  edgemap_deinit(eik->cutset);
  edgemap_dealloc(&eik->cutset);

  utetra_s *utetra;
  for (size_t i = 0; i < array_size(eik->old_updates); ++i) {
    array_get(eik->old_updates, i, &utetra);
    utetra_deinit(utetra);
    utetra_dealloc(&utetra);
  }
  array_deinit(eik->old_updates);
  array_dealloc(&eik->old_updates);
}

static cutedge_s get_cutedge(eik3_s const *eik, size_t l0, size_t l1) {
  cutedge_s cutedge;
  bool found = edgemap_get(eik->cutset, make_edge(l0, l1), &cutedge);
  assert(found);

  // Reverse `cutedge.t` if we need to. Remember that we want (1 -
  // t)*x0 + t*x1 to equal to point where the shadow boundary
  // intersects each edge in the shadow cut.
  if (l0 > l1)
    cutedge.t = 1 - cutedge.t;

  return cutedge;
}

static bool is_shadow_p1_diff(eik3_s const *eik, size_t l0, size_t l) {
  dbl const atol = 1e-14;
  dbl n[3];
  dbl3_cross(&eik->jet[l0].fx, &eik->jet[l].fx, n);
  return dbl3_norm(n) > atol;
}

/**
 * Determine if a node with two parents that was updated from a
 * diffracting edge (l0, l1) is a shadow node.
 *
 * To do this, we compute the cross product between the gradients at
 * nodes l0 and l1 to approximate dZ with a plane. Then, we check
 * whether l0 lies on the far side of that or not.
 *
 * TODO:
 * - [ ] A better way to do this would be to use the normal vectors
 *   (that have the correct orientation!) to check which side l0 is
 *   on.
 */
static bool is_shadow_p2_diff(eik3_s const *eik, size_t l0, size_t const *l) {
  dbl const atol = 1e-14;
  dbl n[3];
  dbl3_cross(&eik->jet[l[0]].fx, &eik->jet[l[1]].fx, n);
  dbl dot = dbl3_dot(&eik->jet[l0].fx, n);
  return fabs(dot) > atol;
}

static bool is_shadow_p2(eik3_s const *eik, par3_s const *parent) {
  size_t const *l = parent->l;

  if (eik->state[l[0]] == eik->state[l[1]])
    return eik->state[l[0]] == SHADOW;

  dbl t = parent->b[1], t_shadow = get_cutedge(eik, l[0], l[1]).t;

  if (eik->state[l[0]] == VALID)
    return t >= t_shadow;
  else
    return t <= t_shadow;
}

static bool is_shadow_p3(eik3_s const *eik, par3_s par) {
  dbl const atol = 1e-14;

  assert(par3_size(&par) == 3);

  size_t num_shadow = 0;
  for (int i = 0; i < 3; ++i)
    num_shadow += eik3_is_shadow(eik, par.l[i]);

  assert(num_shadow <= 3);

  // Handle the easy cases first (if all the nodes are `SHADOW` or
  // `VALID`).
  if (num_shadow == 0 || num_shadow == 3)
    return num_shadow == 3;

  // We use a different "check state" depending on whether two of the
  // nodes are shadow, or only one is. In either case, we move the
  // node with the "check state" to `par.l[0]`. This lets us easily
  // use the standard triangle to check whether the optimum is
  // emanating from the shadow zone or not.
  state_e check_state = num_shadow == 1 ? SHADOW : VALID;

  // First, swap the entries of `par` so that the first entry
  // corresponds to the "check state" node.
  if (eik->state[par.l[1]] == check_state) {
    SWAP(par.b[1], par.b[0]);
    SWAP(par.l[1], par.l[0]);
  }
  if (eik->state[par.l[2]] == check_state) {
    SWAP(par.b[2], par.b[0]);
    SWAP(par.l[2], par.l[0]);
  }

  // Get coefficients of the optimum in the standard triangle after
  // permuting.
  dbl t[2] = {par.b[1], par.b[2]};

  // Get the coordinates of the intersection of dZ with the two shadow
  // cutset edges.
  dbl ts[2] = {
    get_cutedge(eik, par.l[0], par.l[1]).t,
    get_cutedge(eik, par.l[0], par.l[2]).t
  };

  // TODO: deal with this case here to avoid NaNs below
  if (ts[0] == 0 && ts[1] == 0) {
    if (check_state == SHADOW)
      return t[0] == 0 && t[1] == 0;
    else
      return t[0] > 0 || t[1] > 0;
  }

  // Return `true` if the optimum is on the side of the line
  // connecting the intersection points on the shadow cutset nearer to
  // x0.
  if (check_state == SHADOW)
    return t[0] <= ts[0] + atol && t[1] <= ts[1]*(1 - t[0]/ts[0]) + atol;
  else
    return t[0] >= ts[0] - atol || t[1] >= ts[1]*(1 - t[0]/ts[0]) - atol;
}

static bool is_shadow(eik3_s const *eik, size_t l0) {
  if (eik3_is_point_source(eik, l0))
    return false;

  mesh3_s const *mesh = eik3_get_mesh(eik);

  par3_s const *parent = &eik->par[l0];

  int num_parents = par3_size(parent);

  size_t const *l = parent->l;

  // TODO: eventually, we should order these cases in decreasing order
  // of frequency to optimize this

  if (num_parents == 1 && eik3_is_point_source(eik, l[0]))
    return eik3_is_shadow(eik, l[0]);

  if (num_parents == 1 && mesh3_vert_incident_on_diff_edge(mesh, l[0]))
    return is_shadow_p1_diff(eik, l0, l[0]);

  if (num_parents == 1)
    return eik3_is_shadow(eik, l[0]);

  if (num_parents == 2 && mesh3_is_diff_edge(mesh, l))
    return is_shadow_p2_diff(eik, l0, l);

  if (num_parents == 2)
    return is_shadow_p2(eik, parent);

  if (num_parents == 3)
    return is_shadow_p3(eik, *parent);

  die();
}

static bool can_update_from_point(eik3_s const *eik, size_t l) {
  return (eik->state[l] == VALID || eik->state[l] == SHADOW) &&
    !eik3_is_point_source(eik, l);
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
  int nve = mesh3_nve(eik->mesh, l);
  size_t (*ve)[2] = malloc(nve*sizeof(size_t[2]));
  mesh3_ve(eik->mesh, l, ve);

  utri_s **utri = malloc(nve*sizeof(utri_s *));

  size_t l1;
  int nup = 0;
  for (int i = 0; i < nve; ++i) {
    utri[i] = NULL;
    if (ve[i][0] == l0 || ve[i][1] == l0) {
      l1 = ve[i][0] == l0 ? ve[i][1] : ve[i][0];
      size_t f[3] = {l, l0, l1};
      if (can_update_from_point(eik, l1) && !eik3_is_point_source(eik, l1) &&
          mesh3_bdf(eik->mesh, f)) {
        utri_alloc(&utri[i]);
        utri_spec_s spec = utri_spec_from_eik(eik, l, l0, l1);
        if (utri_init(utri[i], &spec)) {
          utri_solve(utri[i]);
          ++nup;
        }
      }
    }
  }
  qsort(utri, nve, sizeof(utri_s *), (compar_t)utri_cmp);

  // An edge on the boundary is incident on two boundary faces
  // (assuming the surface is manifold, which we do). The other
  // vertices of these faces are the `l1`s we find above. So, we can't
  // have done more than two updates.
  assert(nup <= 2);

  // Try to commit a triangle update
  //
  // TODO: as always, this is a bit complicated. Would be nice to
  // simplify this or at least factor it out somewhere else.
  for (int i = 0; i < nup; ++i) {
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
    }
  }

  for (int i = 0; i < nve; ++i) {
    if (utri[i] == NULL) continue;
    utri_deinit(utri[i]);
    utri_dealloc(&utri[i]);
  }

  free(utri);
  free(ve);
}

static bool can_update_from_edge(eik3_s const *eik, size_t l[2]) {
  return can_update_from_point(eik, l[0]) && can_update_from_point(eik, l[1]);
}

static int get_update_tri_fan(eik3_s const *eik, size_t l0, size_t **l1, size_t **l2) {
  // 1. get faces surrounding l0

  int nve = mesh3_nve(eik->mesh, l0);
  size_t (*ve)[2] = malloc(nve*sizeof(size_t[2]));
  mesh3_ve(eik->mesh, l0, ve);

  // 2. count number of cells with exactly three valid or shadow states

  int ntri = 0;
  bool *updateable_tri = malloc(nve*sizeof(bool));
  for (int i = 0; i < nve; ++i)
    ntri += updateable_tri[i] = can_update_from_edge(eik, ve[i]);

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
  size_t num_utetra = get_update_tri_fan(eik, l0, &l1, &l2);
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
    if (!utri_init(u, &spec))
      continue;
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
      if (!utri_init(u, &spec))
        continue;
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
    if (utetra_has_interior_point_solution(utetra[i]) ||
        utetra_has_shadow_boundary_solution(utetra[i])) {
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

static bool edges_overlap(edge_s edge, void const *elt, edge_s const *other) {
  (void)elt;
  size_t const *l = other->l;
  return edge.l[0] == l[0] || edge.l[1] == l[0] ||
    edge.l[0] == l[1] || edge.l[1] == l[1];
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
    if (!utri_init(utri[i], &spec))
      continue;

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
  assert(l0_state == VALID || l0_state == SHADOW);

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

/**
 * Compute the surface normal at a diffraction edge when `l0` has two
 * parents. The index `l0` indicates a new `SHADOW` point, which we
 * use to orient the surface normal, computed as the cross product of
 * the eikonal gradients at indices `l[0]` and `l[1]`. The result goes
 * in `n`.
 */
static void get_diff_edge_surf_normal_p2(eik3_s const *eik, size_t l0,
                                         size_t l[2], dbl n[3]) {
  assert(eik3_is_shadow(eik, l0));

  // Get DT at each parent and compute their cross product.
  dbl DT[2][3];
  eik3_get_DT(eik, l[0], DT[0]);
  eik3_get_DT(eik, l[1], DT[1]);
  dbl3_cross(DT[0], DT[1], n);
  dbl3_normalize(n);

  mesh3_s const *mesh = eik3_get_mesh(eik);

  // Reorient the surface normal by checking which side of the tangent
  // plane at (x[l[0]] + x[l[1]])/2 the point x[l0] is on.
  dbl dx[3], lm[3];
  dbl3_add(mesh3_get_vert_ptr(mesh, l[0]), mesh3_get_vert_ptr(mesh, l[1]), lm);
  dbl3_dbl_div_inplace(lm, 2);
  dbl3_sub(mesh3_get_vert_ptr(mesh, l0), lm, dx);
  if (dbl3_dot(n, dx) > 0)
    dbl3_negate(n);
}

static void set_cutedge_jet_p2(eik3_s const *eik, edge_s edge,
                               size_t const l[2], cutedge_s *cutedge) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  dbl xt[3];
  edge_get_xt(edge, mesh, cutedge->t, xt);

  utri_s *utri;
  utri_alloc(&utri);
  utri_spec_s spec = utri_spec_from_eik_without_l(eik, xt, l[0], l[1]);
  if (!utri_init(utri, &spec)) {
    utri_deinit(utri);
    utri_dealloc(&utri);
    return;
  }

  utri_solve(utri);
  // TODO: check whether ray is physical? ugh

  if (utri_has_interior_point_solution(utri)) {
    par3_s utri_par = utri_get_par(utri);
    assert(!is_shadow_p2(eik, &utri_par));

    utri_get_jet(utri, &cutedge->jet);
    assert(jet3_is_finite(&cutedge->jet));

    utri_deinit(utri);
    utri_dealloc(&utri);

    return;
  }

  /* If we haven't found an interior point solution, we assume that
   * we're on a diffracting edge, and that we need to explore the
   * surrounding edges to find the correct update.
   *
   * First, we find the active endpoint, then find all of the
   * diffracting edges adjacent to it. */

  assert(mesh3_is_diff_edge(mesh, l));

  size_t l_active = utri_get_active_ind(utri);

  utri_deinit(utri);
  utri_dealloc(&utri);

  size_t num_inc = mesh3_get_num_inc_diff_edges(mesh, l_active);
  assert(num_inc > 1);

  size_t (*le)[2] = malloc(num_inc*sizeof(size_t[2]));
  mesh3_get_inc_diff_edges(mesh, l_active, le);

  utri_s **u = malloc(num_inc*sizeof(utri_s *));

  for (size_t i = 0; i < num_inc; ++i) {
    utri_alloc(&u[0]);
    utri_spec_s spec_ = utri_spec_from_eik_without_l(eik, xt, le[i][0], le[i][1]);
    if (utri_init(u[0], &spec_))
      utri_solve(u[0]);
  }

  qsort(u, num_inc, sizeof(utri_s *), (compar_t)utri_cmp);

  for (size_t i = 0; i < num_inc; ++i) {
    if (!utri_is_finite(u[i]))
      assert(false); // should never happen...
    if (utri_has_interior_point_solution(u[i]) ||
        (i + 1 < num_inc &&
         utris_yield_same_update(u[i], u[i + 1]))) {
      // TODO: check whether ray is physical? ugh
      utri_get_jet(u[i], &cutedge->jet);
      assert(jet3_is_finite(&cutedge->jet));
      break;
    }
  }

  for (size_t i = 0; i < num_inc; ++i) {
    if (u[i] == NULL) continue;
    utri_deinit(u[i]);
    utri_dealloc(&u[i]);
  }

  free(u);

  free(le);
}

static void set_cutedge_jet_p3_a1(eik3_s const *eik, dbl const xt[3],
                                  size_t l0, cutedge_s *cutedge) {
  assert(l0 != (size_t)NO_INDEX);

  size_t *l1, *l2;
  size_t num_utetra = get_update_tri_fan(eik, l0, &l1, &l2);
  assert(num_utetra != 0);

  utetra_s **u = malloc(num_utetra*sizeof(utetra_s *));

  utetra_spec_s spec;

  for (size_t i = 0; i < num_utetra; ++i) {
    utetra_alloc(&u[i]);

    spec = utetra_spec_from_eik_without_l(eik, xt, l0, l1[i], l2[i]);
    if (!utetra_init(u[i], &spec))
      continue;

    if (utetra_is_degenerate(u[i]))
      continue;

    utetra_solve(u[i], NULL);
  }

  qsort(u, num_utetra, sizeof(utetra_s *), (compar_t)utetra_cmp);

  for (size_t i = 0; i < num_utetra; ++i) {
    assert(isfinite(utetra_get_value(u[i])));
    if (utetra_has_interior_point_solution(u[i]) ||
        utetra_has_shadow_boundary_solution(u[i])) {
      // TODO: check whether ray is physical? ugh
      utetra_get_jet(u[i], &cutedge->jet);
      assert(jet3_is_finite(&cutedge->jet));
      break;
    } else {
      size_t num_int = utetra_get_num_interior_coefs(u[i]);
      assert(num_int == 1 || num_int == 2);
      size_t num_adj = 4 - num_int;
      if (i + num_adj <= num_utetra &&
          utetras_yield_same_update((utetra_s const **)&u[i], num_adj)) {
        // TODO: check whether ray is physical? ugh
        utetra_get_jet(u[i], &cutedge->jet);
        assert(jet3_is_finite(&cutedge->jet));
        break;
      }
    }
  }

  for (size_t i = 0; i < num_utetra; ++i) {
    utetra_deinit(u[i]);
    utetra_dealloc(&u[i]);
  }

  free(u);

  free(l1);
  free(l2);
}

static void set_cutedge_jet_p3_a2(eik3_s const *eik, dbl const xt[3],
                                  size_t la[2], cutedge_s *cutedge) {
  assert(la[0] != (size_t)NO_INDEX && la[1] != (size_t)NO_INDEX);

  size_t nev = mesh3_nev(eik->mesh, la);
  size_t *ev = malloc(nev*sizeof(size_t));
  mesh3_ev(eik->mesh, la, ev);

  utetra_s *u;
  utetra_alloc(&u);

  utetra_spec_s spec;

  // TODO: switch this over to use the same approach as everywhere
  // else...

  for (size_t i = 0; i < nev; ++i) {
    spec = utetra_spec_from_eik_without_l(eik, xt, la[0], la[1], ev[i]);

    if (!utetra_init(u, &spec))
      continue;

    utetra_solve(u, NULL);

    if (utetra_has_interior_point_solution(u) ||
        utetra_has_shadow_boundary_solution(u)) {
      // TODO: check whether ray is physical? ugh

      utetra_get_jet(u, &cutedge->jet);
      assert(jet3_is_finite(&cutedge->jet));

      goto cleanup;
    }
  }

cleanup:
  utetra_dealloc(&u);
  free(ev);
}

static void set_cutedge_jet_p3(eik3_s const *eik, edge_s edge,
                               size_t const l[3], cutedge_s *cutedge) {
  /* Here's a little research project for this function, which is
   * actually pretty important... If we find an update with a "shadow
   * boundary solution" (i.e. lying on the shadow boundary with
   * nonzero Lagrange multipliers for the active constraints), can we
   * move the cutset point a bit until these constraints are relaxed?
   * Is that a good idea? Seems complicated but could be
   * interesting... */

  dbl xt[3];
  edge_get_xt(edge, eik->mesh, cutedge->t, xt);

  utetra_spec_s spec = utetra_spec_from_eik_without_l(eik, xt, l[0], l[1], l[2]);

  utetra_s *utetra;
  utetra_alloc(&utetra);
  utetra_init(utetra, &spec);
  utetra_solve(utetra, NULL);

  if (utetra_has_interior_point_solution(utetra) ||
      utetra_has_shadow_boundary_solution(utetra)) {
    // TODO: check whether ray is physical? ugh

    utetra_get_jet(utetra, &cutedge->jet);
    assert(jet3_is_finite(&cutedge->jet));

    goto coda;
  }

  size_t la[3];
  switch(utetra_get_active_inds(utetra, la)) {
  case 1:
    set_cutedge_jet_p3_a1(eik, xt, la[0], cutedge);
    break;
  case 2:
    set_cutedge_jet_p3_a2(eik, xt, la, cutedge);
    break;
  default:
    assert(false);
  }

coda:
  utetra_deinit(utetra);
  utetra_dealloc(&utetra);
}

static void set_cutedge_jet(eik3_s const *eik, edge_s edge, cutedge_s *cutedge) {
  dbl t = cutedge->t;

  /* First, check if either of the parents are diffracting edges. */
  par3_s par[2] = {eik3_get_par(eik, edge.l[0]), eik3_get_par(eik, edge.l[1])};
  size_t par_size[2] = {par3_size(&par[0]), par3_size(&par[1])};

  /* Handle a few different cases depending on whether the edge
   * endpoints were updated from diffracting edges or not... */
  par3_s sel;
  if (par_size[0] == 2 && mesh3_is_diff_edge(eik->mesh, par[0].l)) {
    sel = par[0];
  } else if (par_size[1] == 2 && mesh3_is_diff_edge(eik->mesh, par[1].l)) {
    sel = par[1];
  } else {
    /* Get the edge endpoint closer to x(t) and grab that node's parents. */
    sel = t < 0.5 ? par[0] : par[1];
  }

  /* For two parents, compute the cutedge jet using a triangle update;
   * for three, use a tetrahedron update. We don't current handle a
   * single parent. */
  switch (par3_size(&sel)) {
  case 2: return set_cutedge_jet_p2(eik, edge, sel.l, cutedge);
  case 3: return set_cutedge_jet_p3(eik, edge, sel.l, cutedge);
  default: die();
  }
}

static void estimate_cutedge_data_from_incident_cutedges(
  eik3_s const *eik, edge_s edge, cutedge_s *cutedge, edgemap_s const *inc)
{
  assert(!edgemap_is_empty(inc));

  mesh3_s const *mesh = eik->mesh;

  dbl x0[3], dx[3], yt[3];

  edge_get_x0_and_dx(edge, mesh, x0, dx);

  cutedge->t = 0;
  dbl3_zero(cutedge->n);

  edgemap_iter_s *iter;
  edgemap_iter_alloc(&iter);
  edgemap_iter_init(iter, inc);

  edge_s inc_edge;
  cutedge_s inc_cutedge;

  while (edgemap_iter_next(iter, &inc_edge, &inc_cutedge)) {
    // Get the location of the intersection between the shadow
    // boundary and the current cut edge

    edge_get_xt(inc_edge, mesh, inc_cutedge.t, yt);

    dbl t = dbl3_dot(inc_cutedge.n, yt);
    t -= dbl3_dot(inc_cutedge.n, x0);
    if (t != 0) {
      dbl denom = dbl3_dot(inc_cutedge.n, dx);
      assert(denom != 0);
      t /= denom;
    }

    cutedge->t += t;

    dbl3_add_inplace(cutedge->n, inc_cutedge.n);
  }

  edgemap_iter_dealloc(&iter);

  size_t num_inc = edgemap_size(inc);

  cutedge->t /= num_inc;

  // Verify that t is (up to machine precision) in the interval [0,
  // 1]. If it's slightly outside of the interval due to roundoff
  // error, clamp it back before returning.
  assert(-eik->tlim <= cutedge->t && cutedge->t <= 1 + eik->tlim);
  cutedge->t = fmax(0, fmin(1, cutedge->t));

  // TODO: this is the perfect place to do a weighted spherical
  // average instead of just normalizing!
  dbl3_normalize(cutedge->n);

  /* Finally, set the cutedge jet. */
  set_cutedge_jet(eik, edge, cutedge);
}

typedef enum cutedge_status {
  CUTEDGE_CONTINUE,
  CUTEDGE_VALID,
  CUTEDGE_REINSERT,
  CUTEDGE_SKIP
} cutedge_status_e;

/* Orient the normal `n` so that it supports the boundary faces
 * incident on `de`. We assume that `n` is already perpendicular to
 * `de`, but this function ensures that the boundary is on the right
 * side of it. The edge indexed by `de` is assumed to be a diffracting
 * edge. */
static void orient_normal(mesh3_s const *mesh, size_t de[2], dbl n[3]) {
  assert(mesh3_is_diff_edge(mesh, de));

  dbl xm[3], dx[3];
  mesh3_get_edge_centroid(mesh, de, xm);

  size_t nev = mesh3_nev(mesh, de);
  size_t *lv = malloc(nev*sizeof(size_t));
  mesh3_ev(mesh, de, lv);

  /* Go through the vertices of the boundary faces that aren't
   * incident on the diffracting edge and, when we find one that's
   * oriented inconsistently, re-orient. */
  // TODO: if we really wanted to get crazy, seems like we could just
  // check a single, arbitrary boundary face, re-orient, and call it
  // good. In debug mode, though, we don't want to do this. Better to
  // be doubly sure there aren't any weird cases we forgot about.
  for (size_t i = 0; i < nev; ++i) {
    if (!mesh3_bdf(mesh, (size_t[3]) {de[0], de[1], lv[i]}))
      continue;
    dbl3_sub(mesh3_get_vert_ptr(mesh, lv[i]), xm, dx);
    if (dbl3_dot(n, dx) > 0) {
      dbl3_negate(n);
      break;
    }
  }

  /* Verify the normal supports the boundary faces correctly. */
  for (size_t i = 0; i < nev; ++i) {
    if (!mesh3_bdf(mesh, (size_t[3]) {de[0], de[1], lv[i]}))
      continue;
    dbl3_sub(mesh3_get_vert_ptr(mesh, lv[i]), xm, dx);
    assert(dbl3_dot(n, dx) <= 0);
  }

  free(lv);
}

static
bool set_cutedge_for_diff_vert(eik3_s const *eik, edge_s edge, cutedge_s *cutedge) {
  dbl const atol = 1e-13;

  mesh3_s const *mesh = eik->mesh;

  size_t lv = edge_get_valid_index(edge, eik);

  assert(mesh3_vert_incident_on_diff_edge(mesh, lv));

  size_t nde = mesh3_get_num_inc_diff_edges(eik->mesh, lv); // TODO: inefficient, lift!
  size_t (*de)[2] = malloc(nde*sizeof(size_t[2]));
  mesh3_get_inc_diff_edges(eik->mesh, lv, de);

  assert(nde <= 2); // TODO: handle corners later...

  dbl (*N)[3] = malloc(nde*sizeof(dbl[3]));

  dbl DT1[3];
  eik3_get_DT(eik, lv, DT1);

  dbl const *xv = mesh3_get_vert_ptr(mesh, lv);
  dbl dx[3];

  /* Start by compute the normalized cross product between the ray
   * vector at `xv` and each diffracting edge incident on `xv`. */
  for (size_t i = 0, j; i < nde; ++i) {
    assert(de[i][0] == lv || de[i][1] == lv);
    j = de[i][0] == lv ? 1 : 0;
    dbl3_sub(mesh3_get_vert_ptr(mesh, de[i][j]), xv, dx);
    dbl3_cross(DT1, dx, N[i]);
    dbl3_normalize(N[i]);
    orient_normal(mesh, de[i], N[i]);
  }

  free(de);

  // TODO: for now, ensure that all these normals are the same at this
  // point. This is because we assume that the edges are all entirely
  // straight. Later, when we relax this assumption, we'll need to
  // handle this.
  for (size_t i = 1; i < nde; ++i)
    assert(dbl3_dist(N[i], N[i - 1]) < atol);

  /* Set the surface normal (arbitrarily) to be the first of these
   * surface normals, since we're assuming they're all the same at
   * this point, anyway. */
  dbl3_copy(N[0], cutedge->n);

  free(N);

  cutedge->t = lv == edge.l[0] ? 0 : 1;
  cutedge->jet = eik->jet[lv];

  return true;
}

/* In this function, we want to fill in the cutedge information for a
 * node that was updated from a diffracting edge. This is only the
 * case when l0 has two parents. This is the second "base case": check
 * and see if l0's parents are incident on a diffracting edge.
 *
 * What we're checking here:
 * - l1 is one of l0's parents
 * - l0 has exactly two parents (why not just one == l1?)
 * - l2 is the other parent
 * - [l1, l2] is a diffracting edge
 */
static cutedge_status_e
set_cutedge_for_diff_update(eik3_s const *eik, edge_s edge, cutedge_s *cutedge) {
  size_t lv = edge_get_valid_index(edge, eik);
  size_t ls = edge_get_shadow_index(edge, eik);

  par3_s par = eik3_get_par(eik, ls);

  if (par3_size(&par) != 2)
    return CUTEDGE_CONTINUE;

  size_t *l = par.l;

  if (lv != l[0] && lv != l[1])
    return CUTEDGE_CONTINUE;

  if (!mesh3_is_diff_edge(eik->mesh, l))
    return CUTEDGE_CONTINUE;

  // TODO: I think we can simplify this a bit...
  if (lv == l[0] || lv == l[1]) {
    if (lv == l[0])
      SWAP(l[0], l[1]);
  } else if (mesh3_is_diff_edge(eik->mesh, (size_t[2]) {lv, l[1]})) {
    SWAP(l[0], l[1]);
  } else if (!mesh3_is_diff_edge(eik->mesh, (size_t[2]) {lv, l[0]})) {
    log_warn("failed to insert cutedge (%lu, %lu)\n", MIN(lv, ls), MAX(lv, ls));
    return CUTEDGE_SKIP;
  }

  assert(ls != lv && lv != l[0]);

  get_diff_edge_surf_normal_p2(eik, ls, (size_t[2]) {lv, l[0]}, cutedge->n);

  cutedge->t = lv == edge.l[0] ? 0 : 1; // This assumes that l1 is valid

  cutedge->jet = eik->jet[lv];

  return CUTEDGE_VALID;
}

/**
 * Compute the coefficient for the new edge in shadow cutset. This is
 * a double t such that 0 <= t <= 1 and where the shadow boundary
 * (approximately) passes through (1 - t)*x[l0] + t*x[l1].
 *
 * The index l0 corresponds to the node that has just been accepted,
 * and l1 is some neighbor of l0 which has the "opposite" state from
 * l0 (i.e., VALID if l0 is SHADOW and vice versa). So, one of
 * eik->state[l0] and eik->state[l1] is VALID and the other is
 * SHADOW. No assumption is made about which is which.
 */
static
cutedge_status_e set_cutedge(eik3_s const *eik, edge_s edge, cutedge_s *cutedge) {
  cutedge_status_e status;

  status = set_cutedge_for_diff_update(eik, edge, cutedge);
  if (status != CUTEDGE_CONTINUE)
    return status;

  size_t lv = edge_get_valid_index(edge, eik);

  if (mesh3_vert_incident_on_diff_edge(eik->mesh, lv) &&
      set_cutedge_for_diff_vert(eik, edge, cutedge))
    return CUTEDGE_VALID;

  /* "Inductive step": find all other cutset edges incident on `edge`
   * use them to extrapolate the shadow boundary in order to find the
   * intersection point. */

  // TODO: not sure if this is the *best* way to find nearby cutedges,
  // but seems reasonable...

  edgemap_s *inc;
  edgemap_alloc(&inc);
  edgemap_init(inc, sizeof(cutedge_s));
  edgemap_filter(eik->cutset, inc, (edgemap_prop_t)edges_overlap, &edge);

  status = edgemap_is_empty(inc) ? CUTEDGE_REINSERT : CUTEDGE_VALID;

  if (status == CUTEDGE_VALID)
    estimate_cutedge_data_from_incident_cutedges(
      eik, edge, cutedge, inc);

  edgemap_deinit(inc);
  edgemap_dealloc(&inc);

  return status;
}

static void insert_cutedge(eik3_s *eik, edge_s edge, cutedge_s *cutedge) {
  assert(0 <= cutedge->t && cutedge->t <= 1);
  assert(fabs(1 - dbl3_norm(cutedge->n)) < 1e-15);

  cutedge_s old_cutedge;

  /* If the edge is already contained in the cutset, check to see if
   * the cut point indicated by `cutedge` has a lower eikonal value
   * than the cut point already stored in the cutset. If it does,
   * replace it. (What's our rationale for this? We need a way to
   * break ties, and this seems like a reasonable way to break the
   * tie. We might find that there's a more principled way of doing
   * this later. */
  if (edgemap_get(eik->cutset, edge, &old_cutedge)) {
    // TODO: just using linear interpolation here, which is probably
    // OK. Would be good to replace this with cubic interpolation.
    dbl T0 = eik->jet[edge.l[0]].f, T1 = eik->jet[edge.l[1]].f;
    dbl t = cutedge->t, old_t = old_cutedge.t;
    dbl T_old = (1 - old_t)*T0 + old_t*T1, T = (1 - t)*T0 + t*T1;
    if (T < T_old)
      edgemap_set(eik->cutset, edge, cutedge);
  } else {
    edgemap_set(eik->cutset, edge, cutedge);
  }
}

static void update_shadow_cutset(eik3_s *eik, size_t l0) {
  size_t l1;
  edge_s edge;
  cutedge_s cutedge;

  // Cutset edges consist of two nodes with opposite states
  state_e op_state = eik->state[l0] == VALID ? SHADOW : VALID;

  // Get l0's neighbors
  int nvv = mesh3_nvv(eik->mesh, l0);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, vv);

  array_s *l1_arr;
  array_alloc(&l1_arr);
  array_init(l1_arr, sizeof(size_t), nvv);

  // Fill `l1_arr` with neighbors of the opposite state
  for (int i = 0; i < nvv; ++i)
    if (eik->state[l1 = vv[i]] == op_state)
      array_append(l1_arr, &l1);

  // TODO: document me!
  if (eik3_is_valid(eik, l0) &&
      mesh3_vert_incident_on_diff_edge(eik->mesh, l0))
    for (size_t i = 0; i < array_size(l1_arr); ++i) {
      array_get(l1_arr, i, &l1);
      edge = make_edge(l0, l1);
      set_cutedge_for_diff_vert(eik, edge, &cutedge);
      insert_cutedge(eik, edge, &cutedge);
    }

  // We process the nodes neighboring l0 that have the opposite state
  // using a queue. When we try to compute the data for each new
  // cutedge, because of the way we do it right now, we need to
  // compute this data for each node in a certain order. If we fail
  // the first time, we reinsert the node into the queue and try
  // again.
  while (!array_is_empty(l1_arr)) {
    array_pop_front(l1_arr, &l1);
    edge_s edge = make_edge(l0, l1);
    switch (set_cutedge(eik, edge, &cutedge)) {
    case CUTEDGE_CONTINUE:
      die();
    case CUTEDGE_VALID:
      break;
    case CUTEDGE_REINSERT:
      array_append(l1_arr, &l1);
    case CUTEDGE_SKIP:
      continue;
    }
    insert_cutedge(eik, edge, &cutedge);
  }

  array_deinit(l1_arr);
  array_dealloc(&l1_arr);

  free(vv);
}

static void update_statistics(eik3_s *eik) {
  ++eik->num_accepted;
}

size_t eik3_step(eik3_s *eik) {
  size_t l0 = heap_front(eik->heap);
  assert(eik->state[l0] == TRIAL);
  heap_pop(eik->heap);

  /* When we accept a node with an infinite jet, we can conclude
   * immediately that it must be a `SHADOW` node, since this should
   * only happen when a node wasn't reachable by any update. */
  if (isinf(eik->jet[l0].f)) {
    eik->state[l0] = SHADOW;
    update_shadow_cutset(eik, l0);
    update_neighbors(eik, l0);
    return l0;
  }

  /* Once we accept a node, we can clear out the old updates that
   * are targeting it */
  utetra_s *old_utetra;
  for (size_t i = array_size(eik->old_updates); i > 0; --i) {
    array_get(eik->old_updates, i - 1, &old_utetra);
    if (utetra_get_l(old_utetra) == l0) {
      utetra_deinit(old_utetra);
      utetra_dealloc(&old_utetra);
      array_delete(eik->old_updates, i - 1);
    }
  }

  eik->state[l0] = is_shadow(eik, l0) ? SHADOW : VALID;

  /* Reset a `SHADOW` node's jet upon accepting it. We never want to
   * do an update using information from a `SHADOW` node, since it's
   * garbage! Instead, we avoid this by using split updates (see
   * `utetra` and `utri`).  */
  if (eik->state[l0] == SHADOW)
    eik->jet[l0] = jet3_make_empty();

  update_shadow_cutset(eik, l0);
  update_neighbors(eik, l0);
  update_statistics(eik);

  return l0;
}

void eik3_solve(eik3_s *eik) {
  while (heap_size(eik->heap) > 0) {
    (void)eik3_step(eik);
  }
}

void eik3_add_trial(eik3_s *eik, size_t l, jet3 jet) {
  if (eik->state[l] == VALID) {
    return;
  } else if (eik->pos[l] == NO_INDEX) {
    assert(eik->state[l] == FAR);
    eik->jet[l] = jet;
    eik->state[l] = TRIAL;
    heap_insert(eik->heap, l);
  } else if (jet.f < eik->jet[l].f) {
    assert(eik->state[l] == TRIAL);
    eik->jet[l] = jet;
    adjust(eik, l);
  }
}

void eik3_add_valid(eik3_s *eik, size_t l, jet3 jet) {
  // TODO: need to decide what to do if the state is TRIAL
  if (eik->state[l] == TRIAL) {
    abort();
  }

  eik->jet[l] = jet;
  eik->state[l] = VALID;
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

bool eik3_is_shadow(eik3_s const *eik, size_t l) {
  return eik->state[l] == SHADOW;
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

void eik3_get_DT(eik3_s const *eik, size_t l, dbl DT[3]) {
  memcpy(DT, &eik->jet[l].fx, 3*sizeof(dbl));
}

edgemap_s const *eik3_get_cutset(eik3_s const *eik) {
  return eik->cutset;
}

bool eik3_get_cutedge_t(eik3_s const *eik, size_t l0, size_t l1, dbl *t) {
  cutedge_s cutedge;
  if (!edgemap_get(eik->cutset, make_edge(l0, l1), &cutedge))
    return false;
  *t = cutedge.t;
  if (l0 > l1)
    *t = 1 - *t;
  return true;
}

bool eik3_get_cutedge_jet(eik3_s const *eik, size_t l0, size_t l1, jet3 *jet) {
  cutedge_s cutedge;
  if (!edgemap_get(eik->cutset, make_edge(l0, l1), &cutedge))
    return false;
  *jet = cutedge.jet;
  return true;
}
