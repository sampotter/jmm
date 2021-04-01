#include "eik3.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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
  int num_valid;
  edgemap_s *cutset;

  /**
   * In some cases, we'll skip old updates that might be useful at a
   * later stage. We keep track of them here.
   */
  array_s *old_updates;
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

  eik->num_valid = 0;

  edgemap_alloc(&eik->cutset);
  edgemap_init(eik->cutset, sizeof(cutedge_s));

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

  dbl ts = get_cutedge(eik, l[0], l[1]).t;
  dbl t = parent->b[1];
  return eik->state[l[0]] == VALID ? t > ts : t <= ts;
}

static bool is_shadow_p3(eik3_s const *eik, par3_s par, int num_shadow) {
  assert(par3_size(&par) == 3);
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
    return t[0] <= ts[0] && t[1] <= ts[1]*(1 - t[0]/ts[0]);
  else
    return t[0] > ts[0] || t[1] > ts[1]*(1 - t[0]/ts[0]);
}

static bool is_shadow(eik3_s const *eik, size_t l0) {
  if (eik3_is_point_source(eik, l0))
    return false;

  mesh3_s const *mesh = eik3_get_mesh(eik);

  par3_s const *parent = &eik->par[l0];

  int num_parents = par3_size(parent);

  size_t const *l = parent->l;

  int num_valid = 0;
  for (int i = 0; i < num_parents; ++i)
    num_valid += eik3_is_valid(eik, l[i]);

  int num_shadow = 0;
  for (int i = 0; i < num_parents; ++i)
    num_shadow += eik3_is_shadow(eik, l[i]);

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
    return is_shadow_p3(eik, *parent, num_shadow);

  assert(false);
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

  size_t l[3] = {[2] = NO_PARENT};
  utri_get_update_inds(utri, l);

  dbl b[3] = {[2] = NAN};
  utri_get_bary_coords(utri, b);

  eik3_set_par(eik, lhat, make_par3(l, b));

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
    utri_alloc(&utri[i]);
    utri_reset(utri[i]);
    if (ve[i][0] == l0 || ve[i][1] == l0) {
      l1 = ve[i][0] == l0 ? ve[i][1] : ve[i][0];
      size_t f[3] = {l, l0, l1};
      if (can_update_from_point(eik, l1) && !eik3_is_point_source(eik, l1) &&
          mesh3_bdf(eik->mesh, f)) {
        utri_init_from_eik3(utri[i], eik, l, l0, l1);
        utri_solve(utri[i]);
        ++nup;
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
    if (!isfinite(utri_get_value(utri[i])))
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

  for (int i = 0; i < nve; ++i)
    utri_dealloc(&utri[i]);

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

  size_t l[3];
  utetra_get_update_inds(utetra, l);

  dbl b[3];
  utetra_get_bary_coords(utetra, b);

  eik3_set_par(eik, lhat, make_par3(l, b));

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
  int num_utetra = get_update_tri_fan(eik, l0, &l1, &l2);
  if (num_utetra == 0)
    return;

  /**
   * Before doing tetrahedron updates, we want to check if there are
   * any diffracting edges updates that aren't adjacent to `l0`. These
   * won't be covered by `do_all_diff_edge_updates_and_adjust` in
   * `update_neighbors`.
   */

  bool *l_l1_adj = malloc(num_utetra*sizeof(bool));
  for (int i = 0; i < num_utetra; ++i)
    l_l1_adj[i] = mesh3_is_edge(eik->mesh, (size_t[2]) {l, l1[i]});

  bool *l_l2_adj = malloc(num_utetra*sizeof(bool));
  for (int i = 0; i < num_utetra; ++i)
    l_l2_adj[i] = mesh3_is_edge(eik->mesh, (size_t[2]) {l, l2[i]});

  bool *is_diff_edge = malloc(num_utetra*sizeof(bool));

  int num_diff_utri = 0;
  for (int i = 0; i < num_utetra; ++i) {
    if (l_l1_adj[i] || l_l2_adj[i]) {
      is_diff_edge[i] = false;
      continue;
    } else {
      size_t e[2] = {l1[i], l2[i]};
      is_diff_edge[i] = mesh3_is_diff_edge(eik->mesh, e);
      num_diff_utri += is_diff_edge[i];
    }
  }

  utri_s *utri;
  utri_alloc(&utri);

  bool *updated_from_diff_edge = malloc(num_utetra*sizeof(bool));

  for (int i = 0, j = 0; i < num_utetra; ++i) {
    updated_from_diff_edge[i] = false;

    if (!is_diff_edge[i])
      continue;

    utri_reset(utri);
    utri_init_from_eik3(utri, eik, l, l1[j], l2[j]);
    utri_solve(utri);

    if (utri_has_interior_point_solution(utri) &&
        utri_update_ray_is_physical(utri, eik)) {
      commit_tri_update(eik, l, utri);
      updated_from_diff_edge[i] = true;
    }
  }

  utri_dealloc(&utri);

  free(is_diff_edge);
  free(l_l1_adj);
  free(l_l2_adj);

  /**
   * Now we move on to doing tetrahedron updates
   */

  // Allocate tetrahedron updates
  utetra_s **utetra = malloc(num_utetra*sizeof(utetra_s *));

  // Do each tetrahedron update and sort
  for (int i = 0; i < num_utetra; ++i) {
    utetra_alloc(&utetra[i]);
    utetra_reset(utetra[i]);

    // This is a gross hack. What we do here is prioritize a
    // diffracting edge that's incident on this tetrahedron update. It
    // might yield a somewhat higher value, but when we're close to a
    // diffracting edge, it's important to correct the ray to ensure
    // that it emits from the diffracting edge. So, we skip the
    // tetrahedron update here.
    if (updated_from_diff_edge[i])
      continue;

    if (utetra_init_from_eik3(utetra[i], eik, l, l0, l1[i], l2[i]) &&
        !utetra_is_degenerate(utetra[i])) {
      utetra_set_lambda(utetra[i], (dbl[2]) {0, 0});
      utetra_solve(utetra[i]);
    }
  }

  free(updated_from_diff_edge);

  // Go through old updates and append any that have the same update
  // index (`l`) and share an edge with the updates currently in
  // `utetra`. Note: we will have already solved these updates! No
  // need to redo them.
  //
  // TODO: this is a bit of a mess :-(
  size_t copied_utetra = 0;
  for (int i = 0; i < num_utetra; ++i) {
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
  bool *should_dealloc_utetra = malloc(num_utetra*sizeof(bool));
  for (int i = 0; i < num_utetra; ++i)
    should_dealloc_utetra[i] = true;

  // See if we can commit a tetrahedron update
  for (int i = 0; i < num_utetra; ++i) {
    if (!isfinite(utetra_get_value(utetra[i])))
      break;
    if (utetra_has_interior_point_solution(utetra[i]) ||
        utetra_has_shadow_solution(utetra[i], eik)) {
      if (utetra_update_ray_is_physical(utetra[i], eik) &&
          commit_tetra_update(eik, l, utetra[i]))
        break;
    } else {
      int num_int = utetra_get_num_interior_coefs(utetra[i]);
      assert(num_int == 1 || num_int == 2);
      int num_adj = 4 - num_int;
      if (i + num_adj <= num_utetra &&
          utetras_yield_same_update((utetra_s const **)&utetra[i], num_adj) &&
          utetra_update_ray_is_physical(utetra[i], eik) &&
          commit_tetra_update(eik, l, utetra[i])) {
        break;
      } else {
        array_append(eik->old_updates, &utetra[i]);
        should_dealloc_utetra[i] = false;
      }
    }
  }

  for (int i = 0; i < num_utetra; ++i)
    if (should_dealloc_utetra[i])
      utetra_dealloc(&utetra[i]);
  free(utetra);

  free(should_dealloc_utetra);

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

static bool cutedge_is_incident_on_vertex(edge_s edge, void const * elt,
                                          void const *aux) {
  (void)elt;
  size_t l1 = *(size_t *)aux;
  return edge.l[0] == l1 || edge.l[1] == l1;
}

bool should_be_in_shadow_zone(eik3_s const *eik, size_t l0, size_t l) {
  edgemap_s *incident_cutedges;
  edgemap_alloc(&incident_cutedges);
  edgemap_init(incident_cutedges, sizeof(cutedge_s));

  edgemap_filter(
    eik->cutset, incident_cutedges,
    (edgemap_prop_t)cutedge_is_incident_on_vertex, &l0);

  assert(!edgemap_is_empty(incident_cutedges));

  edgemap_iter_s *iter;
  edgemap_iter_alloc(&iter);
  edgemap_iter_init(iter, incident_cutedges);

  mesh3_s const *mesh = eik->mesh;
  dbl dy[3], yt[3];

  dbl const *x = mesh3_get_vert_ptr(mesh, l), *y0;

  // Number of supporting shadow boundary planes that say that `x`
  // should be in the shadow zone.
  size_t num_shadow_votes = 0;

  edge_s edge;
  cutedge_s cutedge;
  while (edgemap_iter_next(iter, &edge, &cutedge)) {
    // Get the location of the intersection between the shadow
    // boundary and the current cut edge
    y0 = mesh3_get_vert_ptr(mesh, edge.l[0]);
    dbl3_sub(mesh3_get_vert_ptr(mesh, edge.l[1]), y0, dy);
    dbl3_saxpy(cutedge.t, dy, y0, yt);

    num_shadow_votes += dbl3_dot(cutedge.n, x) - dbl3_dot(cutedge.n, y0) < 0;
  }

  assert(num_shadow_votes == 0 ||
         num_shadow_votes == edgemap_size(incident_cutedges));

  edgemap_iter_dealloc(&iter);

  edgemap_deinit(incident_cutedges);
  edgemap_dealloc(&incident_cutedges);

  return num_shadow_votes > 0;
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
    // Do a triangle update for each neighbor l of the diffracting
    // edge (l0, l1)...
    utri_alloc(&utri[i]);
    utri_reset(utri[i]);
    array_get(nb, i, &l);
    utri_init_from_eik3(utri[i], eik, l, l0, l1);
    utri_solve(utri[i]);

    // ... and attempt to commit it.
    if (utri_has_interior_point_solution(utri[i]) &&
        utri_update_ray_is_physical(utri[i], eik) &&
        commit_tri_update(eik, l, utri[i]))
      adjust(eik, l);
  }

  // Free triangle updates
  for (int i = 0; i < nnb; ++i)
    utri_dealloc(&utri[i]);
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
 * Compute the surface normal for a cutedge attached to a diffracting
 * edge when `l0` has one parent, and where `l0`'s state is
 * `SHADOW`. Since `l0` only has one parent, this will use the
 * diffracting edge to complete the information required to find `n`.
 */
static void get_diff_edge_surf_normal_p1(eik3_s const *eik,
                                         size_t l0, size_t l1,
                                         size_t (*e)[2], size_t ne,
                                         dbl n[3])
{
  dbl const atol = 1e-13;

  mesh3_s const *mesh = eik->mesh;

  dbl (*N)[3] = malloc(ne*sizeof(dbl[3]));

  dbl DT1[3];
  eik3_get_DT(eik, l1, DT1);

  dbl const *x1 = mesh3_get_vert_ptr(mesh, l1);
  dbl dx[3];

  /* Start by compute the normalized cross product between the ray
   * vector at `x1` and each diffracting edge incident on `x1`. */
  for (size_t i = 0; i < ne; ++i) {
    assert(e[i][0] == l1 || e[i][1] == l1);
    if (e[i][0] == l1)
      dbl3_sub(mesh3_get_vert_ptr(mesh, e[i][1]), x1, dx);
    else
      dbl3_sub(mesh3_get_vert_ptr(mesh, e[i][0]), x1, dx);
    dbl3_cross(DT1, dx, N[i]);
    dbl3_normalize(N[i]);
  }

  /* Orient the surface normals consistently. */
  dbl3_sub(mesh3_get_vert_ptr(mesh, l0), x1, dx);
  for (size_t i = 0; i < ne; ++i)
    if (dbl3_dot(N[i], dx) > 0)
      dbl3_negate(N[i]);

  // TODO: for now, ensure that all these normals are the same at this
  // point. This is because we assume that the edges are all entirely
  // straight. Later, when we relax this assumption, we'll need to
  // handle this.
  for (size_t i = 1; i < ne; ++i)
    assert(dbl3_dist(N[i], N[i - 1]) < atol);

  /* Set the surface normal (arbitrarily) to be the first of these
   * surface normals, since we're assuming they're all the same at
   * this point, anyway. */
  dbl3_copy(N[0], n);

  free(N);
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

static void estimate_cutedge_data_from_incident_cutedges(
  mesh3_s const *mesh, size_t l0, size_t l1, cutedge_s *cutedge,
  edgemap_s const *incident_cutedges)
{
  assert(!edgemap_is_empty(incident_cutedges));

  edge_s edge;
  cutedge_s incident_cutedge;

  dbl const *x0, *y0;
  dbl dx[3], dy[3], yt[3];
  x0 = mesh3_get_vert_ptr(mesh, l0);
  dbl3_sub(mesh3_get_vert_ptr(mesh, l1), x0, dx);

  cutedge->t = 0;
  dbl3_zero(cutedge->n);

  edgemap_iter_s *iter;
  edgemap_iter_alloc(&iter);
  edgemap_iter_init(iter, incident_cutedges);

  while (edgemap_iter_next(iter, &edge, &incident_cutedge)) {
    // Get the location of the intersection between the shadow
    // boundary and the current cut edge

    y0 = mesh3_get_vert_ptr(mesh, edge.l[0]);
    dbl3_sub(mesh3_get_vert_ptr(mesh, edge.l[1]), y0, dy);
    dbl3_saxpy(incident_cutedge.t, dy, y0, yt);

    dbl t = dbl3_dot(incident_cutedge.n, yt);
    t -= dbl3_dot(incident_cutedge.n, x0);
    t /= dbl3_dot(incident_cutedge.n, dx);

    cutedge->t += t;

    dbl3_add_inplace(cutedge->n, incident_cutedge.n);
  }

  edgemap_iter_dealloc(&iter);

  size_t num_incident = edgemap_size(incident_cutedges);

  cutedge->t /= num_incident;

  // TODO: this is the perfect place to do a weighted spherical
  // average instead of just normalizing!
  dbl3_normalize(cutedge->n);
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
static bool get_cut_edge_coef_and_surf_normal(eik3_s const *eik,
                                              size_t l0, size_t l1,
                                              cutedge_s *cutedge) {
  /* TODO: pretty sure this function is next on the chopping block.
   * I really don't like how arbitrary the approach to find nearby
   * cutedges to use to estimate the data for the passed cutedge
   * is. Let's see if we can push it a bit further for now,
   * though.
   *
   * Some observations...
   *
   * 1. The "base case" section and "inductive step" sections should
   *    really be split off into two separate functions. */

  //////////////////////////////////////////////////////////////////////////////
  // BASE CASE /////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  /* Setup */

  mesh3_s const *mesh = eik3_get_mesh(eik);

  dbl const atol = 1e-15;

  dbl *t = &cutedge->t;
  dbl *n = cutedge->n;

  // For convenience, get the index with the VALID state...
  size_t l_valid = eik->state[l0] == VALID ? l0 : l1;
  assert(eik->state[l_valid] == VALID);

  // .. and the one with the SHADOW state.
  size_t l_shadow = eik->state[l0] == SHADOW ? l0 : l1;
  assert(eik->state[l_shadow] == SHADOW);

  /**
   * "Base case": check and see if l0's parents are incident on a
   * diffracting edge.
   *
   * What we're checking here:
   * - l1 is one of l0's parents
   * - l0 has exactly two parents (why not just one == l1?)
   * - l2 is the other parent
   * - [l1, l2] is a diffracting edge
   */
  par3_s par = eik3_get_par(eik, l0);
  int npar = par3_size(&par);
  if (npar == 1) {
    size_t num_inc_diff_edges = mesh3_get_num_inc_diff_edges(mesh, l1);
    if (num_inc_diff_edges > 0) {
      assert(l0 == l_shadow && l1 == l_valid);
      assert(num_inc_diff_edges <= 2); // TODO: corners later...
      size_t (*e)[2] = malloc(num_inc_diff_edges*sizeof(size_t[2]));
      mesh3_get_inc_diff_edges(mesh, l1, e);
      get_diff_edge_surf_normal_p1(eik, l0, l1, e, num_inc_diff_edges, n);
      free(e);
      *t = 1;
      return true;
    }
  } else if (npar == 2) {
    // Get indices
    size_t l2 = par.l[0], l3 = par.l[1];
    // Check if l0 was updated from l1
    if (l2 == l1 || l3 == l1) {
      // Swap l2 and l3 so that l1 and l2 are the parents of l0.
      if (l2 == l1) SWAP(l2, l3);
      // Check if (l1, l2) is a diffracting edge. In this case we
      // assume that l1 is the VALID node (since diff => VALID) and
      // return t = 1.
      if (mesh3_is_diff_edge(mesh, (size_t[2]) {l1, l2})) {
        assert(l0 == l_shadow && l1 == l_valid);
        get_diff_edge_surf_normal_p2(eik, l0, (size_t[2]) {l1, l2}, n);
        *t = 1; // This assumes that l1 is valid
        return true;
      }
    }
  }

  // /**
  //  * "Another base case": if a SHADOW node receives its final value
  //  * from a single node on a diffracting edge (TODO: or, later,
  //  * diffracting vertex!), then we won't have enough information to
  //  * compute the surface normal for the cutedge connecting it and its
  //  * parent, since we'll only have ...
  //  */

  bool updated_cutedge = false;

  edgemap_s *incident_cutedges;
  edgemap_alloc(&incident_cutedges);
  edgemap_init(incident_cutedges, sizeof(cutedge_s));

  edgemap_filter(
    eik->cutset, incident_cutedges,
    (edgemap_prop_t)cutedge_is_incident_on_vertex, &l_valid);

  size_t num_incident = edgemap_size(incident_cutedges);

  /**
   * Next, we check if the valid node is incident on a diffracting
   * edge. If it does, we can steal the surface normals that are
   * already computed there. This can be helpful if we're trying to
   * attach a deeper shadow node to a diffracting edge.
   */

  if (mesh3_vert_incident_on_diff_edge(mesh, l_valid)) {
    assert(!edgemap_is_empty(incident_cutedges));

    edge_s edge;
    cutedge_s incident_cutedge;

    edgemap_iter_s *iter;
    edgemap_iter_alloc(&iter);
    edgemap_iter_init(iter, incident_cutedges);

    dbl3_zero(n);
    while (edgemap_iter_next(iter, &edge, &incident_cutedge)) {
      if (eik3_is_valid(eik, edge.l[0]))
        assert(incident_cutedge.t == 0);
      else
        assert(incident_cutedge.t == 1);
      dbl3_add_inplace(n, incident_cutedge.n);
    }

	edgemap_iter_dealloc(&iter);

    // TODO: this is the perfect place to do a weighted spherical
    // average instead of just normalizing!
    dbl3_normalize(cutedge->n);

    assert(l1 == l_valid);
    *t = 1;

    updated_cutedge = true;

    goto cleanup;
  }

  /**
   * "Inductive step": find all of the cutset edges on l1 (there
   * should be some!). They will be valid and should have a surface
   * normal that we can use to extrapolate the shadow boundary in
   * order to find the intersection point.
   */

  // After running this filter, incident_cutedges will contain all of
  // the cutedges incident on either l_valid or l_shadow.
  edgemap_filter(
    eik->cutset, incident_cutedges,
    (edgemap_prop_t)cutedge_is_incident_on_vertex, &l_shadow);
  assert(num_incident <= edgemap_size(incident_cutedges));

  updated_cutedge = !edgemap_is_empty(incident_cutedges);

  if (updated_cutedge)
    estimate_cutedge_data_from_incident_cutedges(
      mesh, l0, l1, cutedge, incident_cutedges);

cleanup:
  edgemap_deinit(incident_cutedges);
  edgemap_dealloc(&incident_cutedges);

  // Verify that t is (up to machine precision) in the interval [0,
  // 1]. If it's slightly outside of the interval due to roundoff
  // error, clamp it back before returning.
  //
  // TODO: seems like the errors will be bigger than machine precision
  // the way we're doing it now, but clamping back to [0, 1] seems to
  // help.
  if (updated_cutedge) {
    if (*t < -atol || 1 + atol < *t)
      log_warn("bad cutset coef: t = %1.3f\n", *t);
    *t = fmax(0, fmin(1, *t));
  }

  return updated_cutedge;
}

static void update_shadow_cutset(eik3_s *eik, size_t l0) {
  // Determine what state we're looking for to find new edges in the
  // shadow cut
  state_e op_state = eik->state[l0] == VALID ? SHADOW : VALID;

  // Find the vertex neighbors of node l0
  int nvv = mesh3_nvv(eik->mesh, l0);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, vv);

  array_s *l1_arr;
  array_alloc(&l1_arr);
  array_init(l1_arr, sizeof(size_t), nvv);

  size_t l1;
  for (int i = 0; i < nvv; ++i) {
    l1 = vv[i];
    if (eik->state[l1] == op_state)
      array_append(l1_arr, &l1);
  }

  // We process the nodes neighboring l0 that have the opposite state
  // using a queue. When we try to compute the data for each new
  // cutedge, because of the way we do it right now, we need to
  // compute this data for each node in a certain order. If we fail
  // the first time, we reinsert the node into the queue and try
  // again.
  edge_s edge;
  while (!array_is_empty(l1_arr)) {
    array_pop_front(l1_arr, &l1);

    cutedge_s cutedge;
    if (!get_cut_edge_coef_and_surf_normal(eik, l0, l1, &cutedge)) {
      array_append(l1_arr, &l1);
      continue;
    }
    assert(0 <= cutedge.t && cutedge.t <= 1);
    assert(fabs(1 - dbl3_norm(cutedge.n)) < 1e-15);

    if (l0 > l1)
      cutedge.t = 1 - cutedge.t;

    edge = make_edge(l0, l1);
    assert(!edgemap_contains(eik->cutset, edge));

    edgemap_set(eik->cutset, edge, &cutedge);
  }

  array_deinit(l1_arr);
  array_dealloc(&l1_arr);

  free(vv);
}

size_t eik3_step(eik3_s *eik) {
  size_t l0 = heap_front(eik->heap);
  assert(eik->state[l0] == TRIAL);
  heap_pop(eik->heap);

  assert(isfinite(eik->jet[l0].f));

  // Once we accept a node, we can clear out the old updates that are
  // targeting it
  utetra_s *old_utetra;
  for (size_t i = array_size(eik->old_updates); i > 0; --i) {
    array_get(eik->old_updates, i - 1, &old_utetra);
    if (utetra_get_l(old_utetra) == l0) {
      utetra_dealloc(&old_utetra);
      array_delete(eik->old_updates, i - 1);
    }
  }

  eik->state[l0] = is_shadow(eik, l0) ? SHADOW : VALID;
  ++eik->num_valid;

  update_shadow_cutset(eik, l0);
  update_neighbors(eik, l0);

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
