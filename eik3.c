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

struct eik3 {
  mesh3_s const *mesh;
  jet3 *jet;
  state_e *state;
  int *pos;
  par3_s *par;
  heap_s *heap;
  int num_valid;
  edgemap_s *cutset;
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

void eik3_init(eik3_s *eik, mesh3_s const *mesh) {
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
}

static cutedge_s get_cutedge(eik3_s const *eik, size_t l0, size_t l1) {
  cutedge_s cutedge;
  bool found = edgemap_get(eik->cutset, make_edge(l0, l1), &cutedge);
  assert(found);

  // Reverse `cutedge.t` if we need to.
  if (l0 > l1)
    cutedge.t = 1 - cutedge.t;

  return cutedge;
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
  assert(num_shadow >= 0);
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
    if (check_state == VALID)
      return t[0] > 0 || t[1] > 0;
    else
      assert(false);
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
    assert(false);

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
  return eik->state[l] == VALID || eik->state[l] == SHADOW;
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
  eik->jet[l] = jet;
  eik->par[l] = (par3_s) {
    .l = {l0, NO_PARENT, NO_PARENT},
    .b = {1, NAN, NAN}
  };
}

static bool commit_tri_update(eik3_s *eik, size_t lhat, utri_s const *utri) {
  if (utri_get_value(utri) >= eik->jet[lhat].f)
    return false;

  utri_get_jet(utri, &eik->jet[lhat]);

  size_t l[3] = {[2] = NO_PARENT};
  utri_get_update_inds(utri, l);

  dbl b[3] = {[2] = NAN};
  utri_get_bary_coords(utri, b);

  eik->par[lhat] = make_par3(l, b);

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
    return ntri;

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

  utetra_get_jet(utetra, &eik->jet[lhat]);

  size_t l[3];
  utetra_get_update_inds(utetra, l);

  dbl b[3];
  utetra_get_bary_coords(utetra, b);

  eik->par[lhat] = make_par3(l, b);

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

  // TODO: got through l1 and l2 and check if they're point
  // sources. If they are, we need to do the point source update and
  // return.
  for (int i = 0; i < num_utetra; ++i) {
    if (eik3_is_point_source(eik, l1[i])) {
      do_1pt_update(eik, l, l1[i]);
      goto cleanup;
    }
    if (eik3_is_point_source(eik, l2[i])) {
      do_1pt_update(eik, l, l2[i]);
      goto cleanup;
    }
  }

  /**
   * Now we move on to doing tetrahedron updates
   */

  // Allocate tetrahedron updates
  utetra_s **utetra = malloc(num_utetra*sizeof(utetra_s *));

  // Do each tetrahedron update and sort
  for (int i = 0; i < num_utetra; ++i) {
    utetra_alloc(&utetra[i]);
    utetra_reset(utetra[i]);
    if (eik3_is_point_source(eik, l1[i]) || eik3_is_point_source(eik, l2[i]))
      continue;
    if (!utetra_init_from_eik3(utetra[i], eik, l, l0, l1[i], l2[i]))
      continue;
    if (utetra_is_degenerate(utetra[i]))
      continue;
    utetra_set_lambda(utetra[i], (dbl[2]) {0, 0});
    utetra_solve(utetra[i]);
  }
  qsort(utetra, num_utetra, sizeof(utetra_s *), (compar_t)utetra_cmp);

  // See if we can commit a tetrahedron update
  for (int i = 0; i < num_utetra; ++i) {
    if (!isfinite(utetra_get_value(utetra[i])))
      break;
    if (utetra_has_interior_point_solution(utetra[i])) {
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
      }
    }
  }

  for (int i = 0; i < num_utetra; ++i)
    utetra_dealloc(&utetra[i]);
  free(utetra);

cleanup:
  free(l1);
  free(l2);
}

static void adjust(eik3_s *eik, size_t l) {
  assert(eik->state[l] == TRIAL);
  assert(l >= 0);
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

#if 0
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
#endif

#if 0
/**
 * We assume that `l0` is a node in the shadow zone. This function
 * first predicts whether the `TRIAL` node `l` is in the shadow
 * zone. If it is, we look for `FAR` neighbors of both `l0` and `l`
 * that also (seem!) to be in the shadow zone. We then "pull" updates
 * for these neighbors. This may happen recursively.
 */
void pull_shadow_updates(eik3_s *eik, size_t l0, size_t l,
                         bool check_shadow) {
  assert(eik->state[l0] == SHADOW);
  assert(eik->state[l] == TRIAL);

  if (check_shadow && !should_be_in_shadow_zone(eik, l0, l))
    return;

  int nvv = mesh3_nvv(eik->mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(eik->mesh, l, vv);

  size_t l_;
  for (int j = 0; j < nvv; ++j) {
    if (eik->state[l_ = vv[j]] == FAR) {
      eik->state[l_] = SHADOW;
      // TODO: do this recursively!!!
      update(eik, l_, l0);
    }
  }

  free(vv);
}
#endif

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

  int nnb = array_size(nb);
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
#if 0
  else {
    // When we update the TRIAL neighbors of a newly accepted SHADOW
    // node, we want to make sure we have all the information we need
    // to be able to do the update reasonably accurately. We don't
    // stage the FAR neighbors of a newly accepted SHADOW node as a
    // matter of course. Instead, we look for TRIAL neighbors and
    // stage FAR nodes which are neighbors of the TRIAL neighbor *and*
    // the newly accepted SHADOW node.
    for (int i = 0; i < nnb; ++i)
      if (eik->state[l = nb[i]] == TRIAL)
        pull_shadow_updates(eik, l0, l, true);
  }
#endif

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
 * Compute the surface normal at a diffraction edge. The index `l0`
 * indicates a new `SHADOW` point, which we use to orient the surface
 * normal, computed as the cross product of the eikonal gradients at
 * indices `l[0]` and `l[1]`. The result goes in `n`.
 */
static void get_diff_edge_surf_normal(eik3_s const *eik, size_t l0, size_t l[2],
                                      dbl n[3]) {
  assert(eik3_is_shadow(eik, l0));

  // Get DT at each parent and compute their cross product.
  dbl DT[2][3];
  eik3_get_DT(eik, l[0], DT[0]);
  eik3_get_DT(eik, l[1], DT[1]);
  dbl3_cross(DT[0], DT[1], n);
  dbl3_normalize(n);

  mesh3_s const *mesh = eik3_get_mesh(eik);

  // Reorient the surface normal by checking which side of the tangent
  // place at x[l[0]] the point x[l0] is on.
  dbl dx[3];
  dbl3_sub(mesh3_get_vert_ptr(mesh, l[0]), mesh3_get_vert_ptr(mesh, l0), dx);
  if (dbl3_dot(n, dx) < 0)
    dbl3_negate(n);
}

static void estimate_t_and_normal_from_cut_edges(eik3_s const *eik,
                                                 size_t l0, size_t l1,
                                                 edge_s const *edge,
                                                 cutedge_s const *cutedge,
                                                 int num_incident,
                                                 dbl *that, dbl normal[3]) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  dbl *t = malloc(num_incident*sizeof(dbl));

  // TODO: unnecessary! we're just copying this...
  dbl (*n)[3] = malloc(num_incident*sizeof(dbl[3]));

  dbl const *x0, *y0;
  dbl dx[3], dy[3], yt[3];
  x0 = mesh3_get_vert_ptr(mesh, l0);
  dbl3_sub(mesh3_get_vert_ptr(mesh, l1), x0, dx);

  for (int i = 0; i < num_incident; ++i) {
    // Get the location of the intersection between the shadow
    // boundary and the current cut edge
    y0 = mesh3_get_vert_ptr(mesh, edge[i].l[0]);
    dbl3_sub(mesh3_get_vert_ptr(mesh, edge[i].l[1]), y0, dy);
    dbl3_saxpy(cutedge[i].t, dy, y0, yt);

    memcpy(n[i], cutedge[i].n, sizeof(dbl[3]));
    t[i] = (dbl3_dot(n[i], yt) - dbl3_dot(n[i], x0))/dbl3_dot(n[i], dx);
  }

  // printf("\n"
  //        "estimate_t_and_normal_from_cut_edges(%lu, %lu)\n"
  //        "\tl0\tl1\tt (in)\tn[0]\tn[1]\tn[2]\tt (out)\n",
  //        l0, l1);
  // for (int i = 0; i < num_incident; ++i)
  //   printf("\t%lu\t%lu\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n",
  //          edge[i].l[0], edge[i].l[1], cutedge[i].t, cutedge[i].n[0],
  //          cutedge[i].n[1], cutedge[i].n[2], t[i]);

  /**
   * To estimate t and the normal vector, we just compute the average
   * over each normal and intersection point...
   */
  // TODO: we can probably come up with a smarter way to do this

  *that = dblN_mean(t, num_incident);

  normal[0] = normal[1] = normal[2] = 0;
  for (int i = 0; i < num_incident; ++i) {
    normal[0] += n[i][0];
    normal[1] += n[i][1];
    normal[2] += n[i][2];
  }
  normal[0] /= num_incident;
  normal[1] /= num_incident;
  normal[2] /= num_incident;

  // printf("=> normal = (%1.3f, %1.3f, %1.3f), t = %1.3f\n"
  //        "\n",
  //        normal[0], normal[1], normal[2], *that);

  free(n);
  free(t);
}

/**
 * Compute the coefficient for the new edge in shadow cutset. This is
 * a double t such that 0 <= t <= 1 and where the shadow boundary
 * (approximately) passes through (1 - t)*x[l0] + t*x[l1].
 *
 * The index l0 corresponds to the node that has just been accepted,
 * and l1 is some neighbor of l0 which has the "opposite" state from
 * l0 (i.e., VALID if l0 is SHADOW and vice versa).
 *
 * We assume that one of eik->state[l0] and eik->state[l1] is VALID
 * and the other is SHADOW. It doesn't matter which is which.
 */
static dbl get_cut_edge_coef_and_surf_normal(eik3_s const *eik,
                                             size_t l0, size_t l1,
                                             dbl normal[3]) {
  dbl const atol = 1e-15;

  mesh3_s const *mesh = eik3_get_mesh(eik);

  // For convenience, get the index with the VALID state...
  size_t l_valid = eik->state[l0] == VALID ? l0 : l1;
  assert(eik->state[l_valid] == VALID);

  // .. and the one with the SHADOW state.
  size_t l_shadow = eik->state[l0] == SHADOW ? l0 : l1;
  assert(eik->state[l_shadow] == SHADOW);

  /**
   * "Base case": check and see if l0's parents are incident on a
   * diffracting edge.
   */
  par3_s par = eik3_get_par(eik, l0);
  int npar = par3_size(&par);
  if (npar == 2) {
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
        assert(l1 == l_valid);
        get_diff_edge_surf_normal(eik, l0, (size_t[2]) {l1, l2}, normal);
        return 1;
      }
    }
  }

  dbl t;

  edgemap_s *incident_cutedges;
  edgemap_alloc(&incident_cutedges);
  edgemap_init(incident_cutedges, sizeof(cutedge_s));

  /**
   * Next, we check if the valid node is incident on a diffracting
   * edge. If it does, we can steal the surface normals that are
   * already computed there. This can be helpful if we're trying to
   * attach a deeper shadow node to a diffracting edge.
   */

  if (mesh3_vert_incident_on_diff_edge(mesh, l_valid)) {
    edgemap_filter(
      eik->cutset, incident_cutedges,
      (edgemap_prop_t)cutedge_is_incident_on_vertex, &l_valid);

    // TODO: safeguard here for now...
    assert(!edgemap_is_empty(incident_cutedges));

    t = 1;

    edgemap_iter_s *iter;
    edgemap_iter_alloc(&iter);
    edgemap_iter_init(iter, incident_cutedges);

    normal[0] = normal[1] = normal[2] = 0;

    cutedge_s cutedge;
    while (edgemap_iter_next(iter, NULL, &cutedge))
      for (int i = 0; i < 3; ++i)
        normal[i] += cutedge.n[i];

    int num_incident = edgemap_size(incident_cutedges);
    for (int i = 0; i < 3; ++i)
      normal[i] /= num_incident;

    goto cleanup;
  }

  /**
   * "Inductive step": find all of the cutset edges on l1 (there
   * should be some!). They will be valid and should have a surface
   * normal that we can use to extrapolate the shadow boundary in
   * order to find the intersection point.
   */

  /**
   * Try to compute `t` and `normal` using cutset edges incident on
   * `l_shadow`.
   */

  // Filter out the edges that aren't incident on the SHADOW vertex
  edgemap_clear(incident_cutedges);
  edgemap_filter(
    eik->cutset, incident_cutedges,
    (edgemap_prop_t)cutedge_is_incident_on_vertex, &l_shadow);

  if (!edgemap_is_empty(incident_cutedges)) {
    int num_incident = edgemap_size(incident_cutedges);
    edge_s *edge = malloc(num_incident*sizeof(edge_s));
    cutedge_s *cutedge = malloc(num_incident*sizeof(cutedge_s));

    // TODO: replace this with a new "edgemap_map" function
    edgemap_iter_s *iter;
    edgemap_iter_alloc(&iter);
    edgemap_iter_init(iter, incident_cutedges);

    for (int i = 0; i < num_incident; ++i)
      edgemap_iter_next(iter, &edge[i], &cutedge[i]);
    estimate_t_and_normal_from_cut_edges(
      eik, l0, l1, edge, cutedge, num_incident, &t, normal);

    edgemap_iter_dealloc(&iter);

    free(cutedge);
    free(edge);

    goto cleanup;
  }

  /**
   * If that didn't work, try to compute `t` and `normal` using cutset
   * edges incident on `l_valid`.
   */

  // Filter out the edges that aren't incident on the SHADOW vertex
  edgemap_clear(incident_cutedges);
  edgemap_filter(
    eik->cutset, incident_cutedges,
    (edgemap_prop_t)cutedge_is_incident_on_vertex, &l_valid);
  assert(!edgemap_is_empty(incident_cutedges));

  int num_incident = edgemap_size(incident_cutedges);
  edge_s *edge = malloc(num_incident*sizeof(edge_s));
  cutedge_s *cutedge = malloc(num_incident*sizeof(cutedge_s));

  // TODO: replace this with a new "edgemap_map" function
  edgemap_iter_s *iter;
  edgemap_iter_alloc(&iter);
  edgemap_iter_init(iter, incident_cutedges);

  edgemap_iter_init(iter, incident_cutedges);
  for (int i = 0; i < num_incident; ++i)
    edgemap_iter_next(iter, &edge[i], &cutedge[i]);
  estimate_t_and_normal_from_cut_edges(
    eik, l0, l1, edge, cutedge, num_incident, &t, normal);

  edgemap_iter_dealloc(&iter);

  free(cutedge);
  free(edge);

cleanup:
  edgemap_deinit(incident_cutedges);
  edgemap_dealloc(&incident_cutedges);

  // Verify that t is (up to machine precision) in the interval [0,
  // 1]. If it's slightly outside of the interval due to roundoff
  // error, clamp it back before returning.
  if (t < -atol || 1 + atol < t) {
    log_warn("bad cutset coef: l0 = %lu, l1 = %lu, t = %1.3f\n",
             l0, l1, t);
  }
  t = fmax(0, fmin(1, t));

  return t;
}

static void update_shadow_cutset(eik3_s *eik, size_t l0) {
  // Determine what state we're looking for to find new edges in the
  // shadow cut
  state_e op_state = eik->state[l0] == VALID ? SHADOW : VALID;

  // Find the vertex neighbors of node l0
  int nvv = mesh3_nvv(eik->mesh, l0);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, vv);

  size_t l1;
  edge_s edge;
  for (int i = 0; i < nvv; ++i) {
    if (eik->state[l1 = vv[i]] != op_state)
      continue;

    cutedge_s cutedge;
    cutedge.t = get_cut_edge_coef_and_surf_normal(eik, l0, l1, cutedge.n);
    assert(0 <= cutedge.t && cutedge.t <= 1);

    // Reverse the coefficient here if necessary. We want `t = 0` to
    // correspond to the `SHADOW` node and `t = 1` to the `VALID`
    // node.
    if (l0 > l1)
      cutedge.t = 1 - cutedge.t;

    edge = make_edge(l0, l1);
    assert(!edgemap_contains(eik->cutset, edge));

    printf("added shadow cutset edge: l0 = %lu (%s), l1 = %lu (%s), t = %g, ",
           edge.l[0], eik->state[edge.l[0]] == VALID ? "VALID" : "SHADOW",
           edge.l[1], eik->state[edge.l[1]] == VALID ? "VALID" : "SHADOW",
           cutedge.t);
    printf("n = (%f, %f, %f)\n", cutedge.n[0], cutedge.n[1], cutedge.n[2]);

    edgemap_set(eik->cutset, edge, &cutedge);
  }

  free(vv);
}

size_t eik3_step(eik3_s *eik) {
  size_t l0 = heap_front(eik->heap);
  assert(eik->state[l0] == TRIAL);
  heap_pop(eik->heap);

  assert(jet3_is_finite(&eik->jet[l0]));

  eik->state[l0] = is_shadow(eik, l0) ? SHADOW : VALID;
  ++eik->num_valid;

  update_shadow_cutset(eik, l0);

  // TODO: had `update_neighbors` before `update_shadow_cutset`, but
  // trying this now, so that we can use info from the updated shadow
  // cutset to decide when to reach back
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

mesh3_s const *eik3_get_mesh(eik3_s const *eik) {
  return eik->mesh;
}

jet3 eik3_get_jet(eik3_s const *eik, size_t l) {
  return eik->jet[l];
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

void eik3_get_DT(eik3_s const *eik, size_t l, dbl DT[3]) {
  memcpy(DT, &eik->jet[l].fx, 3*sizeof(dbl));
}

edgemap_s const *eik3_get_cutset(eik3_s const *eik) {
  return eik->cutset;
}
