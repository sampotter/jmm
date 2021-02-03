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

static bool is_shadow(eik3_s const *eik, size_t l0) {
  if (eik3_is_point_source(eik, l0)) {
    return false;
  }

  dbl const atol = 1e-14;

  par3_s const *par = &eik->par[l0];

  // Get the number of parents
  int n = par3_size(par);
  assert(n == 1 || n == 2 || n == 3);

  // Get the information about the parents
  struct {
    int i;
    size_t l;
    dbl b;
    bool active, bd, diff;
  } p[n];
  for (int i = 0; i < n; ++i) {
    p[i].i = i;
    p[i].l = par->l[i];
    p[i].b = par->b[i];
    p[i].active = p[i].b > atol;
    p[i].bd = mesh3_bdv(eik->mesh, p[i].l);
    p[i].diff = eik->state[p[i].l] == SHADOW;
  }

  // Sort the parent info structs so that the first entries are
  // active.
  //
  // TODO: check at some point if this is actually necessary...
  if (n == 2) {
    if (!p[0].active && p[1].active) {
      SWAP(p[0], p[1]);
    }
  } else if (n == 3) {
    if (p[0].active) {
      if (!p[1].active) {
        SWAP(p[1], p[2]);
      }
    } else { // !p[0].active
      if (p[1].active) {
        int i = p[2].active ? 2 : 1;
        SWAP(p[0], p[i]);
      } else {
        SWAP(p[0], p[2]);
      }
    }
  }

  int num_active = 0;
  for (int i = 0; i < n; ++i) num_active += p[i].active;

  int num_diff = 0;
  for (int i = 0; i < n; ++i) num_diff += p[i].diff;

  // Check if the active indices are a subset of the boundary indices.
  bool active_is_subset_of_bd = true;
  for (int i = 0; i < num_active; ++i) {
    if (!p[i].bd) {
      active_is_subset_of_bd = false;
      break;
    }
  }

  // Check if the active indices are a subset of the indices whose
  // parents were themselves reached by diffracted rays. If they are,
  // this is definitely a point with a diffracted ray leading to it,
  // and we can return early.
  bool active_is_subset_of_diff = true;
  for (int i = 0; i < num_active; ++i) {
    if (!p[i].diff) {
      active_is_subset_of_diff = false;
      break;
    }
  }
  if (active_is_subset_of_diff) {
    return true;
  }

  // Now we need to go through and check more carefully whether the
  // ray leading to this point has been diffracted or not.
  //
  // TODO: this is a mess... once this stabilizes, we really want to
  // try and extract the essence of what's going on and simplify this
  // code.
  if (active_is_subset_of_bd) {
    if (num_active == 1) {
      if (eik3_is_point_source(eik, p[0].l)) {
        return false;
      } else if (p[0].diff) {
        assert(false);
      } else if (p[0].bd) {
        dbl dot = dbl3_dot(&eik->jet[l0].fx, &eik->jet[p[0].l].fx);
        return fabs(1 - dot) > atol;
      }
    } else if (num_active == 2) {
      size_t const e[2] = {p[0].l, p[1].l};
      bool bde = mesh3_bde(eik->mesh, e);
      if (!bde) {
        return false;
      }
      if (p[0].diff && p[1].diff) {
        return true;
      } else if (p[0].diff) {
        return !ray_and_face_are_coplanar(
          eik->mesh, l0, p[0].l, p[1].l, &eik->jet[p[1].l].fx);
      } else if (p[1].diff) {
        return !ray_and_face_are_coplanar(
          eik->mesh, l0, p[0].l, p[1].l, &eik->jet[p[0].l].fx);
      } else if (mesh3_is_diff_edge(eik->mesh, e)) {
        dbl n[3];
        dbl3_cross(&eik->jet[p[0].l].fx, &eik->jet[p[1].l].fx, n);
        dbl dot = dbl3_dot(&eik->jet[l0].fx, n);
        return fabs(dot) > atol;
      } else {
        return false;
      }
    } else if (num_active == 3) {
      bool bdf = mesh3_bdf(eik->mesh, (size_t[3]) {p[0].l, p[1].l, p[2].l});
      assert(!bdf);
      return num_diff > 0;
    } else {
      assert(false);
    }
  } else {
    return false;
  }

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

  // Do a triangle update for each neighbor l of the diffracting edge
  // (l0, l1) and sort the completed updates
  for (int i = 0; i < nnb; ++i) {
    utri_alloc(&utri[i]);
    utri_reset(utri[i]);
    array_get(nb, i, &l);
    utri_init_from_eik3(utri[i], eik, l, l0, l1);
    utri_solve(utri[i]);
  }
  qsort(utri, nnb, sizeof(utri_s *), (compar_t)utri_cmp);

  // Try to commit the triangle updates
  for (int i = 0; i < nnb; ++i) {
    if (!isfinite(utri_get_value(utri[i])))
      break;
    array_get(nb, i, &l);
    if (utri_has_interior_point_solution(utri[i])) {
      if (utri_update_ray_is_physical(utri[i], eik) &&
          commit_tri_update(eik, l, utri[i])) {
        adjust(eik, l);
        break;
      }
    } else if (i + 1 < nnb &&
               utris_yield_same_update(utri[i], utri[i + 1]) &&
               utri_update_ray_is_physical(utri[i], eik) &&
               commit_tri_update(eik, l, utri[i])) {
      adjust(eik, l);
      break;
    }
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

void update_neighbors(eik3_s *eik, size_t l0, bool stage_neighbors) {
  size_t l; // Node l is a neighbor of node l0

  // Get i0's neighboring nodes.
  int nnb = mesh3_nvv(eik->mesh, l0);
  size_t *nb = malloc(nnb*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, nb);

  // Set FAR nodes to TRIAL and insert them into the heap.
  if (stage_neighbors) {
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

static bool cutedge_is_incident_on_vertex(edge_s edge, void const * elt,
                                          void const *aux) {
  (void)elt;
  size_t l1 = *(size_t *)aux;
  return edge.l[0] == l1 || edge.l[1] == l1;
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
  // mesh3_s const *mesh = eik3_get_mesh(eik);

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
      if (mesh3_is_diff_edge(eik->mesh, (size_t[2]) {l1, l2})) {
        assert(l1 == l_valid);
        get_diff_edge_surf_normal(eik, l0, (size_t[2]) {l1, l2}, normal);
        return 1;
      }
    }
  }

  dbl t;

  /**
   * "Inductive step": find all of the cutset edges on l1 (there
   * should be some!). They will be valid and should have a surface
   * normal that we can use to extrapolate the shadow boundary in
   * order to find the intersection point.
   */
  edgemap_s *incident_cutedges;
  edgemap_alloc(&incident_cutedges);
  edgemap_init(incident_cutedges, sizeof(cutedge_s));

  // Filter out the edges that aren't incident on the SHADOW vertex
  edgemap_filter(
    eik->cutset, incident_cutedges,
    (edgemap_prop_t)cutedge_is_incident_on_vertex, &l_shadow);

  int num_incident = edgemap_size(incident_cutedges);
  assert(num_incident > 0);

  edgemap_iter_s *iter;
  edgemap_iter_alloc(&iter);
  edgemap_iter_init(iter, incident_cutedges);

  edge_s *edge = malloc(num_incident*sizeof(edge_s));
  cutedge_s *cutedge = malloc(num_incident*sizeof(cutedge_s));
  for (int i = 0; i < num_incident; ++i)
    edgemap_iter_next(iter, &edge[i], &cutedge[i]);
  estimate_t_and_normal_from_cut_edges(
    eik, l0, l1, edge, cutedge, num_incident, &t, normal);

  free(edge);
  free(cutedge);
  edgemap_iter_dealloc(&iter);

  edgemap_deinit(incident_cutedges);
  edgemap_dealloc(&incident_cutedges);

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

    // Reverse the coefficient here if necessary.
    edge = make_edge(l0, l1);
    if (edge.l[0] != l0)
      cutedge.t = 1 - cutedge.t;

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

  assert(isfinite(eik->jet[l0].f));

  eik->state[l0] = is_shadow(eik, l0) ? SHADOW : VALID;
  ++eik->num_valid;

  bool should_stage = eik->state[l0] == VALID;
  update_neighbors(eik, l0, should_stage);

  update_shadow_cutset(eik, l0);

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
