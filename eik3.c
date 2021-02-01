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

// static bool can_update_from_face(eik3_s const *eik, size_t const l[3]) {
//   return can_update_from_point(eik, l[0]) &&
//     can_update_from_point(eik, l[1]) &&
//     can_update_from_point(eik, l[2]);
// }

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

// static bool find_l1(eik3_s *eik, size_t l, size_t l0, size_t *l1) {
//   int nvv0 = mesh3_nvv(eik->mesh, l0);
//   size_t *vv0 = malloc(nvv0*sizeof(size_t));
//   mesh3_vv(eik->mesh, l0, vv0);

//   utri_s *utri;
//   utri_alloc(&utri);

//   dbl T = INFINITY;
//   int min_i1 = NO_INDEX;
//   for (int i1 = 0; i1 < nvv0; ++i1) {
//     size_t l1 = vv0[i1];
//     if (!can_update_from_point(eik, l1)) {
//       continue;
//     }
//     if (eik3_is_point_source(eik, l1)) {
//       do_1pt_update(eik, l, l1);
//       utri_dealloc(&utri);
//       free(vv0);
//       return false;
//     }
//     utri_init_from_eik3(utri, eik, l, l0, l1);
//     utri_solve(utri);
//     // TODO: Masha says a better criterion is to set l1 to be the tri
//     // update with an interior point solution (is it unique?)... not
//     // so sure about this...
//     dbl Tnew = utri_get_value(utri);
//     if (Tnew < T) {
//       T = Tnew;
//       min_i1 = i1;
//     }
//   }

//   utri_dealloc(&utri);
//   *l1 = vv0[min_i1];
//   free(vv0);
//   return min_i1 != NO_INDEX;
// }

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

// /**
//  * TODO: We *REALLY* need to reduce some of the redundancy between
//  * `do_tetra_updates` and `do_tetra_update` above.
//  */
// static void do_tetra_updates(eik3_s *eik, size_t l, size_t l0, size_t l1,
//                              size_t const *l2, int n) {
//   assert(!eik3_is_point_source(eik, l0));
//   assert(!eik3_is_point_source(eik, l1));

//   utetra_s **utetra = malloc(n*sizeof(utetra_s *));
//   memset(utetra, 0x0, n*sizeof(utetra_s *));

//   dbl lam[2], b[3];

//   utetra_s *cf;
//   for (int i = 0; i < n; ++i) {
//     if (!can_update_from_point(eik, l2[i])) {
//       continue;
//     }
//     if (eik3_is_point_source(eik, l2[i])) {
//       do_1pt_update(eik, l, l2[i]);
//       goto cleanup;
//     }
//     utetra_alloc(&utetra[i]);
//     cf = utetra[i];
//     utetra_init_from_eik3(cf, eik, l, l0, l1, l2[i]);
//     // TODO: more efficient and simpler to check if x, x0, x1, and x2
//     // are coplanar *before* calling utetra_init_from_eik3 (then we
//     // don't need to dealloc below(
//     if (utetra_is_degenerate(cf)) {
//       utetra_dealloc(&utetra[i]);
//       utetra[i] = NULL;
//       continue;
//     }
//     lam[0] = lam[1] = 0;
//     utetra_set_lambda(cf, lam);
//     utetra_solve(cf);
//   }

//   // Sort the resulting tetrahedron updates. This sorts the update by
//   // their eikonal value.
//   qsort(utetra, n, sizeof(utetra_s *), (compar_t)utetra_cmp);

//   // Traverse the sorted updates, looking for minimum updates. There
//   // may be multiple minimizers. In this case, we look further up in
//   // the order and check if the minimizers and minimizing arguments
//   // match.
//   for (int i = 0; i < n; ++i) {
//     cf = utetra[i];
//     if (cf == NULL || utetra_get_value(cf) > eik->jet[l].f) {
//       break;
//     }

//     if (!utetra_update_ray_is_physical(cf, eik, (size_t[3]) {l0, l1, l2[i]}))
//       continue;

//     utetra_get_bary_coords(cf, b);
//     // TODO: handle an update that is optimal and goes through *three*
//     // faces (i.e. a "vertex" solution)
//     if (utetra_has_interior_point_solution(cf)) {
//       utetra_get_jet(cf, &eik->jet[l]);
//       eik->par[l] = make_par3((size_t[3]) {l0, l1, l2[i]}, b);
//       goto cleanup;
//     }
//     if (i + 1 < n && utetra[i + 1] != NULL &&
//         utetra_adj_are_optimal(cf, utetra[i + 1])) {
//       utetra_get_jet(cf, &eik->jet[l]);
//       eik->par[l] = make_par3((size_t[3]) {l0, l1, NO_PARENT}, b);
//       goto cleanup;
//     }
//   }

// cleanup:
//   for (int i = 0; i < n; ++i) {
//     if (utetra[i]) {
//       utetra_dealloc(&utetra[i]);
//     }
//   }
//   free(utetra);
// }


// static utetra_s *do_tetra_update(eik3_s *eik, size_t l,
//                                  size_t l0, size_t l1, size_t l2) {
//   utetra_s *utetra;
//   utetra_alloc(&utetra);

//   utetra_init_from_eik3(utetra, eik, l, l0, l1, l2);
//   if (utetra_is_degenerate(utetra)) {
//     utetra_dealloc(&utetra);
//     return NULL;
//   }

//   utetra_set_lambda(utetra, (dbl[2]) {0, 0});
//   utetra_solve(utetra);

//   if (utetra_get_value(utetra) >= eik->jet[l].f)
//     goto cleanup;

//   if (!utetra_has_interior_point_solution(utetra))
//     goto cleanup;

//   size_t L[3] = {l0, l1, l2};

//   if (!utetra_update_ray_is_physical(utetra, eik, L))
//     goto cleanup;

//   utetra_get_jet(utetra, &eik->jet[l]);

//   dbl b[3];
//   utetra_get_bary_coords(utetra, b);
//   eik->par[l] = make_par3(L, b);
// }

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
   * Next, we want to see if l0 is incident on any diffracting edges
   * and do these updates (we assume that the only two-point updates
   * that can be done are ones emanating from diffracting edges). Like
   * below, we need to do all of the updates, sort them, and check for
   * boundary cases with nonzero Lagrange multipliers.
   */

  array_s *l1_diff;
  array_alloc(&l1_diff);
  array_init(l1_diff, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  get_valid_incident_diff_edges(eik, l0, l1_diff);

  int num_utri = array_size(l1_diff);
  utri_s **utri = malloc(num_utri*sizeof(utri_s *));

  // Do a triangle update for each diffracting edge and sort
  for (int i = 0; i < num_utri; ++i) {
    utri_alloc(&utri[i]);
    utri_reset(utri[i]);
    size_t l1;
    array_get(l1_diff, i, &l1);
    assert(!eik3_is_point_source(eik, l1));
    utri_init_from_eik3(utri[i], eik, l, l0, l1);
    utri_solve(utri[i]);
  }
  qsort(utri, num_utri, sizeof(utri_s *), (compar_t)utri_cmp);

  // Try to commit a triangle update
  for (int i = 0; i < num_utri; ++i) {
    size_t L[2] = {l0, l1[i]};
    if (!isfinite(utri_get_value(utri[i])))
      break;
    if (utri_has_interior_point_solution(utri[i])) {
      if (utri_update_ray_is_physical(utri[i], eik, L) &&
          commit_tri_update(eik, l, utri[i]))
        break;
    } else if (i + 1 < num_utri &&
               utris_yield_same_update(utri[i], utri[i + 1]) &&
               utri_update_ray_is_physical(utri[i], eik, L) &&
               commit_tri_update(eik, l, utri[i])) {
      break;
    }
  }

  // Free triangle updates
  for (int i = 0; i < num_utri; ++i)
    utri_dealloc(&utri[i]);
  free(utri);

  // Free diff edge indices
  array_deinit(l1_diff);
  array_dealloc(&l1_diff);

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

  // size_t l1;
  // if (!find_l1(eik, l, l0, &l1)) {
  //   return;
  // }

  // // TODO: probably best to just simplify this, using "mesh3_ev"
  // // instead...
  // int nee = mesh3_nee(eik->mesh, (size_t[2]) {l0, l1});
  // size_t (*ee)[2] = malloc(nee*sizeof(size_t[2]));
  // mesh3_ee(eik->mesh, (size_t[2]) {l0, l1}, ee);

  // // TODO: not using l3 right now, but might want to use it later if
  // // we want to try out "volume updates"
  // size_t *l2 = malloc(nee*sizeof(size_t));
  // for (int i = 0; i < nee; ++i)
  //   l2[i] = ee[i][0];

  // // TODO: use initial guess for lambda taken from `do_2pt_updates` to
  // // use as a warm start in `do_tetra_updates`
  // do_tetra_updates(eik, l, l0, l1, l2, nee);

  // // Do any updates corresponding to valid faces surrounding l and
  // // including l0. These may have been missed when doing the preceding
  // // hierarchical update.
  // int nvf = mesh3_nvf(eik->mesh, l);
  // size_t (*vf)[3] = malloc(3*nvf*sizeof(size_t));
  // mesh3_vf(eik->mesh, l, vf);
  // for (int i = 0; i < nvf; ++i) {
  //   if ((l0 == vf[i][0] || l0 == vf[i][1] || l0 == vf[i][2]) &&
  //       can_update_from_face(eik, vf[i])) {
  //     do_tetra_update(eik, l, vf[i]);
  //   }
  // }

  // free(vf);
  // free(l2);
  // free(ee);
}

// static void do_tri_update(eik3_s *eik, size_t l, size_t l0, size_t l1) {
//   dbl const atol = 1e-15;

//   if (!can_update_from_point(eik, l0) || !can_update_from_point(eik, l1)) {
//     return;
//   }

//   // TODO: this stuff here assumes this is a triangle update immersed
//   // in the boundary. Should we really check that in here, or do it
//   // elsewhere?
//   if (!mesh3_bdf(eik->mesh, (size_t[3]) {l, l0, l1})) {
//     return;
//   }
//   assert(mesh3_bdv(eik->mesh, l));
//   assert(mesh3_bdv(eik->mesh, l0));
//   assert(mesh3_bdv(eik->mesh, l1));

//   if (eik3_is_point_source(eik, l0)) {
//     do_1pt_update(eik, l, l0);
//     return;
//   }

//   if (eik3_is_point_source(eik, l1)) {
//     do_1pt_update(eik, l, l1);
//     return;
//   }

//   utri_s *utri;
//   utri_alloc(&utri);
//   utri_init_from_eik3(utri, eik, l, l0, l1);
//   utri_solve(utri);

//   if (utri_get_value(utri) >= eik->jet[l].f)
//     goto cleanup;

//   if (utri_get_lag_mult(utri) > atol)
//     goto cleanup;

//   if (!utri_update_ray_is_physical(utri, eik, (size_t[2]) {l0, l1}))
//     goto cleanup;

//   utri_get_jet(utri, &eik->jet[l]);

//   dbl b[3] = {[2] = 0};
//   utri_get_bary_coords(utri, b);
//   eik->par[l] = make_par3((size_t[3]) {l0, l1, NO_PARENT}, b);

// cleanup:
//   utri_dealloc(&utri);
// }

// /**
//  * Updating boundary points is much trickier than updating interior
//  * points. We handle this case separately.
//  */
// static void update_bd(eik3_s *eik, size_t l, size_t l0) {
//   if (eik3_is_point_source(eik, l0)) {
//     do_1pt_update(eik, l, l0);
//     return;
//   }

//   // Find the adjacent faces
//   //
//   // TODO: we could compute this information on the fly from vc if we
//   // wanted to
//   int nvf = mesh3_nvf(eik->mesh, l);
//   size_t (*vf)[3] = malloc(3*nvf*sizeof(size_t)), *lf;
//   mesh3_vf(eik->mesh, l, vf);

//   // Do adjacent two-point updates that are immersed in the boundary.
//   //
//   // TODO: ideally, we'd just coalesce this into the tetrahedron
//   // updates below with a suitable check for diffracting boundary
//   // solutions
//   for (int i = 0; i < nvf; ++i) {
//     lf = vf[i];
//     do_2pt_bd_update(eik, l, lf[0], lf[1]);
//     do_2pt_bd_update(eik, l, lf[1], lf[2]);
//     do_2pt_bd_update(eik, l, lf[2], lf[0]);
//   }

//   // TODO: okay... trying to think about this a bit more. Here's a new
//   // idea that might work a bit better to combine the bd/non-bd
//   // cases. To determine the updates we want to do, we of course first
//   // have the new VALID index l0. So, let's find the "VALID triangle
//   // fan" surrounding l0...
//   assert(false);

//   /**
//    * Do adjacent tetrahedron updates.
//    */
//   for (int i = 0; i < nvf; ++i) {
//     lf = vf[i];
//     if (point_in_face(l0, lf) && can_update_from_face(eik, lf))
//       do_tetra_update(eik, l, lf);
//   }

//   // Find the adjacent cells
//   int nvc = mesh3_nvc(eik->mesh, l);
//   size_t *vc = malloc(nvc*sizeof(size_t));
//   mesh3_vc(eik->mesh, l, vc);

//   /**
//    * Look for tetrahedron updates on the opposite sides of the first
//    * layer of faces, and do them now
//    */
//   size_t l3;
//   for (int i = 0; i < nvf; ++i) {
//     lf = vf[i];
//     if (!point_in_face(l0, lf))
//       continue;
//     // Get point opposite face
//     if (!mesh3_ccfv(eik->mesh, vc[i], lf, &l3)) {
//       // If this fails, it should be because we tried to find a point
//       // on the other side of a boundary face
//       assert(mesh3_bdf(eik->mesh, lf));
//       continue;
//     }
//     // If pair of tetras is convex
//     if (!adj_tetra_pair_is_convex(eik->mesh, l, lf, l3))
//       continue;
//     // Do each valid tetra on the other side
//     for (int j = 0; j < 3; ++j) {
//       if (lf[j] == l0) // l0 will get swapped out...
//         continue; // ... so move to the next update
//       SWAP(lf[j], l3);
//       if (can_update_from_face(eik, lf))
//         do_tetra_update(eik, l, lf);
//       SWAP(lf[j], l3);
//     }
//   }

//   free(vc);
//   free(vf);
// }

static void adjust(eik3_s *eik, size_t l) {
  assert(eik->state[l] == TRIAL);
  assert(l >= 0);
  assert(l < mesh3_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l]);
}

size_t eik3_peek(eik3_s const *eik) {
  return heap_front(eik->heap);
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

  // Update neighboring nodes.
  for (int i = 0; i < nnb; ++i) {
    if (eik->state[l = nb[i]] == TRIAL) {
      update(eik, l, l0);
      adjust(eik, l);
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

  /**
   * "Inductive step": find all of the cutset edges on l1 (there
   * should be some!). They will be valid and should have a surface
   * normal that we can use to extrapolate the shadow boundary in
   * order to find the intersection point.
   */
  edgemap_s *incident_cutedges;
  edgemap_alloc(&incident_cutedges);
  edgemap_init(incident_cutedges, sizeof(cutedge_s));

  edgemap_filter(eik->cutset, incident_cutedges,
                 (edgemap_prop_t)cutedge_is_incident_on_vertex, &l1);

  assert(!edgemap_is_empty(incident_cutedges));

  edgemap_iter_s *iter;
  edgemap_iter_alloc(&iter);
  edgemap_iter_init(iter, incident_cutedges);

  edge_s edge;
  cutedge_s cutedge;
  while (edgemap_iter_next(iter, &edge, &cutedge)) {
    assert(false);
  }

  edgemap_iter_dealloc(&iter);

  edgemap_deinit(incident_cutedges);
  edgemap_dealloc(&incident_cutedges);

  assert(false);
  return -1;

  // par3_s par = eik3_get_par(eik, l0);
  // int npar = par3_size(&par);


  // // Get all of the cells adjacent to (l0, l1)
  // int nec = mesh3_nec(mesh, l0, l1);
  // size_t *ec = malloc(nec*sizeof(size_t));
  // mesh3_ec(mesh, l0, l1, ec);

  // array_s *valid_edges;
  // array_alloc(&valid_edges);
  // array_init(valid_edges, sizeof(edge_s), ARRAY_DEFAULT_CAPACITY);

  // // Try to find any neighboring valid edges
  // size_t e[2];
  // for (int i = 0; i < nec; ++i) {
  //   mesh3_cee(mesh, ec[i], (size_t[2]) {l0, l1}, e);
  //   if (can_update_from_edge(eik, e)) {
  //     edge_s edge = make_edge(e[0], e[1]);
  //     assert(!array_contains(valid_edges, &edge));
  //     array_append(valid_edges, &edge);
  //   }
  // }

  // int num_valid_edges = array_size(valid_edges);
  // assert(num_valid_edges > 0);

  // // First, check for diffracting edges...
  // int num_valid_diff_edges = 0;
  // for (int i = 0; i < num_valid_edges; ++i) {

  // }

  // array_deinit(valid_edges);
  // array_dealloc(&valid_edges);

  // free(ec);

  // // TODO: this is probably not what we want to do... comment it out
  // // for now
  // // if (mesh3_vert_incident_on_diff_edge(eik->mesh, l_valid))
  // //   return l_valid == l0 ? 0 : 1;

  // par3_s par = eik3_get_par(eik, l0);
  // int npar = par3_size(&par);

  // /**
  //  * In this next rather complicated section we try to figure out what
  //  * t by looking at the parents of l0.
  //  */
  // dbl t;
  // if (npar == 1) {
  //   size_t l2 = par.l[0];
  //   (void)l2;

  //   // TODO: implement this

  //   assert(false);
  // } else if (npar == 2) {
  //   // Get indices
  //   size_t l2 = par.l[0], l3 = par.l[1];

  //   // Check if l0 was updated from l1
  //   if (l2 == l1 || l3 == l1) {
  //     if (l2 == l1) SWAP(l2, l3);

  //     // Check if (l1, l2) is a diffracting edge. In this case we
  //     // assume that l1 is the VALID node (since diff => VALID) and
  //     // return t = 1.
  //     if (mesh3_is_diff_edge(eik->mesh, (size_t[2]) {l1, l2})) {
  //       assert(l1 == l_valid);
  //       return 1;
  //     } else {
  //       // TODO: not 100% sure what this case corresponds to.
  //       assert(false);
  //     }
  //   } else {
  //     // First, get grad T at each parent and compute their cross
  //     // product.
  //     dbl DT2[3], DT3[3], n[3];
  //     eik3_get_DT(eik, l2, DT2);
  //     eik3_get_DT(eik, l3, DT3);
  //     dbl3_cross(DT2, DT3, n);

  //     // Compute some intermediate vector quantities that we need to
  //     // find where the plane defined by n and x2 intersects [x0, x1].
  //     dbl const *x0 = mesh3_get_vert_ptr(eik->mesh, l0);
  //     dbl dx1[3], dx2[3];
  //     dbl3_sub(mesh3_get_vert_ptr(eik->mesh, l1), x0, dx1);
  //     dbl3_sub(mesh3_get_vert_ptr(eik->mesh, l2), x0, dx2);

  //     // Find the coefficient parametrizing the intersection between
  //     // [x0, x1] and {x: n'*(x - x2) = 0}.
  //     t = dbl3_dot(n, dx2)/dbl3_dot(n, dx1);
  //     assert(0 <= t && t <= 1);
  //     return t;
  //   }
  // } else if (npar == 3) {


  //   assert(false);
  // } else {
  //   assert(false);
  // }
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

    // TODO: we want to not only insert t here, but also the surface
    // normal at this point. This will give us a better jumping-off
    // point when trying to insert new cutset edges.
    //
    // Here's an idea for an algorithm:
    // 1) We found a new cutset edge (l0, l1)
    // 2) Find the cutset edges (l1, l2) & (xt, n(xt))
    // 3) All cutset edges should be upwind (assert)
    // 4) Try to extrapolate the surface along from (xt, n(xt))
    // 5) Can we do some kind of averaging to make this more stable?
    //    Is there some nonlinear equation that a point on the shadow
    //    boundary should satisfy?
    //
    // Also, see page 123 of Farin's Tri BB paper for more about
    // computing parametric 9-parameter interpolants from three points
    // and three normal vectors
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
