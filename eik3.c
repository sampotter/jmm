#include "eik3.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bb.h"
#include "heap.h"
#include "macros.h"
#include "mat.h"
#include "mesh3.h"
#include "update.h"
#include "vec.h"

struct eik3 {
  mesh3_s const *mesh;
  jet3 *jet;
  state_e *state;
  int *pos;
  heap_s *heap;

  int num_valid;

  int *full_update;
  int num_full_updates;
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

  /**
   * When we compute the initial heap capacity, we want to estimate
   * the number of nodes that could comprise the expanding numerical
   * front at any one time. We can't know this ahead of time, so we
   * set it to a constant multiple times (# nodes)^(1/d). In this
   * case, d=3. Even if this is an underestimate, well still reduce
   * the number of times the heap needs to be expanded at solve time.
   */
  int capacity = (int) 3*cbrt(nverts);

  heap_alloc(&eik->heap);
  heap_init(eik->heap, capacity, value, setpos, (void *)eik);

  eik->num_valid = 0;

  eik->full_update = malloc(nverts*sizeof(int));
  for (size_t l = 0; l < nverts; ++l) {
    eik->full_update[l] = false;
  }

  eik->num_full_updates = 0;
}

void eik3_deinit(eik3_s *eik) {
  free(eik->jet);
  eik->jet = NULL;

  free(eik->state);
  eik->state = NULL;

  free(eik->pos);
  eik->pos = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);

  free(eik->full_update);
  eik->full_update = NULL;
}

static jet3 solve_2pt_bvp(eik3_s const *eik, size_t l, size_t l0) {
  dbl const *x = mesh3_get_vert_ptr(eik->mesh, l);
  dbl const *x0 = mesh3_get_vert_ptr(eik->mesh, l0);
  dbl n0[3];
  dbl3_sub(x, x0, n0);
  dbl L = dbl3_normalize(n0);
  return (jet3) {.f = L, .fx = n0[0], .fy = n0[1], .fz = n0[2]};
}

static void do_1pt_update(eik3_s *eik, size_t l, size_t l0) {
  jet3 jet = solve_2pt_bvp(eik, l, l0);
  assert(jet.f <= eik->jet[l].f);
  eik->jet[l] = jet;
}

static bool do_2pt_updates(eik3_s *eik, size_t l, size_t l0, size_t *l1) {
  int nvv0 = mesh3_nvv(eik->mesh, l0);
  size_t *vv0 = malloc(nvv0*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, vv0);

  utri_s *utri;
  utri_alloc(&utri);

  dbl T = INFINITY;
  int min_i1 = NO_INDEX;
  for (int i1 = 0; i1 < nvv0; ++i1) {
    size_t l1 = vv0[i1];
    if (eik->state[l1] != VALID) {
      continue;
    }
    if (eik3_is_point_source(eik, l1)) {
      do_1pt_update(eik, l, l1);
      utri_dealloc(&utri);
      free(vv0);
      return false;
    }
    utri_init(utri, eik->mesh, eik->jet, l, l0, l1);
    utri_solve(utri);
    dbl Tnew = utri_get_value(utri);
    if (Tnew < T) {
      T = Tnew;
      min_i1 = i1;
    }
  }

  utri_dealloc(&utri);
  free(vv0);
  *l1 = vv0[min_i1];
  return min_i1 != NO_INDEX;
}

static void get_opposite_cell_edge(mesh3_s const *mesh,
                                   size_t c,
                                   size_t l0, size_t l1,
                                   size_t *l2, size_t *l3) {
  size_t v[4];
  mesh3_cv(mesh, c, v);
  int k = 0;
  for (int i = 0; i < 4; ++i) {
    if (v[i] == l0 || v[i] == l1) {
      continue;
    }
    if (k == 0) {
      *l2 = v[i];
      ++k;
    } else if (k == 1) {
      *l3 = v[i];
      ++k;
    } else {
      assert(false);
    }
  }
}

static void sort_and_orient(size_t *l2, size_t *l3, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      // We don't want any duplicate entries among the l2 or l3
      // indices, but this can easily happen. This is the "orient"
      // part of this function.
      if (l2[j] == l2[i] || l3[j] == l3[i]) {
        SWAP(l2[j], l3[j]);
      }
      if (l3[j] == l2[i]) {
        SWAP(l2[i + 1], l2[j]);
        SWAP(l3[i + 1], l3[j]);
      }
    }
  }
}

static int get_l2_and_l3(eik3_s const *eik, size_t l0, size_t l1,
                         size_t **l2_ptr, size_t **l3_ptr) {
  int nec = mesh3_nec(eik->mesh, l0, l1);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(eik->mesh, l0, l1, ec);

  size_t *l2 = *l2_ptr = malloc(nec*sizeof(size_t));
  size_t *l3 = *l3_ptr = malloc(nec*sizeof(size_t));

  for (int i = 0; i < nec; ++i) {
    size_t v[4];
    mesh3_cv(eik->mesh, ec[i], v);
    get_opposite_cell_edge(eik->mesh, ec[i], l0, l1, &l2[i], &l3[i]);
  }

  sort_and_orient(l2, l3, nec);

  free(ec);

  return nec;
}

static void do_tetra_updates(eik3_s *eik, size_t l, size_t l0, size_t l1,
                             size_t const *l2, int n) {
  utetra_s *utetra;
  utetra_alloc(&utetra);

  dbl lam[2], alpha[3];
  jet3 jet;

  for (int i = 0; i < n; ++i) {
    if (eik->state[l2[i]] != VALID) {
      continue;
    }
    assert(!eik3_is_point_source(eik, l0));
    assert(!eik3_is_point_source(eik, l1));
    if (eik3_is_point_source(eik, l2[i])) {
      do_1pt_update(eik, l, l2[i]);
      goto cleanup;
    }
    utetra_init(utetra, eik->mesh, eik->jet, l, l0, l1, l2[i]);
    if (utetra_is_degenerate(utetra)) {
      continue;
    }
    lam[0] = lam[1] = 0;
    utetra_set_lambda(utetra, lam);
    utetra_solve(utetra);
    utetra_get_lag_mults(utetra, alpha);
    if (dbl3_maxnorm(alpha) <= 1e-15) {
      utetra_get_jet(utetra, &jet);
      if (jet.f < eik->jet[l].f) {
        eik->jet[l] = jet;
      }
    }
  }

cleanup:
  utetra_dealloc(&utetra);
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

  size_t l1;
  if (!do_2pt_updates(eik, l, l0, &l1)) {
    return;
  }

  // TODO: not using l3 right now, but might want to use it later if
  // we want to try out "volume updates"
  size_t *l2, *l3;
  int n = get_l2_and_l3(eik, l0, l1, &l2, &l3);

  do_tetra_updates(eik, l, l0, l1, l2, n);

  free(l2);
  free(l3);
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

static void full_update(eik3_s *eik, size_t l) {
  assert(isinf(eik->jet[l].f));

  int nvv = mesh3_nvv(eik->mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(eik->mesh, l, vv);

  int lvalid_capacity = nvv;
  size_t *lvalid = malloc(lvalid_capacity*sizeof(size_t));

  int num_valid = 0;
  for (int i = 0; i < nvv; ++i) {
    if (eik->state[vv[i]] == VALID) {
      if (num_valid == lvalid_capacity) {
        lvalid_capacity *= 2;
        lvalid = realloc(lvalid, lvalid_capacity*sizeof(size_t));
      }
      lvalid[num_valid++] = vv[i];
    }
  }
  // Now, VALID get neighbors of neighbors
  for (int i = 0; i < nvv; ++i) {
    int nvv_ = mesh3_nvv(eik->mesh, vv[i]);
    size_t *vv_ = malloc(nvv_*sizeof(size_t));
    mesh3_vv(eik->mesh, vv[i], vv_);
    for (int j = 0; j < nvv_; ++j) {
      // Move on if this vertex isn't VALID
      if (eik->state[vv_[j]] != VALID) {
        continue;
      }
      // Check if we've already found this VALID vertex
      bool repeat = false;
      for (int k = 0; k < num_valid; ++k) {
        if (vv_[j] == lvalid[k]) {
          repeat = true;
          break;
        }
      }
      if (repeat) {
        continue;
      }
      // Append this vertex
      if (num_valid == lvalid_capacity) {
        lvalid_capacity *= 2;
        lvalid = realloc(lvalid, lvalid_capacity*sizeof(size_t));
      }
      lvalid[num_valid++] = vv_[j];
    }
    free(vv_);
  }

  assert(num_valid > 3);

  size_t l0, l1, *l2;
  for (int i0 = 0; i0 < num_valid; ++i0) {
    l0 = lvalid[i0];
    for (int i1 = i0 + 1; i1 < num_valid; ++i1) {
      l1 = lvalid[i1];
      l2 = &lvalid[i1 + 1];
      do_tetra_updates(eik, l, l0, l1, l2, num_valid - i1 - 1);
    }
  }

  free(lvalid);
  free(vv);
}

size_t eik3_step(eik3_s *eik) {
  size_t l, l0 = heap_front(eik->heap);
  assert(eik->state[l0] == TRIAL);
  heap_pop(eik->heap);

  if (!isfinite(eik->jet[l0].f)) {
    full_update(eik, l0);
    assert(isfinite(eik->jet[l0].f));
    eik->full_update[l0] = true;
    ++eik->num_full_updates;
  }

  eik->state[l0] = VALID;

  ++eik->num_valid;

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
      adjust(eik, l);
    }
  }

  free(nb);

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
  return eik->state[l] == VALID && isfinite(eik->jet[l].f)
    && isnan(eik->jet[l].fx) && isnan(eik->jet[l].fy) && isnan(eik->jet[l].fz);
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

jet3 *eik3_get_jet_ptr(eik3_s const *eik) {
  return eik->jet;
}

state_e *eik3_get_state_ptr(eik3_s const *eik) {
  return eik->state;
}

int eik3_get_num_full_updates(eik3_s const *eik) {
  return eik->num_full_updates;
}

int *eik3_get_full_update_ptr(eik3_s const *eik) {
  return eik->full_update;
}
