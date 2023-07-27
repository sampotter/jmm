#include <jmm/eik2mp.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <jmm/heap.h>
#include <jmm/mat.h>
#include <jmm/utri21.h>

#include "log.h"

struct eik2mp {
  mesh22_s const *mesh;
  jet22t *jet; // eik_vals, grads from eik_grid plus the hessian
  state_e *state; // current_states from eik_grid but actually using FAR, VALID, TRIAL instead of 0, 1, 2
  size_t *pos; // queue_index in the Priority_queue struct in  priority+queue.c
  heap_s *heap; // subset of Priority_queue struct in priority_queue.c
  size_t nvalid; // nothing like this
};

static dbl value(eik2mp_s const *eik, int l) {
  // returns the eikonal value stored at the l-th node
  assert(l >= 0);
  assert(l < (int)mesh22_nverts(eik->mesh));
  return eik->jet[l].f;
}

static void setpos(eik2mp_s const *eik, int l, int pos) {
  // changes the index in the priority_queue (all these methods can be
  // found in priority_queue.c 
  eik->pos[l] = pos;
}

void eik2mp_alloc(eik2mp_s **eik) {
  *eik = malloc(sizeof(eik2mp_s));
}

void eik2mp_dealloc(eik2mp_s **eik) {
  free(*eik);
}

void eik2mp_init(eik2mp_s *eik, mesh22_s const *mesh) {
  eik->mesh = mesh;

  size_t nverts = mesh22_nverts(eik->mesh);

  eik->jet = malloc(nverts*sizeof(jet22t));
  for (size_t l = 0; l < nverts; ++l)
    eik->jet[l] = jet22t_make_empty();

  eik->state = malloc(nverts*sizeof(jet22t));
  for (size_t l = 0; l < nverts; ++l)
    eik->state[l] = FAR;

  eik->pos = malloc(nverts*sizeof(jet22t));
  for (size_t l = 0; l < nverts; ++l)
    eik->pos[l] = NO_INDEX;

  heap_alloc(&eik->heap);
  heap_init(eik->heap, 3*sqrt(nverts), (value_f)value, (setpos_f)setpos, eik);

  eik->nvalid = 0;
}

void eik2mp_deinit(eik2mp_s *eik) {
  eik->mesh = NULL;

  free(eik->jet);
  eik->jet = NULL;

  free(eik->state);
  eik->state = NULL;

  free(eik->pos);
  eik->pos = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);
}

size_t eik2mp_peek(eik2mp_s const *eik) {
  // returns the value in the 0th place in the heap 
  return heap_front(eik->heap);
}

static void tri(eik2mp_s *eik, size_t l, size_t l0, size_t l1) {
  mesh22_s const *mesh = eik->mesh;

  dbl2 xhat, x[2];
  mesh22_get_vert(mesh, l, xhat);
  mesh22_get_vert(mesh, l0, x[0]);
  mesh22_get_vert(mesh, l1, x[1]);

  jet22t jet[2] = {eik->jet[l0], eik->jet[l1]};

  utri21_s utri;
  utri21_init(&utri, xhat, x, jet);

  dbl lam;
  if (!utri21_solve(&utri, &lam))
    return;

  // check if tangent vector of ray lies in cone spanned by DT0 and DT1
  // MAYBE DO THIS AFTER CONSIDERING DT0 AND DT1 USING SNELLS???

  dbl22 DT0_and_DT1;
  dbl2_copy(eik->jet[l0].Df, DT0_and_DT1[0]);
  dbl2_copy(eik->jet[l1].Df, DT0_and_DT1[1]);
  dbl22_transpose(DT0_and_DT1);

  dbl2 cone_coefs;
  dbl22_dbl2_solve(DT0_and_DT1, utri.jet.Df, cone_coefs);
  if (cone_coefs[0] < 0 || cone_coefs[1] < 0)
    return;

  // update jet
  eik->jet[l] = utri.jet;
}

static void get_op_edge(size_t const lf[3], size_t l, size_t le[2]) {
  size_t j = 0;
  for (size_t i = 0; i < 3; ++i)
    if (lf[i] != l)
      le[j++] = lf[i];
}

// static void update(eik2mp_s *eik, size_t l) {
//   log_debug("update(l = %lu)", l);

//   size_t nvf = mesh22_nvf(eik->mesh, l);
//   size_t *vf = malloc(nvf*sizeof(size_t));
//   mesh22_vf(eik->mesh, l, vf);

//   size_t lv[3], le[2];

//   for (size_t i = 0; i < nvf; ++i) {
//     mesh22_fv(eik->mesh, vf[i], lv);
//     get_op_edge(lv, l, le);
//     log_debug("l = %lu, le[0] = %lu, le[1] = %lu", l, le[0], le[1]);
//     if (eik->state[le[0]] == VALID && eik->state[le[1]] == VALID)
//       tri(eik, l, le[0], le[1]);
//   }

//   free(vf);
// }

static void adjust(eik2mp_s *eik, size_t l0) {
  // adjusts the heap once a better update is found
  assert(eik->state[l0] == TRIAL); // do this ONLY for trial nodes
  assert(l0 < mesh22_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l0]); // sends up the heap
}

size_t eik2mp_step(eik2mp_s *eik) {
  size_t l0 = heap_front(eik->heap);
  assert(eik->state[l0] == TRIAL);
  heap_pop(eik->heap); // remove x0 from the heap
  eik->state[l0] = VALID; // set its state to valid

  size_t nvv = mesh22_nvv(eik->mesh, l0); // number of neighbors
  size_t* vv = malloc(nvv*sizeof(size_t)); // indices of the neighbors
  mesh22_vv(eik->mesh, l0, vv);

  for (size_t i = 0; i < nvv; ++i) {
    size_t l = vv[i]; // neighbor of x0
    if (eik->state[l] == FAR) { 
      eik->state[l] = TRIAL;
      heap_insert(eik->heap, l);
    }
  }

  // Return early if all of l0's neighbors are VALID:

  size_t ntrial = 0;
  for (size_t i = 0; i < nvv; ++i)
    ntrial += eik->state[vv[i]] == TRIAL;
  if (ntrial == 0)
    goto cleanup;

  // First, find the edges on the VALID front that are incident on l0:

  size_t nvf = mesh22_nvf(eik->mesh, l0);
  size_t *vf = malloc(nvf*sizeof(size_t));
  mesh22_vf(eik->mesh, l0, vf);

  size_t lv[3], le[2], l1[2] = {NO_INDEX, NO_INDEX}, nvalid;
  for (size_t i = 0; i < nvf; ++i) {
    if (l1[0] != (size_t)NO_INDEX && l1[1] != (size_t)NO_INDEX)
      break;

    mesh22_fv(eik->mesh, vf[i], lv);
    get_op_edge(lv, l0, le);

    nvalid = (eik->state[le[0]] == VALID) + (eik->state[le[1]] == VALID);

    if (nvalid == 0 || nvalid == 2)
      continue;

    size_t l_valid = eik->state[le[0]] == VALID ? le[0] : le[1];

    if (l1[0] == (size_t)NO_INDEX)
      l1[0] = l_valid;
    if (l1[0] != l_valid && l1[1] == (size_t)NO_INDEX)
      l1[1] = l_valid;
  }

  free(vf);

  assert(l1[0] != (size_t)NO_INDEX || l1[1] != (size_t)NO_INDEX);

  // Update each TRIAL node l neighboring l0 from the edges on the
  // VALID front:

  for (size_t i = 0; i < nvv; ++i) {
    size_t l = vv[i];

    if (eik->state[l] != TRIAL)
      continue;

    if (l1[0] != (size_t)NO_INDEX)
      tri(eik, l, l0, l1[0]);

    if (l1[1] != (size_t)NO_INDEX)
      tri(eik, l, l0, l1[1]);

    adjust(eik, l);
  }

cleanup:
  free(vv);

  return l0;
}

void eik2mp_solve(eik2mp_s *eik) {
  while (heap_size(eik->heap) > 0) {
    // size_t l0 = eik2mp_step(eik);
    eik2mp_step(eik);
    ++eik->nvalid;
  }
}

void eik2mp_add_trial(eik2mp_s *eik, size_t l, jet22t jet) {
  eik->jet[l] = jet;
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = TRIAL;
  heap_insert(eik->heap, l);
}

void eik2mp_add_valid(eik2mp_s *eik, size_t l, jet22t jet) {
  eik->jet[l] = jet;
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = VALID;
  ++eik->nvalid;
}

bool eik2mp_is_valid(eik2mp_s const *eik, size_t l) {
  return eik->state[l] == VALID;
}

bool eik2mp_is_trial(eik2mp_s const *eik, size_t l) {
  return eik->state[l] == TRIAL;
}

bool eik2mp_is_far(eik2mp_s const *eik, size_t l) {
  return eik->state[l] == FAR;
}

state_e const *eik2mp_get_state_ptr(eik2mp_s const *eik) {
  return eik->state;
}

jet22t const *eik2mp_get_jet_ptr(eik2mp_s const *eik) {
  return eik->jet;
}
