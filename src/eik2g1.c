#include "eik2g1.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "bb.h"
#include "heap.h"
#include "mat.h"
#include "utri21.h"
#include "vec.h"

struct eik2g1 {
  grid2_s const *grid;
  grid2info_s grid_info;
  jet21t *jet;
  state_e *state;
  size_t *pos;
  heap_s *heap;
  par2_s *par;
};

static dbl value(eik2g1_s const *eik, int l) {
  assert(l >= 0);
  assert(l < (int)grid2_nind(eik->grid));
  return eik->jet[l].f;
}

static void setpos(eik2g1_s const *eik, int l, int pos) {
  eik->pos[l] = pos;
}

void eik2g1_alloc(eik2g1_s **eik) {
  *eik = malloc(sizeof(eik2g1_s));
}

void eik2g1_dealloc(eik2g1_s **eik) {
  free(*eik);
}

void eik2g1_init(eik2g1_s *eik, grid2_s const *grid) {
  eik->grid = grid;

  grid2info_init(&eik->grid_info, grid);

  size_t num_nodes = grid2_nind(grid);

  eik->jet = malloc(num_nodes*sizeof(jet21t));
  for (size_t i = 0; i < num_nodes; ++i)
    eik->jet[i] = jet21t_make_empty();

  eik->state = malloc(num_nodes*sizeof(jet21t));
  for (size_t i = 0; i < num_nodes; ++i)
    eik->state[i] = FAR;

  eik->pos = malloc(num_nodes*sizeof(jet21t));
  for (size_t i = 0; i < num_nodes; ++i)
    eik->pos[i] = NO_INDEX;

  heap_alloc(&eik->heap);
  heap_init(eik->heap,3*sqrt(num_nodes),(value_f)value,(setpos_f)setpos,eik);

  eik->par = malloc(num_nodes*sizeof(par2_s));
  for (size_t i = 0; i < num_nodes; ++i)
    par2_init_empty(&eik->par[i]);
}

void eik2g1_deinit(eik2g1_s *eik) {
  eik->grid = NULL;

  free(eik->jet);
  eik->jet = NULL;

  free(eik->state);
  eik->state = NULL;

  free(eik->pos);
  eik->pos = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);

  free(eik->par);
  eik->par = NULL;
}

size_t eik2g1_peek(eik2g1_s const *eik) {
  return heap_front(eik->heap);
}

static void tri(eik2g1_s *eik, size_t l, size_t l0, size_t l1) {
  dbl2 xhat, x[2];
  grid2_l2xy(eik->grid, l, xhat);
  grid2_l2xy(eik->grid, l0, x[0]);
  grid2_l2xy(eik->grid, l1, x[1]);

  jet21t jet[2] = {eik->jet[l0], eik->jet[l1]};

  utri21_s utri;
  utri21_init(&utri, xhat, x, jet);

  dbl lam;
  if (!utri21_solve(&utri, &lam))
    return;

  // update jet
  eik->jet[l] = utri.jet;

  // set parent
  eik->par[l] = (par2_s) {.l = {l0, l1}, .b = {1 - lam, lam}};
}

static void update(eik2g1_s *eik, int l) {
  grid2_s const *grid = eik->grid;

  bool inbounds[GRID2_NUM_NB + 1];
  grid2_get_inbounds(grid, &eik->grid_info, l, inbounds);

  for (int i0 = 1; i0 < 8; i0 += 2) {
    if (!inbounds[i0])
      continue;

    size_t l0 = l + eik->grid_info.nb_dl[i0];
    if (eik->state[l0] != VALID)
      continue;

    if (inbounds[i0 - 1]) {
      size_t l1 = l + eik->grid_info.nb_dl[i0 - 1];
      if (eik->state[l1] == VALID)
        tri(eik, l, l0, l1);
    }

    if (inbounds[i0 + 1]) {
      size_t l1 = l + eik->grid_info.nb_dl[i0 + 1];
      if (eik->state[l1] == VALID)
        tri(eik, l, l0, l1);
    }
  }
}

static void adjust(eik2g1_s *eik, size_t l0) {
  assert(eik->state[l0] == TRIAL);
  assert(l0 < grid2_nind(eik->grid));

  heap_swim(eik->heap, eik->pos[l0]);
}

size_t eik2g1_step(eik2g1_s *eik) {
  size_t l0 = heap_front(eik->heap);
  assert(eik->state[l0] == TRIAL);
  heap_pop(eik->heap);
  eik->state[l0] = VALID;

  bool inbounds[GRID2_NUM_NB + 1];
  grid2_get_inbounds(eik->grid, &eik->grid_info, l0, inbounds);

  for (size_t i = 0; i < GRID2_NUM_NB; ++i) {
    if (inbounds[i]) {
      size_t l = l0 + eik->grid_info.nb_dl[i];
      if (eik->state[l] == FAR) {
        eik->state[l] = TRIAL;
        heap_insert(eik->heap, l);
      }
    }
  }

  for (size_t i = 0; i < GRID2_NUM_NB; ++i) {
    if (inbounds[i]) {
      size_t l = l0 + eik->grid_info.nb_dl[i];
      if (eik->state[l] == TRIAL) {
        update(eik, l);
        adjust(eik, l);
      }
    }
  }

  return l0;
}

void eik2g1_solve(eik2g1_s *eik) {
  while (heap_size(eik->heap) > 0)
    eik2g1_step(eik);
}

void eik2g1_add_trial(eik2g1_s *eik, int2 const ind, jet21t jet) {
  size_t l = grid2_ind2l(eik->grid, ind);
  eik->jet[l] = jet;
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = TRIAL;
  heap_insert(eik->heap, l);
}

void eik2g1_add_valid(eik2g1_s *eik, int2 const ind, jet21t jet) {
  size_t l = grid2_ind2l(eik->grid, ind);
  eik->jet[l] = jet;
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = VALID;
}

bool eik2g1_is_valid(eik2g1_s const *eik, int2 const ind) {
  size_t l = grid2_ind2l(eik->grid, ind);
  return eik->state[l] == VALID;
}

bool eik2g1_is_trial(eik2g1_s const *eik, int2 const ind) {
  size_t l = grid2_ind2l(eik->grid, ind);
  return eik->state[l] == TRIAL;
}

bool eik2g1_is_far(eik2g1_s const *eik, int2 const ind) {
  size_t l = grid2_ind2l(eik->grid, ind);
  return eik->state[l] == FAR;
}

state_e const *eik2g1_get_state_ptr(eik2g1_s const *eik) {
  return eik->state;
}

jet21t const *eik2g1_get_jet_ptr(eik2g1_s const *eik) {
  return eik->jet;
}

par2_s eik2g1_get_par(eik2g1_s const *eik, int2 const ind) {
  size_t l = grid2_ind2l(eik->grid, ind);
  return eik->par[l];
}

par2_s const *eik2g1_get_par_ptr(eik2g1_s const *eik) {
  return eik->par;
}
