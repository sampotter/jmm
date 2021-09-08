#include "eik2g1.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "bb.h"
#include "def.h"
#include "heap.h"
#include "hybrid.h"
#include "mat.h"
#include "vec.h"

struct eik2g1 {
  grid2_s const *grid;
  grid2info_s grid_info;
  jet22t *jet;
  state_e *state;
  size_t *pos;
  heap_s *heap;
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

  eik->jet = malloc(num_nodes*sizeof(jet22t));
  for (size_t i = 0; i < num_nodes; ++i)
    eik->jet[i] = jet22t_make_empty();

  eik->state = malloc(num_nodes*sizeof(jet22t));
  for (size_t i = 0; i < num_nodes; ++i)
    eik->state[i] = FAR;

  eik->pos = malloc(num_nodes*sizeof(jet22t));
  for (size_t i = 0; i < num_nodes; ++i)
    eik->pos[i] = NO_INDEX;

  heap_alloc(&eik->heap);
  heap_init(eik->heap,3*sqrt(num_nodes),(value_f)value,(setpos_f)setpos,eik);
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
}

size_t eik2g1_peek(eik2g1_s const *eik) {
  return heap_front(eik->heap);
}

typedef struct tri_wkspc {
  dbl2 xhat, x[2], dx;
  bb31 T;
} tri_wkspc_s;

dbl tri_F(dbl lam, tri_wkspc_s const *wkspc) {
  dbl T = bb31_f(&wkspc->T, (dbl2) {1 - lam, lam});

  dbl2 xlam;
  dbl2_saxpy(lam, wkspc->dx, wkspc->x[0], xlam);

  dbl2 xhat_minus_xlam;
  dbl2_sub(wkspc->xhat, xlam, xhat_minus_xlam);

  dbl L = dbl2_norm(xhat_minus_xlam);

  return T + L;
}

dbl tri_F_lam(dbl lam, tri_wkspc_s const *wkspc) {
  dbl T_lam = bb31_df(&wkspc->T, (dbl2) {1 - lam, lam}, (dbl2) {1, -1});

  dbl2 xlam;
  dbl2_saxpy(lam, wkspc->dx, wkspc->x[0], xlam);

  dbl2 xhat_minus_xlam;
  dbl2_sub(wkspc->xhat, xlam, xhat_minus_xlam);

  dbl L = dbl2_norm(xhat_minus_xlam);
  dbl L_lam = -dbl2_dot(wkspc->dx, xhat_minus_xlam)/L;

  return T_lam + L_lam;
}

static void tri(eik2g1_s *eik, size_t l, size_t l0, size_t l1) {
  grid2_s const *grid = eik->grid;

  tri_wkspc_s wkspc;

  grid2_l2xy(grid, l, wkspc.xhat);
  grid2_l2xy(grid, l0, wkspc.x[0]);
  grid2_l2xy(grid, l1, wkspc.x[1]);

  dbl2_sub(wkspc.x[1], wkspc.x[0], wkspc.dx);

  bb31_init_from_jet22t(
    &wkspc.T, (jet22t[2]) {eik->jet[l0], eik->jet[l1]}, wkspc.x);

  dbl lam;
  bool found = hybrid((hybrid_cost_func_t)tri_F_lam, 0, 1, &wkspc, &lam);

  dbl T;
  if (found) {
    T = tri_F(lam, &wkspc);
  } else {
    dbl T0 = tri_F(0, &wkspc);
    dbl T1 = tri_F(1, &wkspc);
    T = fmin(T0, T1);
    lam = T0 < T1 ? 0 : 1;
  }

  assert(T > eik->jet[l0].f);
  assert(T > eik->jet[l1].f);

  jet22t *jet = &eik->jet[l];

  if (T >= jet->f)
    return;

  dbl2 xlam, DT;
  dbl2_saxpy(lam, wkspc.dx, wkspc.x[0], xlam);
  dbl2_sub(wkspc.xhat, xlam, DT);
  dbl L = dbl2_norm(DT);
  dbl2_dbl_div_inplace(DT, L);

  dbl22 eye, tt, D2T;
  dbl22_eye(eye);
  dbl2_outer(DT, DT, tt);
  dbl22_sub(eye, tt, D2T);

  jet->f = T;
  jet->fx = DT[0];
  jet->fy = DT[1];
  jet->fxx = D2T[0][0];
  jet->fxy = D2T[0][1];
  jet->fyy = D2T[1][1];
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
  assert(l0 >= 0);
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

void eik2g1_add_trial(eik2g1_s *eik, int2 const ind, jet22t jet) {
  size_t l = grid2_ind2l(eik->grid, ind);
  eik->jet[l] = jet;
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = TRIAL;
  heap_insert(eik->heap, l);
}

void eik2g1_add_valid(eik2g1_s *eik, int2 const ind, jet22t jet) {
  size_t l = grid2_ind2l(eik->grid, ind);
  eik->jet[l] = jet;
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = VALID;
}

state_e const *eik2g1_get_state_ptr(eik2g1_s const *eik) {
  return eik->state;
}

jet22t const *eik2g1_get_jet_ptr(eik2g1_s const *eik) {
  return eik->jet;
}
