#include <jmm/eik2g1.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <jmm/bb.h>
#include <jmm/heap.h>
#include <jmm/mat.h>
#include <jmm/utri21.h>
#include <jmm/vec.h>

#include "util.h"

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

bool eik2g1_has_par(eik2g1_s const *eik, int2 const ind) {
  size_t l = grid2_ind2l(eik->grid, ind);
  return !par2_is_empty(&eik->par[l]);
}

par2_s eik2g1_get_par(eik2g1_s const *eik, int2 const ind) {
  size_t l = grid2_ind2l(eik->grid, ind);
  return eik->par[l];
}

par2_s const *eik2g1_get_par_ptr(eik2g1_s const *eik) {
  return eik->par;
}

eik2g1_sol_info_s eik2g1_get_sol_info(eik2g1_s const *eik, int2 const ind) {
  eik2g1_sol_info_s sol_info;

  size_t l = grid2_ind2l(eik->grid, ind);
  size_t l0 = eik->par[l].l[0], l1 = eik->par[l].l[1];

  sol_info.lam_T = eik->par[l].b[1];

  dbl2 xhat, x[2], dx;
  grid2_l2xy(eik->grid, l, xhat);
  grid2_l2xy(eik->grid, l0, x[0]);
  grid2_l2xy(eik->grid, l1, x[1]);
  dbl2_sub(x[1], x[0], dx);

  jet22t jet_tau[2];

  jet_tau[0].f = dbl2_norm(x[0]);
  jet_tau[1].f = dbl2_norm(x[1]);

  dbl2_normalized(x[0], jet_tau[0].Df);
  dbl2_normalized(x[1], jet_tau[1].Df);

  dbl22 eye;
  dbl22_eye(eye);

  for (size_t i = 0; i < 2; ++i) {
    dbl22 DT_otimes_DT;
    dbl2_outer(jet_tau[i].Df, jet_tau[i].Df, DT_otimes_DT);

    dbl22_sub(eye, DT_otimes_DT, jet_tau[i].D2f);

    dbl22_dbl_div_inplace(jet_tau[i].D2f, jet_tau[i].f);
  }

  utri21_s utri_tau;
  utri21_init(&utri_tau, xhat, x, jet_tau);
  utri21_solve(&utri_tau, &sol_info.lam_tau);

  dbl22 Dtau_otimes_Dtau;
  dbl2_outer(xhat, xhat, Dtau_otimes_Dtau);
  dbl22_dbl_div_inplace(Dtau_otimes_Dtau, dbl2_dot(xhat, xhat));

  dbl22 D2tau;
  dbl22_sub(eye, Dtau_otimes_Dtau, D2tau);
  dbl22_dbl_div_inplace(D2tau, dbl2_norm(xhat));

  dbl numer = dbl2_wnormsq(D2tau, x[0]);
  dbl denom = dbl2_wnormsq(D2tau, dx);
  sol_info.lam_star = sqrt(clamp(numer/denom, 0, 1));

  dbl T[2] = {eik->jet[l0].f, eik->jet[l1].f};
  dbl DT[2] = {dbl2_dot(dx, eik->jet[l0].Df), dbl2_dot(dx, eik->jet[l1].Df)};
  bb31 bb_T;
  bb31_init_from_1d_data(&bb_T, T, DT, (dbl2) {0, 1});

  sol_info.E0 = eik->jet[l].f - dbl2_norm(xhat);

  dbl2 xlamT;
  dbl2_saxpy(sol_info.lam_T, dx, x[0], xlamT);

  sol_info.That = bb31_f(&bb_T, (dbl2) {1 - sol_info.lam_T, sol_info.lam_T}) +
    dbl2_dist(xhat, xlamT);

  jet22t jet_T[2] = {eik->jet[l0], eik->jet[l1]};

  utri21_s utri_T;
  utri21_init(&utri_T, xhat, x, jet_T);
  utri21_solve(&utri_T, &sol_info.lam_T_check);

  sol_info.FT_lamT = utri21_F(&utri_T, sol_info.lam_T);
  sol_info.Ftau_lamT = utri21_F(&utri_tau, sol_info.lam_T);
  sol_info.Ftau_lamtau = utri21_F(&utri_tau, sol_info.lam_tau);

  sol_info.E0_check0 = sol_info.That - dbl2_norm(xhat);
  sol_info.E0_check1 = sol_info.FT_lamT - dbl2_norm(xhat);
  sol_info.E0_check2 = sol_info.FT_lamT - sol_info.Ftau_lamT
    + sol_info.Ftau_lamT - dbl2_norm(xhat);
  sol_info.E0_check3 = sol_info.FT_lamT - sol_info.Ftau_lamT
    + sol_info.Ftau_lamT - sol_info.Ftau_lamtau
    + sol_info.Ftau_lamtau - dbl2_norm(xhat);

  sol_info.E0_term1 = sol_info.FT_lamT - sol_info.Ftau_lamT;
  sol_info.E0_term2 = sol_info.Ftau_lamT - sol_info.Ftau_lamtau;
  sol_info.E0_term3 = sol_info.Ftau_lamtau - dbl2_norm(xhat);

  jet22t jet_E[2];
  jet22t_sub(&jet_T[0], &jet_tau[0], &jet_E[0]);
  jet22t_sub(&jet_T[1], &jet_tau[1], &jet_E[1]);

  bb31 bb_E;
  bb31_init_from_jet21t(&bb_E, jet_E, x);
  sol_info.E0_term1_v2 = bb31_f(&bb_E, (dbl2) {1 - sol_info.lam_T, sol_info.lam_T});

  sol_info.T0_error = jet_E[0].f;
  sol_info.T1_error = jet_E[1].f;
  sol_info.DT0_error = dbl2_dot(dx, jet_E[0].Df);
  sol_info.DT1_error = dbl2_dot(dx, jet_E[1].Df);

  return sol_info;
}
