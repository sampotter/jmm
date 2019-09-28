#include "sjs.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "heap.h"

#define NUM_NB 8
#define NUM_CELL_VERTS 4

typedef struct {
  dbl f, fx, fy, fxy;
} jet;

struct sjs {
  ivec2 shape;
  dbl h;
  int nnodes, ncells;
  int nb_ind_offsets[NUM_NB + 1];
  int tri_cell_ind_offsets[NUM_NB];
  int cell_vert_ind_offsets[NUM_CELL_VERTS];
  int nb_cell_ind_offsets[NUM_CELL_VERTS];
  sfield s;
  vfield grad_s;
  bicubic *bicubics;
  jet *jets;
  state *states;
  int *parents;
  int *positions;
  heap_s *heap;
};

void sjs_alloc(sjs_s **sjs) {
  *sjs = malloc(sizeof(sjs_s));
  assert(*sjs != NULL);
}

void sjs_dealloc(sjs_s **sjs) {
  free(*sjs);
  *sjs = NULL;
}

int sjs_lindexi(sjs_s *sjs, ivec2 ind) {
  return (sjs->shape.i + 2)*ind.j + ind.i;
}

int sjs_lindexe(sjs_s *sjs, ivec2 ind) {
  return (sjs->shape.i + 2)*(ind.j + 1) + ind.i + 1;
}

void sjs_vindex(sjs_s *sjs, int l, int *i, int *j) {
  int mpad = sjs->shape.i + 2;
  *i = l/mpad - 1;
  *j = l%mpad - 1;
}

ivec2 offsets[NUM_NB + 1] = {
  {.i = -1, .j = -1},
  {.i = -1, .j =  0},
  {.i = -1, .j =  1},
  {.i =  0, .j =  1},
  {.i =  1, .j =  1},
  {.i =  1, .j =  0},
  {.i =  1, .j = -1},
  {.i =  0, .j = -1},
  {.i = -1, .j = -1}
};

void sjs_set_nb_ind_offsets(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    sjs->nb_ind_offsets[i] = sjs_lindexi(sjs, offsets[i]);
  }
}

ivec2 tri_cell_offsets[NUM_NB] = {
  {.i = -2, .j = -1},
  {.i = -2, .j =  0},
  {.i = -1, .j =  1},
  {.i =  0, .j =  1},
  {.i =  1, .j =  0},
  {.i =  1, .j = -1},
  {.i =  0, .j = -2},
  {.i = -1, .j = -2}
};

void sjs_set_tri_cell_ind_offsets(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB; ++i) {
    sjs->tri_cell_ind_offsets[i] = sjs_lindexi(sjs, tri_cell_offsets[i]);
  }
}

ivec2 cell_vert_offsets[NUM_CELL_VERTS] = {
  {.i = 0, .j = 0},
  {.i = 1, .j = 0},
  {.i = 0, .j = 1},
  {.i = 1, .j = 1}
};

void sjs_set_cell_vert_ind_offsets(sjs_s *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->cell_vert_ind_offsets[i] = sjs_lindexi(sjs, cell_vert_offsets[i]);
  }
}

ivec2 nb_cell_offsets[NUM_CELL_VERTS] = {
  {.i = -1, .j = -1},
  {.i =  0, .j = -1},
  {.i = -1, .j =  0},
  {.i =  0, .j = -1}
};

void sjs_set_nb_cell_ind_offsets(sjs_s *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->nb_cell_ind_offsets[i] = sjs_lindexi(sjs, nb_cell_offsets[i]);
  }
}

void sjs_init(sjs_s *sjs, ivec2 shape, dbl h, sfield s, vfield grad_s) {
  sjs->shape = shape;
  sjs->ncells = (shape.i + 1)*(shape.j + 1);
  sjs->nnodes = (shape.i + 2)*(shape.j + 2);
  sjs->h = h;
  sjs->s = s;
  sjs->bicubics = malloc(sjs->ncells*sizeof(bicubic));
  sjs->jets = malloc(sjs->nnodes*sizeof(jet));
  sjs->states = malloc(sjs->nnodes*sizeof(state));
  sjs->parents = malloc(sjs->nnodes*sizeof(int));
  sjs->positions = malloc(sjs->nnodes*sizeof(int));

  assert(sjs->bicubics != NULL);
  assert(sjs->jets != NULL);
  assert(sjs->states != NULL);
  assert(sjs->parents != NULL);
  assert(sjs->positions != NULL);

#ifndef NDEBUG
  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->positions[l] = NO_INDEX;
  }
#endif

  int capacity = (int) 3*sqrt(sjs->shape.i*sjs->shape.j);
  heap_init(sjs->heap, capacity);

  sjs_set_nb_ind_offsets(sjs);
  sjs_set_tri_cell_ind_offsets(sjs);
  sjs_set_cell_vert_ind_offsets(sjs);
  sjs_set_nb_cell_ind_offsets(sjs);

  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->jets[l].f = INFINITY;
    sjs->states[l] = FAR;
  }

  for (int lc = 0; lc < sjs->ncells; ++lc) {
    sjs->parents[lc] = NO_PARENT;
  }
}

void sjs_deinit(sjs_s *sjs) {
  free(sjs->bicubics);
  free(sjs->jets);
  free(sjs->states);
}

dvec2 sjs_xy(sjs_s *sjs, int l) {
  int mpad = sjs->shape.i + 2;
  dvec2 xy = {
    .x = sjs->h*(l/mpad - 1),
    .y = sjs->h*(l%mpad - 1)
  };
  return xy;
}

dvec2 sjs_cell_center(sjs_s *sjs, int lc) {
  dvec2 xyc = sjs_xy(sjs, lc);
  dbl hhalf = sjs->h/2;
  xyc.x += hhalf;
  xyc.y += hhalf;
  return xyc;
}

int sjs_add_fac_pt_src(sjs_s *sjs, ivec2 ind0, dbl r0) {
  int ncells = (sjs->shape.i + 1)*(sjs->shape.j + 1);

  int nf = 0;
  int l0 = sjs_lindexe(sjs, ind0);
  dvec2 xy0 = sjs_xy(sjs, l0);
  for (int lc = 0; lc < ncells; ++lc) {
    dvec2 xyc = sjs_cell_center(sjs, lc);
    sjs->parents[lc] = dvec2_dist(xyc, xy0) <= r0 ? l0 : -1;

    if (sjs->parents[lc] != NO_PARENT) {
      ++nf;
    }
  }

  jet *J = &sjs->jets[l0];
  J->f = J->fx = J->fy = J->fxy = 0;
  sjs->states[l0] = TRIAL;
  heap_insert(sjs->heap, l0);

  return nf;
}

dbl sjs_get_s(sjs_s *sjs, int l) {
  return sjs->s(sjs_xy(sjs, l));
}

dbl sjs_T(sjs_s *sjs, int l) {
  return sjs->jets[l].f;
}

bicubic_variable tri_bicubic_vars[NUM_NB] = {
  MU, MU, LAMBDA, LAMBDA, MU, MU, LAMBDA, LAMBDA
};

int tri_edges[NUM_NB] = {1, 1, 0, 0, 0, 0, 1, 1};

typedef struct {
  sjs_s *sjs;
  bicubic_variable var;
  cubic cubic;
  dvec2 xy0, xy1;
  int parent;
} F_data;

dbl F(F_data *data, dbl t) {
  dvec2 xyt = dvec2_ccomb(data->xy0, data->xy1, t);
  dbl T = cubic_f(&data->cubic, t);
  dbl s = data->sjs->s(xyt);
  dbl L = sqrt(1 + t*t);
  dbl tmp = T + data->sjs->h*s*L;
  if (data->parent != NO_PARENT) {
    dvec2 xyf = sjs_xy(data->sjs, data->parent);
    tmp += dvec2_dist(xyt, xyf);
  }
  return tmp;
}

dbl dF_dt(F_data *data, dbl t) {
  dvec2 xyt = dvec2_ccomb(data->xy0, data->xy1, t);
  dbl s = data->sjs->s(xyt);
  dvec2 ds = data->sjs->grad_s(xyt);
  dbl ds_dt = data->var == LAMBDA ? ds.x : ds.y;
  dbl dT_dt = cubic_df(&data->cubic, t);
  dbl L = sqrt(1 + t*t);
  dbl dL_dt = t/L;
  dbl tmp = dT_dt + data->sjs->h*(ds_dt*L + s*dL_dt);
  if (data->parent != NO_PARENT) {
    dvec2 xyf = sjs_xy(data->sjs, data->parent);
    dbl rf = dvec2_dist(xyt, xyf);
    tmp += data->var == LAMBDA ?
      data->sjs->h*(xyt.x + data->sjs->h*t - xyt.x)/pow(rf, 3) :
      data->sjs->h*(xyt.y + data->sjs->h*t - xyt.y)/pow(rf, 3);
  }
  return tmp;
}

int sgn(dbl x) {
  if (x > 0) {
    return 1;
  } else if (x < 0) {
    return -1;
  } else {
    return 0;
  }
}

bool sjs_tri(sjs_s *sjs, int l, int l0, int l1, int i0) {
  int lc = l + sjs->tri_cell_ind_offsets[i0];

  bicubic *bicubic = &sjs->bicubics[lc];

  F_data data;
  data.sjs = sjs;
  data.var = tri_bicubic_vars[i0];
  data.cubic = bicubic_restrict(bicubic, data.var, tri_edges[i0]);
  data.xy0 = sjs_xy(sjs, l0);
  data.xy1 = sjs_xy(sjs, l1);
  data.parent = sjs->parents[lc];

  dbl lam, a, b, c, d, fa, fb, fc, fd, dm, df, ds, dd, tmp;

  fa = dF_dt(&data, 0);
  if (fabs(fa) <= EPS) {
    lam = 0;
    goto found;
  }

  fb = dF_dt(&data, 1);
  if (fabs(fb) <= EPS) {
    lam = 1;
    goto found;
  }

  if (sgn(fa) == sgn(fb)) {
    lam = sgn(fa) == 1 ? 0 : 1;
    goto found;
  }

  c = 0;
  fc = fa;
  for (;;) {
    if (fabs(fc) < fabs(fb)) {
      tmp = b; b = c; c = tmp;
      tmp = fb; fb = fc; fc = tmp;
      a = c;
      fa = fc;
    }
    if (fabs(b - c) <= EPS) {
      break;
    }
    dm = (c - b)/2;
    df = fa - fb;
    ds = df == 0 ? dm : -fb*(a - b)/df;
    dd = sgn(ds) != sgn(dm) || fabs(ds) > fabs(dm) ? dm : ds;
    if (fabs(dd) < EPS) {
      dd = EPS*sgn(dm)/2;
    }
    d = b + dd;
    fd = dF_dt(&data, d);
    if (fd == 0) {
      c = d;
      b = c;
      fc = fd;
      fb = fc;
      break;
    }
    a = b;
    b = d;
    fa = fb;
    fb = fd;
    if (sgn(fb) == sgn(fc)) {
      c = a;
      fc = fa;
    }
  }
  lam = (b + c)/2;

  found: {
    dbl T = F(&data, lam);
    jet *J = &sjs->jets[l];
    if (T < J->f) {
      J->f = T;
      dvec2 xy = sjs_xy(sjs, l);
      dvec2 xylam = dvec2_ccomb(data.xy0, data.xy1, lam);
      dbl L = sqrt(1 + lam*lam);
      J->fx = sjs_get_s(sjs, l)*(xy.x - xylam.x)/L;
      J->fy = sjs_get_s(sjs, l)*(xy.y - xylam.y)/L;
      return true;
    } else {
      return false;
    }
  }
}

bool sjs_line(sjs_s *sjs, int l, int l0, int i0) {
  dbl s = sjs_get_s(sjs, l), s0 = sjs_get_s(sjs, l0);
  dbl T0 = sjs_T(sjs, l0);
  dbl T = T0 + sjs->h*(s + s0)/2;
  jet *J = &sjs->jets[l];
  if (T < J->f) {
    J->f = T;
    dbl dist = i0 % 2 == 0 ? SQRT2 : 1;
    J->fx = s*offsets[i0].i/dist;
    J->fy = s*offsets[i0].j/dist;
    return true;
  } else {
    return false;
  }
}

bool sjs_valid_cell(sjs_s *sjs, int lc) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    if (sjs->states[lc + sjs->cell_vert_ind_offsets[i]] != VALID) {
      return false;
    }
  }
  return true;
}

dbl sjs_est_fxy(sjs_s *sjs, int l, int lc) {
  dbl fx[NUM_CELL_VERTS], fy[NUM_CELL_VERTS];

  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc + sjs->cell_vert_ind_offsets[i];
    fx[i] = sjs->jets[l].fx;
    fy[i] = sjs->jets[l].fy;
  }

  dbl fxy[NUM_CELL_VERTS] = {
    (fy[1] - fy[0])/sjs->h, // left
    (fx[3] - fx[1])/sjs->h, // bottom
    (fx[2] - fx[0])/sjs->h, // top
    (fy[3] - fy[2])/sjs->h // right
  };

  static dbl lams[4] = {-1/2, 1/2, 1/2, 3/2};
  static dbl mus[4] = {1/2, -1/2, 3/2, 1/2};

  int i = 0, dl = l - lc;
  while (dl != sjs->cell_vert_ind_offsets[i]) ++i;

  dbl lam = lams[i], mu = mus[i];

  return (1 - mu)*((1 - lam)*fxy[0] + lam*fxy[1]) +
    mu*((1 - lam)*fxy[2] + lam*fxy[3]);
}

dbl sjs_rf(sjs_s *sjs, int l, int lf) {
  dvec2 xy = sjs_xy(sjs, l), xyf = sjs_xy(sjs, lf);
  return dvec2_dist(xy, xyf);
}

dbl sjs_rfx(sjs_s *sjs, int l, int lf) {
  dvec2 xy = sjs_xy(sjs, l), xyf = sjs_xy(sjs, lf);
  return (xy.x - xyf.x)/dvec2_dist(xy, xyf);
}

dbl sjs_rfy(sjs_s *sjs, int l, int lf) {
  dvec2 xy = sjs_xy(sjs, l), xyf = sjs_xy(sjs, lf);
  return (xy.y - xyf.y)/dvec2_dist(xy, xyf);
}

dbl sjs_rfxy(sjs_s *sjs, int l, int lf) {
  dvec2 xy = sjs_xy(sjs, l), xyf = sjs_xy(sjs, lf);
  return -(xy.x - xyf.x)*(xy.y - xyf.y)/pow(dvec2_dist(xy, xyf), 3);
}

void sjs_update_cell(sjs_s *sjs, int lc) {
  int l[4];
  jet *J[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    l[i] = lc + sjs->cell_vert_ind_offsets[i];
    J[i] = &sjs->jets[l[i]];
  }

  dbl data[4][4];

  // TODO: this is a screwy way to do this---we want to lay this out
  // in memory so it will eventually be easy to just do this using
  // SIMD instructions
  //
  // TODO: don't forget to scale by h below...
  int lf = sjs->parents[lc];
  if (lf == NO_PARENT) {
    data[0][0] = J[0]->f;
    data[0][1] = J[2]->f;
    data[0][2] = J[0]->fy;
    data[0][3] = J[2]->fy;
    data[1][0] = J[1]->f;
    data[1][1] = J[3]->f;
    data[1][2] = J[1]->fy;
    data[1][3] = J[3]->fy;
    data[2][0] = J[0]->fx;
    data[2][1] = J[2]->fx;
    data[2][2] = J[0]->fxy;
    data[2][3] = J[2]->fxy;
    data[3][0] = J[1]->fx;
    data[3][1] = J[3]->fx;
    data[3][2] = J[1]->fxy;
    data[3][3] = J[3]->fxy;
  } else {
    data[0][0] = J[0]->f - sjs_rf(sjs, l[0], lf);
    data[0][1] = J[2]->f - sjs_rf(sjs, l[2], lf);
    data[0][2] = J[0]->fy - sjs_rfy(sjs, l[0], lf);
    data[0][3] = J[2]->fy - sjs_rfy(sjs, l[2], lf);
    data[1][0] = J[1]->f - sjs_rf(sjs, l[1], lf);
    data[1][1] = J[3]->f - sjs_rf(sjs, l[3], lf);
    data[1][2] = J[1]->fy - sjs_rfy(sjs, l[1], lf);
    data[1][3] = J[3]->fy - sjs_rfy(sjs, l[3], lf);
    data[2][0] = J[0]->fx - sjs_rfx(sjs, l[0], lf);
    data[2][1] = J[2]->fx - sjs_rfx(sjs, l[2], lf);
    data[2][2] = J[0]->fxy - sjs_rfxy(sjs, l[0], lf);
    data[2][3] = J[2]->fxy - sjs_rfxy(sjs, l[2], lf);
    data[3][0] = J[1]->fx - sjs_rfx(sjs, l[1], lf);
    data[3][1] = J[3]->fx - sjs_rfx(sjs, l[3], lf);
    data[3][2] = J[1]->fxy - sjs_rfxy(sjs, l[1], lf);
    data[3][3] = J[3]->fxy - sjs_rfxy(sjs, l[3], lf);
  }

  bicubic_set_A(&sjs->bicubics[lc], data);
}

void sjs_update_adj_cells(sjs_s *sjs, int l) {
  int lc[NUM_CELL_VERTS];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    lc[i] = l + sjs->nb_cell_ind_offsets[i];
  }

  int nvalid = 0;
  bool valid[NUM_CELL_VERTS];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    valid[i] = sjs_valid_cell(sjs, lc[i]);
    if (valid[i]) {
      ++nvalid;
    }
  }

  dbl fxy_mean = 0;
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    fxy_mean += sjs_est_fxy(sjs, l, lc[i]);
  }
  fxy_mean /= nvalid;

  sjs->jets[l].fxy = fxy_mean;

  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    if (valid[i]) {
      sjs_update_cell(sjs, lc[i]);
    }
  }
}

void sjs_update(sjs_s *sjs, int l) {
  // TODO: need to incorporate factoring...

  bool done[NUM_NB], updated = false;
  memset(done, 0x0, NUM_NB*sizeof(bool));
  for (int i = 1, l0, l1; i < 8; i += 2) {
    l0 = l + sjs->nb_ind_offsets[i];
    if (sjs->states[l0] == VALID) {
      l1 = l + sjs->nb_ind_offsets[i - 1];
      if (sjs->states[l1] == VALID) {
        updated |= sjs_tri(sjs, l, l0, l1, i);
        done[i] = done[i - 1] = true;
      }
      l1 = l + sjs->nb_ind_offsets[i + 1];
      if (sjs->states[l1] == VALID) {
        updated |= sjs_tri(sjs, l, l0, l1, i);
        done[i] = done[(i + 1) % NUM_NB] = true;
      }
    }
  }
  for (int i = 0, l0; i < 8; ++i) {
    l0 = l + sjs->nb_ind_offsets[i];
    if (!done[i] && sjs->states[l0] == VALID) {
      updated |= sjs_line(sjs, l, l0, i);
    }
  }

  if (updated) {
    sjs_update_adj_cells(sjs, l);
  }
}

void sjs_adjust(sjs_s *sjs, int l0) {
  assert(sjs->states[l0] == TRIAL);
  assert(l0 >= 0);
  assert(l0 < sjs->nnodes);

  heap_swim(sjs->heap, sjs->positions[l0]);
}

bool sjs_inbounds(sjs_s *sjs, int l) {
  return 0 <= l && l < sjs->nnodes;
}

void sjs_step(sjs_s *sjs) {
  int l0 = heap_front(sjs->heap);
  heap_pop(sjs->heap);
  assert(sjs->states[l0] == TRIAL);
  sjs->states[l0] = VALID;

  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_ind_offsets[i];
    if (sjs_inbounds(sjs, l) && sjs->states[l] == FAR) {
      sjs->states[l] = TRIAL;
      heap_insert(sjs->heap, l);
    }
  }

  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_ind_offsets[i];
    if (sjs_inbounds(sjs, l) && sjs->states[l] == TRIAL) {
      sjs_update(sjs, l);
      sjs_adjust(sjs, l);
    }
  }
}

void sjs_solve(sjs_s *sjs) {
  while (heap_size(sjs->heap) > 0) {
    sjs_step(sjs);
  }
}

dbl sjs_value(sjs_s *sjs, int l) {
  return sjs->jets[l].f;
}

bicubic *sjs_bicubic(sjs_s *sjs, ivec2 ind) {
  return &sjs->bicubics[sjs_lindexe(sjs, ind)];
}

int sjs_est_heap_cap(sjs_s *sjs) {
  return 3*sqrt(sjs->shape.i*sjs->shape.j);
}
