#include "sjs.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "heap.h"

#define NUM_CELL_VERTS 4
#define NUM_NB 8
#define NUM_NEARBY_CELLS 16

typedef struct {
  dbl f, fx, fy, fxy;
} jet;

struct sjs {
  ivec2 shape;
  dvec2 xymin;
  dbl h;
  int nnodes, ncells;
  int nb_dl[NUM_NB + 1];
  int vert_dl[NUM_CELL_VERTS];
  int tri_dlc[NUM_NB];
  int nb_dlc[NUM_CELL_VERTS];
  int nearby_dlc[NUM_NEARBY_CELLS];
  sfield s;
  vfield grad_s;
  bicubic *bicubics;
  jet *jets;
  state_e *states;
  int *cell_parents;
  int *node_parents;
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

int sjs_clindexi(sjs_s *sjs, ivec2 ind) {
  return (sjs->shape.i + 1)*ind.j + ind.i;
}

int sjs_clindexe(sjs_s *sjs, ivec2 cind) {
  return (sjs->shape.i + 1)*(cind.j + 1) + cind.i + 1;
}

int sjs_l2lc(sjs_s *sjs, int l) {
  return l - l/(sjs->shape.i + 2);
}

int sjs_lc2l(sjs_s *sjs, int lc) {
  return lc + lc/(sjs->shape.i + 1);
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

void sjs_set_nb_dl(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    sjs->nb_dl[i] = sjs_lindexi(sjs, offsets[i]);
  }
}

ivec2 cell_vert_offsets[NUM_CELL_VERTS] = {
  {.i = 0, .j = 0},
  {.i = 1, .j = 0},
  {.i = 0, .j = 1},
  {.i = 1, .j = 1}
};

void sjs_set_vert_dl(sjs_s *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->vert_dl[i] = sjs_lindexi(sjs, cell_vert_offsets[i]);
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

void sjs_set_tri_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB; ++i) {
    sjs->tri_dlc[i] = sjs_clindexi(sjs, tri_cell_offsets[i]);
  }
}

ivec2 nb_cell_offsets[NUM_CELL_VERTS] = {
  {.i = -1, .j = -1},
  {.i =  0, .j = -1},
  {.i = -1, .j =  0},
  {.i =  0, .j =  0}
};

void sjs_set_nb_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->nb_dlc[i] = sjs_clindexi(sjs, nb_cell_offsets[i]);
  }
}

ivec2 nearby_cell_offsets[NUM_NEARBY_CELLS] = {
  {.i = -2, .j = -2}, {.i = -2, .j = -1}, {.i = -2, .j = 0}, {.i = -2, .j = 1},
  {.i = -1, .j = -2}, {.i = -1, .j = -1}, {.i = -1, .j = 0}, {.i = -1, .j = 1},
  {.i =  0, .j = -2}, {.i =  0, .j = -1}, {.i =  0, .j = 0}, {.i =  0, .j = 1},
  {.i =  1, .j = -2}, {.i =  1, .j = -1}, {.i =  1, .j = 0}, {.i =  1, .j = 1},
};

void sjs_set_nearby_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_NEARBY_CELLS; ++i) {
    sjs->nearby_dlc[i] = sjs_clindexi(sjs, nearby_cell_offsets[i]);
  }
}

dvec2 sjs_xy(sjs_s *sjs, int l) {
  int mpad = sjs->shape.i + 2;
  dvec2 xy = {
    .x = sjs->xymin.x + sjs->h*(l%mpad - 1),
    .y = sjs->xymin.y + sjs->h*(l/mpad - 1)
  };
  return xy;
}

// TODO: since the margins are BOUNDARY nodes, we actually don't need
// to allocate an extra margin of cells, since they will never be
// initialized (i.e., they will never have all of their vertex nodes
// become VALID because of the margin...)
void sjs_init(sjs_s *sjs, ivec2 shape, dvec2 xymin, dbl h, sfield s,
              vfield grad_s) {
  sjs->shape = shape;
  sjs->ncells = (shape.i + 1)*(shape.j + 1);
  sjs->nnodes = (shape.i + 2)*(shape.j + 2);
  sjs->xymin = xymin;
  sjs->h = h;
  sjs->s = s;
  sjs->grad_s = grad_s;
  sjs->bicubics = malloc(sjs->ncells*sizeof(bicubic));
  sjs->jets = malloc(sjs->nnodes*sizeof(jet));
  sjs->states = malloc(sjs->nnodes*sizeof(state_e));
  sjs->cell_parents = malloc(sjs->ncells*sizeof(int));
  sjs->node_parents = malloc(sjs->nnodes*sizeof(int));
  sjs->positions = malloc(sjs->nnodes*sizeof(int));

  assert(sjs->bicubics != NULL);
  assert(sjs->jets != NULL);
  assert(sjs->states != NULL);
  assert(sjs->cell_parents != NULL);
  assert(sjs->positions != NULL);

#ifndef NDEBUG
  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->positions[l] = NO_INDEX;
  }
#endif

  heap_alloc(&sjs->heap);

  value_b value = ^(int l) {
    dbl T = sjs->jets[l].f;
    int lf = sjs->node_parents[l];
    if (lf != NO_PARENT) {
      dvec2 xy = sjs_xy(sjs, l), xyf = sjs_xy(sjs, lf);
      dbl rf = dvec2_dist(xy, xyf);
      T += rf;
    }
    return T;
  };

  setpos_b setpos = ^(int l, int pos) {
    sjs->positions[l] = pos;
  };

  int capacity = (int) 3*sqrt(sjs->shape.i*sjs->shape.j);
  heap_init(sjs->heap, capacity, value, setpos);

  sjs_set_nb_dl(sjs);
  sjs_set_vert_dl(sjs);
  sjs_set_tri_dlc(sjs);
  sjs_set_nb_dlc(sjs);
  sjs_set_nearby_dlc(sjs);

  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->jets[l].f = INFINITY;
    sjs->states[l] = FAR;
  }

  for (int i = 0, l; i < shape.i + 2; ++i) {
    l = i;
    sjs->states[l] = BOUNDARY;
  }
  for (int j = 1, l; j < shape.i + 1; ++j) {
    l = (shape.i + 2)*j;
    sjs->states[l] = BOUNDARY;
    l += shape.i + 1;
    sjs->states[l] = BOUNDARY;
  }
  for (int i = 0, l; i < shape.i + 2; ++i) {
    l = (shape.j + 1)*(shape.i + 2) + i;
    sjs->states[l] = BOUNDARY;
  }

  for (int lc = 0; lc < sjs->ncells; ++lc) {
    sjs->cell_parents[lc] = NO_PARENT;
  }

  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->node_parents[l] = NO_PARENT;
  }
}

void sjs_deinit(sjs_s *sjs) {
  free(sjs->bicubics);
  free(sjs->jets);
  free(sjs->states);
  free(sjs->cell_parents);
  free(sjs->node_parents);
  free(sjs->positions);

  sjs->bicubics = NULL;
  sjs->jets = NULL;
  sjs->states = NULL;
  sjs->cell_parents = NULL;
  sjs->node_parents = NULL;
  sjs->positions = NULL;

  heap_deinit(sjs->heap);
  heap_dealloc(&sjs->heap);
}

dvec2 sjs_cell_center(sjs_s *sjs, int lc) {
  int l = sjs_lc2l(sjs, lc);
  dvec2 xyc = sjs_xy(sjs, l);
  dbl hhalf = sjs->h/2;
  xyc.x += hhalf;
  xyc.y += hhalf;
  return xyc;
}

void sjs_add_fac_pt_src(sjs_s *sjs, ivec2 indf, dbl rmax, int *nf, int *nfc) {
  if (nf != NULL) {
    *nf = 0;
  }
  if (nfc != NULL) {
    *nfc = 0;
  }

  int lf = sjs_lindexe(sjs, indf);
  dvec2 xyf = sjs_xy(sjs, lf);
  dbl r;
  for (int l = 0; l < sjs->nnodes; ++l) {
    r = dvec2_dist(sjs_xy(sjs, l), xyf);
    sjs->node_parents[l] = r <= rmax + EPS ? lf : NO_PARENT;
    if (nf != NULL && sjs->node_parents[l] != NO_PARENT) {
      ++*nf;
    }
  }
  for (int lc = 0; lc < sjs->ncells; ++lc) {
    dvec2 xyc = sjs_cell_center(sjs, lc);
    r = dvec2_dist(xyc, xyf);
    sjs->cell_parents[lc] = r <= rmax + EPS ? lf : NO_PARENT;
    if (nfc != NULL && sjs->cell_parents[lc] != NO_PARENT) {
      ++*nfc;
    }
  }

  jet *J = &sjs->jets[lf];
  J->f = J->fx = J->fy = J->fxy = 0;
  sjs->states[lf] = TRIAL;
  heap_insert(sjs->heap, lf);
}

dbl sjs_get_s(sjs_s *sjs, int l) {
  return sjs->s(sjs_xy(sjs, l));
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

static bool tri(sjs_s *sjs, int l, int l0, int l1, int i0) {
  assert(i0 >= 0);
  assert(i0 < NUM_NB);

  int lc = sjs_l2lc(sjs, l) + sjs->tri_dlc[i0];
  assert(lc >= 0);
  assert(lc < sjs->ncells);

  bicubic *bicubic = &sjs->bicubics[lc];

  F_data data;
  data.sjs = sjs;
  data.var = tri_bicubic_vars[i0];
  data.cubic = bicubic_restrict(bicubic, data.var, tri_edges[i0]);
  data.xy0 = sjs_xy(sjs, l0);
  data.xy1 = sjs_xy(sjs, l1);
  data.parent = sjs->cell_parents[lc];

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

  bool updated = false;
  found: {
    dbl T = F(&data, lam);
    jet *J = &sjs->jets[l];
    if (T < J->f) {
      updated = true;

      J->f = T;
      dvec2 xy = sjs_xy(sjs, l);
      dvec2 xylam = dvec2_ccomb(data.xy0, data.xy1, lam);
      dbl L = sqrt(1 + lam*lam);
      J->fx = sjs_get_s(sjs, l)*(xy.x - xylam.x)/L;
      J->fy = sjs_get_s(sjs, l)*(xy.y - xylam.y)/L;

      int lf = sjs->node_parents[l];
      if (lf != NO_PARENT) {
        dvec2 xyf = sjs_xy(sjs, lf);
        dbl rf = dvec2_dist(xy, xyf);
        J->f -= rf;
        J->fx -= (xy.x - xyf.x)/rf;
        J->fy -= (xy.y - xyf.y)/rf;
      }
    }
  }
  return updated;
}

static bool line(sjs_s *sjs, int l, int l0, int i0) {
  dbl s = sjs_get_s(sjs, l), s0 = sjs_get_s(sjs, l0);
  dbl T0 = sjs->jets[l0].f;
  dbl dist = i0 % 2 == 0 ? SQRT2 : 1;
  dbl T = T0 + sjs->h*dist*(s + s0)/2;
  jet *J = &sjs->jets[l];

  bool updated = false;
  if (T < J->f) {
    updated = true;

    J->f = T;
    J->fx = -s*offsets[i0].i/dist;
    J->fy = -s*offsets[i0].j/dist;

    int lf = sjs->node_parents[l];
    if (lf != NO_PARENT) {
      dvec2 xy = sjs_xy(sjs, l), xyf = sjs_xy(sjs, lf);
      dbl rf = dvec2_dist(xy, xyf);
      J->f -= rf;
      J->fx -= (xy.x - xyf.x)/rf;
      J->fy -= (xy.y - xyf.y)/rf;
    }
  }
  return updated;
}

bool sjs_cell_is_ready(sjs_s *sjs, int lc) {
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = sjs_lc2l(sjs, lc);
    state_e state = sjs->states[l + sjs->vert_dl[i]];
    if (state != TRIAL && state != VALID) {
      return false;
    }
  }
  return true;
}

void sjs_est_Txy_values(sjs_s *sjs, int lc, dbl Txy[4]) {
  dbl fx[NUM_CELL_VERTS], fy[NUM_CELL_VERTS];

  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = sjs_lc2l(sjs, lc) + sjs->vert_dl[i];
    fx[i] = sjs->jets[l].fx;
    fy[i] = sjs->jets[l].fy;
  }

  dbl fxy[NUM_CELL_VERTS] = {
    (fy[1] - fy[0])/sjs->h, // left
    (fx[3] - fx[1])/sjs->h, // bottom
    (fx[2] - fx[0])/sjs->h, // top
    (fy[3] - fy[2])/sjs->h // right
  };

  static dbl lams[4] = {-0.5, 0.5, 0.5, 1.5};
  static dbl mus[4] = {0.5, -0.5, 1.5, 0.5};

  dbl lam, mu;

  for (int i = 0; i < 4; ++i) {
    lam = lams[i];
    mu = mus[i];
    Txy[i] = (1 - mu)*((1 - lam)*fxy[0] + lam*fxy[1]) +
      mu*((1 - lam)*fxy[2] + lam*fxy[3]);
  }
}

void sjs_update_cell(sjs_s *sjs, int lc) {
  int l[4];
  jet *J[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    l[i] = sjs_lc2l(sjs, lc) + sjs->vert_dl[i];
    J[i] = &sjs->jets[l[i]];
  }

  dbl data[4][4];

  dbl h = sjs->h, h_sq = h*h;

  // TODO: this is a screwy way to do this---we want to lay this out
  // in memory so it will eventually be easy to just do this using
  // SIMD instructions
  data[0][0] = J[0]->f;
  data[0][1] = J[2]->f;
  data[0][2] = h*J[0]->fy;
  data[0][3] = h*J[2]->fy;
  data[1][0] = J[1]->f;
  data[1][1] = J[3]->f;
  data[1][2] = h*J[1]->fy;
  data[1][3] = h*J[3]->fy;
  data[2][0] = h*J[0]->fx;
  data[2][1] = h*J[2]->fx;
  data[2][2] = h_sq*J[0]->fxy;
  data[2][3] = h_sq*J[2]->fxy;
  data[3][0] = h*J[1]->fx;
  data[3][1] = h*J[3]->fx;
  data[3][2] = h_sq*J[1]->fxy;
  data[3][3] = h_sq*J[3]->fxy;

  bicubic_set_A(&sjs->bicubics[lc], data);
}

static bool update(sjs_s *sjs, int l) {
  bool done[NUM_NB], updated = false;
  memset(done, 0x0, NUM_NB*sizeof(bool));
  for (int i = 1, l0, l1; i < 8; i += 2) {
    l0 = l + sjs->nb_dl[i];
    if (sjs->states[l0] == VALID) {
      l1 = l + sjs->nb_dl[i - 1];
      if (sjs->states[l1] == VALID) {
        updated |= tri(sjs, l, l0, l1, i);
        done[i] = done[i - 1] = true;
      }
      l1 = l + sjs->nb_dl[i + 1];
      if (sjs->states[l1] == VALID) {
        updated |= tri(sjs, l, l0, l1, i);
        done[i] = done[(i + 1) % NUM_NB] = true;
      }
    }
  }
  for (int i = 0, l0; i < 8; ++i) {
    l0 = l + sjs->nb_dl[i];
    if (!done[i] && sjs->states[l0] == VALID) {
      updated |= line(sjs, l, l0, i);
    }
  }
  return updated;
}

void adjust(sjs_s *sjs, int l0) {
  assert(sjs->states[l0] == TRIAL);
  assert(l0 >= 0);
  assert(l0 < sjs->nnodes);

  heap_swim(sjs->heap, sjs->positions[l0]);
}

static bool inbounds(sjs_s *sjs, int l) {
  return 0 <= l && l < sjs->nnodes;
}

#define _ NO_INDEX
static int dirty2updated[NUM_NEARBY_CELLS][NUM_CELL_VERTS - 1] = {
  {0, _, _}, {0, 1, _}, {1, 2, _}, {2, _, _},
  {0, 7, _}, {0, 1, 7}, {1, 2, 3}, {2, 3, _},
  {6, 7, _}, {5, 6, 7}, {3, 4, 5}, {3, 4, _},
  {6, _, _}, {5, 6, _}, {4, 5, _}, {4, _, _}
};
#undef _

static int updated2dirty[NUM_NB][NUM_CELL_VERTS] = {
  {0, 1, 4, 5},
  {1, 2, 5, 6},
  {2, 3, 6, 7},
  {6, 7, 10, 11},
  {10, 11, 14, 15},
  {9, 10, 13, 14},
  {8, 9, 12, 13},
  {4, 5, 8, 9}
};

void sjs_step(sjs_s *sjs) {
  int l0 = heap_front(sjs->heap);
  heap_pop(sjs->heap);

  // TODO: this is where we could try rebuilding the adjacent cells,
  // but I'm still not convinced this should be necessary if we do
  // things correctly below...

  // Set FAR nodes to TRIAL and insert them into the heap.
  assert(sjs->states[l0] == TRIAL);
  sjs->states[l0] = VALID;
  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_dl[i];
    if (inbounds(sjs, l) && sjs->states[l] == FAR) {
      sjs->states[l] = TRIAL;
      heap_insert(sjs->heap, l);
    }
  }

  bool updated[NUM_NB];
  memset(updated, 0x0, NUM_NB*sizeof(bool));

  // Update neighboring nodes and adjust their position in the
  // heap. Keep track of which nodes were updated.
  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_dl[i];
    if (inbounds(sjs, l) && sjs->states[l] == TRIAL) {
      updated[i] = update(sjs, l);
      adjust(sjs, l);
    }
  }

  // Find adjacent "dirty" cells ("ready" cells which have at least
  // one updated vertex).
  bool dirty[NUM_NEARBY_CELLS];
  memset(dirty, 0x0, NUM_NEARBY_CELLS*sizeof(bool));
  for (int i = 0, lc; i < NUM_NEARBY_CELLS; ++i) {
    lc = sjs_l2lc(sjs, l0) + sjs->nearby_dlc[i];
    if (!sjs_cell_is_ready(sjs, lc)) {
      continue;
    }
    for (int j = 0, k; j < NUM_CELL_VERTS; ++j) {
      k = dirty2updated[i][j];
      dirty[i] |= k != NO_INDEX && updated[k];
    }
  }

  // Estimate Txy values at the vertices of each dirty cell.
  //
  // TODO: this is overkill---we won't actually use all these
  // values... but it may be more work to try and find the values we
  // don't need...
  dbl Txy[NUM_NEARBY_CELLS][4];
  for (int i = 0, lc; i < NUM_NEARBY_CELLS; ++i) {
    if (dirty[i]) {
      lc = sjs_l2lc(sjs, l0) + sjs->nearby_dlc[i];
      sjs_est_Txy_values(sjs, lc, Txy[i]);
    }
  }

  // Set new Txy values at updated nodes to be the average of the Txy
  // values estimated at each dirty cell.
  dbl Txy_sum;
  int ndirty;
  for (int i = 0, l; i < NUM_NB; ++i) {
    Txy_sum = 0;
    ndirty = 0;
    for (int j = 0, k; j < NUM_CELL_VERTS; ++j) {
      k = updated2dirty[i][j];
      if (dirty[k]) {
        Txy_sum += Txy[k][j];
        ++ndirty;
      }
    }
    if (ndirty > 0) {
      l = l0 + sjs->nb_dl[i];
      sjs->jets[l].fxy = Txy_sum/ndirty;
    }
  }

  // Update dirty cells.
  for (int i = 0, lc; i < NUM_NEARBY_CELLS; ++i) {
    if (dirty[i]) {
      lc = sjs_l2lc(sjs, l0) + sjs->nearby_dlc[i];
      sjs_update_cell(sjs, lc);
    }
  }
}

void sjs_solve(sjs_s *sjs) {
  while (heap_size(sjs->heap) > 0) {
    sjs_step(sjs);
  }
}

int sjs_xy_to_lc_and_cc(sjs_s *sjs, dvec2 xy, dvec2 *cc) {
  assert(cc != NULL);
  cc->x = (xy.x - sjs->xymin.x)/sjs->h;
  cc->y = (xy.y - sjs->xymin.y)/sjs->h;
  ivec2 ind = {(int) floor(cc->x), (int) floor(cc->y)};
  cc->x = fmod(cc->x, 1.0);
  cc->y = fmod(cc->y, 1.0);
  if (ind.i == sjs->shape.i - 1) {
    --ind.i;
    cc->x = 1.0;
  }
  if (ind.j == sjs->shape.j - 1) {
    --ind.j;
    cc->y = 1.0;
  }
  int l = sjs_lindexe(sjs, ind);
  return sjs_l2lc(sjs, l);
}

dbl sjs_T(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = sjs_xy_to_lc_and_cc(sjs, xy, &cc);
  bicubic *bicubic = &sjs->bicubics[lc];
  dbl T = bicubic_f(bicubic, cc);
  int lf = sjs->cell_parents[lc];
  if (lf != NO_PARENT) {
    dvec2 xyf = sjs_xy(sjs, lf);
    T += dvec2_dist(xy, xyf);
  }
  return T;
}

dbl sjs_Tx(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = sjs_xy_to_lc_and_cc(sjs, xy, &cc);
  bicubic *bicubic = &sjs->bicubics[lc];
  dbl Tx = bicubic_fx(bicubic, cc);
  int lf = sjs->cell_parents[lc];
  if (lf != NO_PARENT) {
    dvec2 xyf = sjs_xy(sjs, lf);
    Tx += (xy.x - xyf.x)/dvec2_dist(xy, xyf);
  }
  return Tx;
}

dbl sjs_Ty(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = sjs_xy_to_lc_and_cc(sjs, xy, &cc);
  bicubic *bicubic = &sjs->bicubics[lc];
  dbl Ty = bicubic_fy(bicubic, cc);
  int lf = sjs->cell_parents[lc];
  if (lf != NO_PARENT) {
    dvec2 xyf = sjs_xy(sjs, lf);
    Ty += (xy.y - xyf.y)/dvec2_dist(xy, xyf);
  }
  return Ty;
}

dbl sjs_Txy(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = sjs_xy_to_lc_and_cc(sjs, xy, &cc);
  bicubic *bicubic = &sjs->bicubics[lc];
  dbl Txy = bicubic_fxy(bicubic, cc);
  int lf = sjs->cell_parents[lc];
  if (lf != NO_PARENT) {
    dvec2 xyf = sjs_xy(sjs, lf);
    Txy += -(xy.x - xyf.x)*(xy.y - xyf.y)/pow(dvec2_dist(xy, xyf), 3);
  }
  return Txy;
}

bicubic *sjs_bicubic(sjs_s *sjs, ivec2 cind) {
  int lc = sjs_clindexe(sjs, cind);
  return &sjs->bicubics[lc];
}
