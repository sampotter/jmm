#include "sjs_eik.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "index.h"
#include "jet.h"
#include "math.h"

#define NUM_CELL_VERTS 4
#define NUM_NB 8
#define NUM_NEARBY_CELLS 16

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
  sfield_f s;
  vfield_f grad_s;
  bicubic_s *bicubics;
  jet_s *jets;
  state_e *states;
  int *positions;
  heap_s *heap;
  void *context;
};

typedef struct {
  sjs_s *sjs;
  bicubic_variable var;
  cubic_s cubic;
  dvec2 xy0, xy1;
} F_data;

static ivec2 offsets[NUM_NB + 1] = {
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

static ivec2 cell_vert_offsets[NUM_CELL_VERTS] = {
  {.i = 0, .j = 0},
  {.i = 1, .j = 0},
  {.i = 0, .j = 1},
  {.i = 1, .j = 1}
};

static ivec2 tri_cell_offsets[NUM_NB] = {
  {.i = -2, .j = -1},
  {.i = -2, .j =  0},
  {.i = -1, .j =  1},
  {.i =  0, .j =  1},
  {.i =  1, .j =  0},
  {.i =  1, .j = -1},
  {.i =  0, .j = -2},
  {.i = -1, .j = -2}
};

static ivec2 nb_cell_offsets[NUM_CELL_VERTS] = {
  {.i = -1, .j = -1},
  {.i =  0, .j = -1},
  {.i = -1, .j =  0},
  {.i =  0, .j =  0}
};

static ivec2 nearby_cell_offsets[NUM_NEARBY_CELLS] = {
  {.i = -2, .j = -2}, {.i = -2, .j = -1}, {.i = -2, .j = 0}, {.i = -2, .j = 1},
  {.i = -1, .j = -2}, {.i = -1, .j = -1}, {.i = -1, .j = 0}, {.i = -1, .j = 1},
  {.i =  0, .j = -2}, {.i =  0, .j = -1}, {.i =  0, .j = 0}, {.i =  0, .j = 1},
  {.i =  1, .j = -2}, {.i =  1, .j = -1}, {.i =  1, .j = 0}, {.i =  1, .j = 1},
};

static bicubic_variable tri_bicubic_vars[NUM_NB] = {
  MU, MU, LAMBDA, LAMBDA, MU, MU, LAMBDA, LAMBDA
};

static int tri_edges[NUM_NB] = {1, 1, 0, 0, 0, 0, 1, 1};

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

static void set_nb_dl(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    sjs->nb_dl[i] = ind2l(sjs->shape, offsets[i]);
  }
}

static void set_vert_dl(sjs_s *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->vert_dl[i] = ind2l(sjs->shape, cell_vert_offsets[i]);
  }
}

static void set_tri_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB; ++i) {
    sjs->tri_dlc[i] = ind2lc(sjs->shape, tri_cell_offsets[i]);
  }
}

static void set_nb_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->nb_dlc[i] = ind2lc(sjs->shape, nb_cell_offsets[i]);
  }
}

static void set_nearby_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_NEARBY_CELLS; ++i) {
    sjs->nearby_dlc[i] = ind2lc(sjs->shape, nearby_cell_offsets[i]);
  }
}

static dvec2 get_xy(sjs_s *sjs, int l) {
  int mpad = sjs->shape.i + 2;
  dvec2 xy = {
    .x = sjs->xymin.x + sjs->h*(l%mpad - 1),
    .y = sjs->xymin.y + sjs->h*(l/mpad - 1)
  };
  return xy;
}

static dbl get_s(sjs_s *sjs, int l) {
  return sjs->s(sjs->context, get_xy(sjs, l));
}

static dbl F(F_data *data, dbl t) {
  sjs_s *sjs = data->sjs;
  dvec2 xyt = dvec2_ccomb(data->xy0, data->xy1, t);
  dbl T = cubic_f(&data->cubic, t);
  dbl s = sjs->s(sjs->context, xyt);
  dbl L = sqrt(1 + t*t);
  return T + sjs->h*s*L;
}

static dbl dF_dt(F_data *data, dbl t) {
  sjs_s *sjs = data->sjs;
  dvec2 xyt = dvec2_ccomb(data->xy0, data->xy1, t);
  dbl s = sjs->s(sjs->context, xyt);
  dvec2 ds = sjs->grad_s(sjs->context, xyt);
  dbl ds_dt = data->var == LAMBDA ? ds.x : ds.y;
  dbl dT_dt = cubic_df(&data->cubic, t);
  dbl L = sqrt(1 + t*t);
  dbl dL_dt = t/L;
  return dT_dt + sjs->h*(ds_dt*L + s*dL_dt);
}

static bool line(sjs_s *sjs, int l, int l0, int i0) {
  dbl s = get_s(sjs, l), s0 = get_s(sjs, l0);
  dbl T0 = sjs->jets[l0].f;
  dbl dist = i0 % 2 == 0 ? SQRT2 : 1;
  dbl T = T0 + sjs->h*dist*(s + s0)/2;
  jet_s *J = &sjs->jets[l];

  bool updated = false;
  if (T < J->f) {
    updated = true;

    J->f = T;
    J->fx = -s*offsets[i0].i/dist;
    J->fy = -s*offsets[i0].j/dist;
  }
  return updated;
}

static bool tri(sjs_s *sjs, int l, int l0, int l1, int i0) {
  assert(i0 >= 0);
  assert(i0 < NUM_NB);

  int lc = l2lc(sjs->shape, l) + sjs->tri_dlc[i0];
  assert(lc >= 0);
  assert(lc < sjs->ncells);

  bicubic_s *bicubic = &sjs->bicubics[lc];

  F_data data;
  data.sjs = sjs;
  data.var = tri_bicubic_vars[i0];
  data.cubic = bicubic_restrict(bicubic, data.var, tri_edges[i0]);
  data.xy0 = get_xy(sjs, l0);
  data.xy1 = get_xy(sjs, l1);

  dbl lam, a, b, c, d, fa, fb, fc, fd, dm, df, ds, dd, tmp;

  bool updated = false;

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
    jet_s *J = &sjs->jets[l];
    if (T < J->f) {
      updated = true;

      J->f = T;
      dvec2 xy = get_xy(sjs, l);
      dvec2 xylam = dvec2_ccomb(data.xy0, data.xy1, lam);
      dbl L = sqrt(1 + lam*lam);
      J->fx = get_s(sjs, l)*(xy.x - xylam.x)/L;
      J->fy = get_s(sjs, l)*(xy.y - xylam.y)/L;
    }
  }
  return updated;
}

static bool can_build_cell(sjs_s const *sjs, int lc) {
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(sjs->shape, lc);
    state_e state = sjs->states[l + sjs->vert_dl[i]];
    if (state != TRIAL && state != VALID) {
      return false;
    }
  }
  return true;
}

static void est_Txy_values(sjs_s *sjs, int lc, dbl Txy[4]) {
  dbl fx[NUM_CELL_VERTS], fy[NUM_CELL_VERTS];

  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(sjs->shape, lc) + sjs->vert_dl[i];
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

/**
 * Build the cell at index `lc`. This uses the jets at the cells
 * vertices to assemble the data matrix for the bicubic
 * interpolant. This function doesn't make any assumptions about the
 * state of the nodes, so it's assumed that the caller has already
 * made sure this is a reasonable thing to try to do.
 */
static void build_cell(sjs_s *sjs, int lc) {
  /* Get linear indices of cell vertices */
  int l[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    l[i] = lc2l(sjs->shape, lc) + sjs->vert_dl[i];
  }

  /* Get jet at each cell vertex */
  jet_s *J[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    J[i] = &sjs->jets[l[i]];
  }

  /* Precompute scaling factors for partial derivatives */
  dbl h = sjs->h, h_sq = h*h;

  /* Compute cell data from jets and scaling factors */
  dmat44 data;
  data.data[0][0] = J[0]->f;
  data.data[1][0] = J[1]->f;
  data.data[0][1] = J[2]->f;
  data.data[1][1] = J[3]->f;
  data.data[2][0] = h*J[0]->fx;
  data.data[3][0] = h*J[1]->fx;
  data.data[2][1] = h*J[2]->fx;
  data.data[3][1] = h*J[3]->fx;
  data.data[0][2] = h*J[0]->fy;
  data.data[1][2] = h*J[1]->fy;
  data.data[0][3] = h*J[2]->fy;
  data.data[1][3] = h*J[3]->fy;
  data.data[2][2] = h_sq*J[0]->fxy;
  data.data[3][2] = h_sq*J[1]->fxy;
  data.data[2][3] = h_sq*J[2]->fxy;
  data.data[3][3] = h_sq*J[3]->fxy;

  /* Set cell data */
  bicubic_set_data(&sjs->bicubics[lc], data);
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

static void adjust(sjs_s *sjs, int l0) {
  assert(sjs->states[l0] == TRIAL);
  assert(l0 >= 0);
  assert(l0 < sjs->nnodes);

  heap_swim(sjs->heap, sjs->positions[l0]);
}

/**
 * Checks if the linear index `l` is valid for `sjs`.
 */
static bool inbounds(sjs_s *sjs, int l) {
  return 0 <= l && l < sjs->nnodes;
}

static dbl value(void *vp, int l) {
  sjs_s *sjs = (sjs_s *)vp;
  dbl T = sjs->jets[l].f;
  return T;
}

static void setpos(void *vp, int l, int pos) {
  sjs_s *sjs = (sjs_s *)vp;
  sjs->positions[l] = pos;
}

void sjs_alloc(sjs_s **sjs) {
  *sjs = malloc(sizeof(sjs_s));
  assert(*sjs != NULL);
}

void sjs_dealloc(sjs_s **sjs) {
  free(*sjs);
  *sjs = NULL;
}

// TODO: since the margins are BOUNDARY nodes, we actually don't need
// to allocate an extra margin of cells, since they will never be
// initialized (i.e., they will never have all of their vertex nodes
// become VALID because of the margin...)
void sjs_init(sjs_s *sjs, ivec2 shape, dvec2 xymin, dbl h, sfield_f s,
              vfield_f grad_s, void *context) {
  sjs->shape = shape;
  sjs->ncells = (shape.i - 1)*(shape.j - 1);
  sjs->nnodes = shape.i*shape.j;
  sjs->xymin = xymin;
  sjs->h = h;
  sjs->s = s;
  sjs->grad_s = grad_s;
  sjs->bicubics = malloc(sjs->ncells*sizeof(bicubic_s));
  sjs->jets = malloc(sjs->nnodes*sizeof(jet_s));
  sjs->states = malloc(sjs->nnodes*sizeof(state_e));
  sjs->positions = malloc(sjs->nnodes*sizeof(int));
  sjs->context = context;

  assert(sjs->bicubics != NULL);
  assert(sjs->jets != NULL);
  assert(sjs->states != NULL);
  assert(sjs->positions != NULL);

#ifndef NDEBUG
  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->positions[l] = NO_INDEX;
  }
#endif

  heap_alloc(&sjs->heap);

  int capacity = (int) 3*sqrt(sjs->shape.i*sjs->shape.j);
  heap_init(sjs->heap, capacity, value, setpos, (void *)sjs);

  set_nb_dl(sjs);
  set_vert_dl(sjs);
  set_tri_dlc(sjs);
  set_nb_dlc(sjs);
  set_nearby_dlc(sjs);

  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->jets[l].f = INFINITY;
    sjs->states[l] = FAR;
  }
}

void sjs_deinit(sjs_s *sjs) {
  free(sjs->bicubics);
  free(sjs->jets);
  free(sjs->states);
  free(sjs->positions);

  sjs->bicubics = NULL;
  sjs->jets = NULL;
  sjs->states = NULL;
  sjs->positions = NULL;

  heap_deinit(sjs->heap);
  heap_dealloc(&sjs->heap);
}

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
    lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[i];
    if (!can_build_cell(sjs, lc)) {
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
      lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[i];
      est_Txy_values(sjs, lc, Txy[i]);
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

  // Rebuild dirty cells.
  for (int i = 0, lc; i < NUM_NEARBY_CELLS; ++i) {
    if (dirty[i]) {
      lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[i];
      build_cell(sjs, lc);
    }
  }
}

void sjs_solve(sjs_s *sjs) {
  while (heap_size(sjs->heap) > 0) {
    sjs_step(sjs);
  }
}

void sjs_add_trial(sjs_s *sjs, ivec2 ind, jet_s jet) {
  int l = ind2l(sjs->shape, ind);
  sjs->jets[l] = jet;
  sjs->states[l] = TRIAL;
  heap_insert(sjs->heap, l);
}

void sjs_add_valid(sjs_s *sjs, ivec2 ind, jet_s jet) {
  int l = ind2l(sjs->shape, ind);
  sjs->jets[l] = jet;
  sjs->states[l] = VALID;
}

void sjs_make_bd(sjs_s *sjs, ivec2 ind) {
  int l = ind2l(sjs->shape, ind);
  sjs->states[l] = BOUNDARY;
}

ivec2 sjs_get_shape(sjs_s const *sjs) {
  return sjs->shape;
}

jet_s sjs_get_jet(sjs_s *sjs, ivec2 ind) {
  int l = ind2l(sjs->shape, ind);
  return sjs->jets[l];
}

state_e sjs_get_state(sjs_s const *sjs, ivec2 ind) {
  int l = ind2l(sjs->shape, ind);
  return sjs->states[l];
}

state_e *sjs_get_states_ptr(sjs_s const *sjs) {
  return sjs->states;
}

dbl sjs_T(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = xy_to_lc_and_cc(sjs->shape, sjs->xymin, sjs->h, xy, &cc);
  if (!can_build_cell(sjs, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &sjs->bicubics[lc];
  return bicubic_f(bicubic, cc);
}

dbl sjs_Tx(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = xy_to_lc_and_cc(sjs->shape, sjs->xymin, sjs->h, xy, &cc);
  if (!can_build_cell(sjs, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &sjs->bicubics[lc];
  return bicubic_fx(bicubic, cc);
}

dbl sjs_Ty(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = xy_to_lc_and_cc(sjs->shape, sjs->xymin, sjs->h, xy, &cc);
  if (!can_build_cell(sjs, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &sjs->bicubics[lc];
  return bicubic_fy(bicubic, cc);
}

dbl sjs_Txy(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = xy_to_lc_and_cc(sjs->shape, sjs->xymin, sjs->h, xy, &cc);
  if (!can_build_cell(sjs, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &sjs->bicubics[lc];
  return bicubic_fxy(bicubic, cc);
}

bool sjs_can_build_cell(sjs_s const *sjs, ivec2 indc) {
  int lc = indc2lc(sjs->shape, indc);
  return can_build_cell(sjs, lc);
}

void sjs_build_cells(sjs_s *sjs) {
  for (int lc = 0; lc < sjs->ncells; ++lc) {
    if (can_build_cell(sjs, lc)) {
      build_cell(sjs, lc);
    }
  }
}

bicubic_s sjs_get_bicubic(sjs_s const *sjs, ivec2 indc) {
  int lc = indc2lc(sjs->shape, indc);
  return sjs->bicubics[lc];
}

bicubic_s *sjs_get_bicubics_ptr(sjs_s const *sjs) {
  return sjs->bicubics;
}

heap_s *sjs_get_heap(sjs_s const *sjs) {
  return sjs->heap;
}
