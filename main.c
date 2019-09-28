#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define EPS 1e-13
#define NO_INDEX -1
#define NO_PARENT -1
#define SQRT2 1.414213562373095

typedef double dbl;

typedef struct {
  dbl f, fx, fy, fxy;
} jet;

typedef struct {
  dbl x;
  dbl y;
} dvec2;

dvec2 dvec2_ccomb(dvec2 v0, dvec2 v1, dbl t) {
  dvec2 vt = {(1 - t)*v0.x + t*v1.x, (1 - t)*v0.y + t*v1.y};
  return vt;
}

dbl dvec2_dist(dvec2 v0, dvec2 v1) {
  dbl dx = v1.x - v0.x, dy = v1.y - v0.y;
  return sqrt(dx*dx + dy*dy);
}

typedef struct {
  int i;
  int j;
} ivec2;

typedef enum {FAR, TRIAL, VALID, BOUNDARY} state;

#define UNFACTORED -1

typedef struct {
  dbl (*f)(dvec2);
  dvec2 (*df)(dvec2);
} func;

typedef struct {
  dbl A[4][4];
} bicubic;

dbl V_inv[4][4] = {
  { 1,  0,  0,  0},
  { 0,  0,  1,  0},
  {-3,  3, -2, -1},
  { 2, -2,  1,  1}
};

void bicubic_set_A(bicubic *bicubic, dbl data[4][4]) {
  dbl tmp[4][4];
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      tmp[i][j] = 0;
      for (int k = 0; k < 4; ++k) {
        tmp[i][j] += V_inv[i][k]*data[k][j];
      }
    }
  }

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      bicubic->A[i][j] = 0;
      for (int k = 0; k < 4; ++k) {
        bicubic->A[i][j] += tmp[i][k]*V_inv[i][k];
      }
    }
  }
}

typedef struct {
  dbl a[4];
} cubic;

dbl cubic_f(cubic *cubic, dbl lam) {
  dbl *a = cubic->a;
  return a[0] + lam*(a[1] + lam*(a[2] + lam*a[3]));
}

dbl cubic_df(cubic *cubic, dbl lam) {
  dbl *a = cubic->a;
  return a[1] + lam*(2*a[2] + 3*lam*a[3]);
}

typedef enum {LAMBDA, MU} bicubic_variable;

cubic bicubic_restrict(bicubic *bicubic, bicubic_variable var, int edge) {
  cubic cubic;
  if (var == LAMBDA) {
    if (edge == 0) {
      for (int alpha = 0; alpha < 4; ++alpha) {
        cubic.a[alpha] = bicubic->A[alpha][0];
      }
    } else {
      for (int alpha = 0; alpha < 4; ++alpha) {
        cubic.a[alpha] = 0;
        for (int beta = 0; beta < 4; ++beta) {
          cubic.a[alpha] += bicubic->A[alpha][beta];
        }
      }
    }
  } else {
    if (edge == 0) {
      for (int beta = 0; beta < 4; ++beta) {
        cubic.a[beta] = bicubic->A[0][beta];
      }
    } else {
      for (int beta = 0; beta < 4; ++beta) {
        cubic.a[beta] = 0;
        for (int alpha = 0; alpha < 4; ++alpha) {
          cubic.a[beta] += bicubic->A[alpha][beta];
        }
      }
    }
  }
  return cubic;
}

struct sjs;

typedef struct heap {
  int capacity;
  int size;
  int* inds;
  struct sjs *sjs;
} heap_s;

#define NUM_NB 8
#define NUM_CELL_VERTS 4

typedef struct sjs {
  ivec2 shape;
  dbl h;
  int nnodes, ncells;
  int nb_ind_offsets[NUM_NB + 1];
  int tri_cell_ind_offsets[NUM_NB];
  int cell_vert_ind_offsets[NUM_CELL_VERTS];
  int nb_cell_ind_offsets[NUM_CELL_VERTS];
  func *s;
  bicubic *bicubics;
  jet *jets;
  state *states;
  int *parents;
  int *positions;
  heap_s heap;
} sjs_s;

void heap_init(heap_s *heap, sjs_s *sjs) {
  heap->capacity = (int) 3*sqrt(sjs->shape.i*sjs->shape.j);
  heap->size = 0;
  heap->inds = malloc(heap->capacity*sizeof(int));
  assert(heap->inds != NULL);
#ifndef NDEBUG
  for (int i = 0; i < heap->capacity; ++i) {
    heap->inds[i] = NO_INDEX;
  }
#endif
  heap->sjs = sjs;
}

void heap_grow(heap_s *heap) {
  heap->capacity *= 2;
  heap->inds = realloc(heap->inds, heap->capacity);
  assert(heap->inds != NULL);
#ifndef NDEBUG
  for (int i = heap->size; i < heap->capacity; ++i) {
    heap->inds[i] = NO_INDEX;
  }
#endif
}

int left(int pos) {
  return 2*pos + 1;
}

int right(int pos) {
  return 2*pos + 2;
}

int parent(int pos) {
  return (pos - 1)/2;
}

dbl value(heap_s *heap, int pos) {
  assert(pos >= 0);
  assert(pos < heap->size);

  int ind = heap->inds[pos];
  assert(ind != NO_INDEX);

  return heap->sjs->jets[heap->inds[pos]].f;
}

void heap_set(heap_s *heap, int pos, int ind) {
  assert(pos >= 0);
  assert(pos < heap->size);
  assert(ind >= 0);
  assert(ind < heap->sjs->nnodes);

  heap->inds[pos] = ind;
  heap->sjs->positions[ind] = pos;
}

void heap_swap(heap_s *heap, int pos1, int pos2) {
  assert(pos1 >= 0);
  assert(pos1 < heap->size);
  assert(pos2 >= 0);
  assert(pos2 < heap->size);

  int tmp = heap->inds[pos1];
  heap->inds[pos1] = heap->inds[pos2];
  heap->inds[pos2] = tmp;

  heap_set(heap, pos1, heap->inds[pos1]);
  heap_set(heap, pos2, heap->inds[pos2]);
}

// TODO: this calls `value` and `heap_set` about 2x as many times as
// necessary
void heap_swim(heap_s *heap, int pos) {
  assert(pos >= 0);
  assert(pos < heap->size);

  int par = parent(pos);
  while (pos > 0 && value(heap, par) > value(heap, pos)) {
    heap_swap(heap, par, pos);
    pos = par;
    par = parent(pos);
  }
}

void heap_insert(heap_s *heap, int ind) {
  assert(ind >= 0);
  assert(ind < heap->sjs->nnodes);

  if (heap->size == heap->capacity) {
    heap_grow(heap);
  }

  int pos = heap->size++;
  heap_set(heap, pos, ind);
  heap_swim(heap, pos);
}

int heap_front(heap_s *heap) {
#ifndef NDEBUG
  int ind = heap->inds[0];
  assert(ind >= 0);
  assert(ind < heap->sjs->nnodes);
  return ind;
#else
  return heap->inds[0];
#endif
}

void heap_sink(heap_s *heap, int pos) {
  assert(pos >= 0);
  assert(pos < heap->size);

  int ch = left(pos), next = ch + 1, n = heap->size;
  dbl cval, nval;
  while (ch < n) {
    cval = value(heap, ch);
    if (next < n) {
      nval = value(heap, next);
      if (cval > nval) {
        ch = next;
        cval = nval;
      }
    }
    if (value(heap, pos) > cval) {
      heap_swap(heap, pos, ch);
    }
    pos = ch;
    ch = left(pos);
    next = ch + 1;
  }
}

void heap_pop(heap_s *heap) {
#ifndef NDEBUG
  heap->sjs->positions[heap->inds[0]] = NO_INDEX;
#endif
  heap_swap(heap, 0, heap->size - 1);
  if (--heap->size > 0) {
    heap_sink(heap, 0);
  }
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

void sjs_init(sjs_s *sjs, ivec2 shape, dbl h, func *s) {
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

  heap_init(&sjs->heap, sjs);

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
  heap_insert(&sjs->heap, l0);

  return nf;
}

dbl sjs_get_s(sjs_s *sjs, int l) {
  return sjs->s->f(sjs_xy(sjs, l));
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
  dbl s = data->sjs->s->f(xyt);
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
  dbl s = data->sjs->s->f(xyt);
  dvec2 ds = data->sjs->s->df(xyt);
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

  heap_swim(&sjs->heap, sjs->positions[l0]);
}

bool sjs_inbounds(sjs_s *sjs, int l) {
  return 0 <= l && l < sjs->nnodes;
}

void sjs_step(sjs_s *sjs) {
  int l0 = heap_front(&sjs->heap);
  printf("%d\n", l0);
  heap_pop(&sjs->heap);
  assert(sjs->states[l0] == TRIAL);
  sjs->states[l0] = VALID;

  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_ind_offsets[i];
    if (sjs_inbounds(sjs, l) && sjs->states[l] == FAR) {
      sjs->states[l] = TRIAL;
      heap_insert(&sjs->heap, l);
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
  while (sjs->heap.size > 0) {
    sjs_step(sjs);
  }
}

dbl f(dvec2 p) {
  return 1.0 + 0.3*p.x - 0.2*p.y;
}

dvec2 df(dvec2 p) {
  (void) p;
  static dvec2 v = {.x = 0.3, .y = -0.2};
  return v;
}

void usage() {
  printf("hello\n");
}

int main(int argc, char *argv[]) {
  int m = 101, n = 101, i = -1, j = -1;
  double h = -1, r = 0.1;
  char *path = NULL;
  bool verbose = false;

  char c;
  while ((c = getopt(argc, argv, "m:n:i:j:h:r:o:vH")) != -1) {
    switch (c) {
    case 'm':
      m = atoi(optarg);
      break;
    case 'n':
      n = atoi(optarg);
      break;
    case 'i':
      i = atoi(optarg);
      break;
    case 'j':
      j = atoi(optarg);
      break;
    case 'h':
      h = atof(optarg);
      break;
    case 'r':
      r = atof(optarg);
      break;
    case 'o':
      path = malloc(strlen(optarg) + 1);
      assert(path != NULL);
      strcpy(path, optarg);
      break;
    case 'v':
      verbose = true;
      break;
    case 'H':
      usage();
      exit(EXIT_SUCCESS);
    }
  }

  if (i == -1) {
    i = m/2;
  }
  if (j == -1) {
    j = n/2;
  }
  if (h == -1) {
    h = 1.0/(m - 1);
  }

  if (verbose) {
    printf("(m, n) = (%d, %d), ", m, n);
    printf("h = %g, ", h);
    printf("(i, j) = (%d, %d), ", i, j);
    printf("r = %g, ", r);
    printf("path: %s\n", path ? path : "A.bin");
  }

  ivec2 ind = {i, j};
  ivec2 shape = {m, n};
  func s = {.f = f, .df = df};

  sjs_s sjs;
  sjs_init(&sjs, shape, h, &s);

  int nf = sjs_add_fac_pt_src(&sjs, ind, r);
  if (nf == 0) {
    printf("ERROR: nf = 0\n");
    exit(EXIT_FAILURE);
  } else {
    if (verbose) {
      printf("nf = %d\n", nf);
    }
  }

  sjs_solve(&sjs);

  FILE *f = fopen(path ? path : "A.bin", "w");
  for (int i = 0; i < m; ++i) {
    for (int j = 0, l; j < n; ++j) {
      ivec2 ind = {i, j};
      l = sjs_lindexe(&sjs, ind);
      fwrite(sjs.bicubics[l].A, sizeof(dbl), 16, f);
    }
  }
  fclose(f);

  if (path != NULL) {
    free(path);
  }

  free(sjs.bicubics);
  free(sjs.jets);
  free(sjs.states);
}
