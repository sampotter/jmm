#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define EPS 1e-13
#define SQRT2 1.414213562373095

typedef double dbl;

typedef struct {
  dbl f, fx, fy, fxy;
} jet;

typedef struct {
  dbl x;
  dbl y;
} dvec2;

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

struct sjs_;

typedef struct {
  int capacity;
  int size;
  int* inds;
  struct sjs_ *sjs;
} heap;

#define NUM_NB 8
#define NUM_CELL_VERTS 4

typedef struct sjs_ {
  ivec2 shape;
  dbl h;
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
  heap heap;
} sjs;

void heap_init(heap *heap, int capacity) {
  heap->capacity = capacity;
  heap->size = 0;
  heap->inds = malloc(heap->capacity*sizeof(int));
}

void heap_grow(heap *heap) {
  heap->capacity *= 2;
  heap->inds = realloc(heap->inds, heap->capacity);
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

dbl value(heap *heap, int pos) {
  return heap->sjs->jets[heap->inds[pos]].f;
}

void heap_set(heap *heap, int pos, int ind) {
  heap->inds[pos] = ind;
  heap->sjs->positions[ind] = pos;
}

void heap_swap(heap *heap, int pos1, int pos2) {
  int tmp = heap->inds[pos1];
  heap->inds[pos1] = heap->inds[pos2];
  heap->inds[pos2] = tmp;

  heap_set(heap, pos1, heap->inds[pos1]);
  heap_set(heap, pos2, heap->inds[pos2]);
}

void heap_swim(heap *heap, int pos) {
  int par = parent(pos);
  // TODO: this calls `value` and `heap_set` about 2x as many times as
  // necessary
  while (pos > 0 && value(heap, par) > value(heap, pos)) {
    heap_swap(heap, par, pos);
    pos = par;
    par = parent(pos);
  }
}

void heap_insert(heap *heap, int ind) {
  if (heap->size == heap->capacity) {
    heap_grow(heap);
  }
  int pos = heap->size++;
  heap_set(heap, pos, ind);
  heap_swim(heap, pos);
}

int heap_front(heap *heap) {
  return heap->inds[0];
}

void heap_sink(heap *heap, int pos) {
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

void heap_pop(heap *heap) {
  if (--heap->size > 0) {
    heap_swap(heap, 0, heap->size);
    heap_sink(heap, 0);
  }
}

int sjs_lindex(sjs *sjs, ivec2 ind) {
  return (sjs->shape.i + 2)*(ind.j + 1) + ind.i + 1;
}

void sjs_vindex(sjs *sjs, int l, int *i, int *j) {
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

void sjs_set_nb_ind_offsets(sjs *sjs) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    sjs->nb_ind_offsets[i] = sjs_lindex(sjs, offsets[i]);
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

void sjs_set_tri_cell_ind_offsets(sjs *sjs) {
  for (int i = 0; i < NUM_NB; ++i) {
    sjs->tri_cell_ind_offsets[i] = sjs_lindex(sjs, tri_cell_offsets[i]);
  }
}

ivec2 cell_vert_offsets[NUM_CELL_VERTS] = {
  {.i = 0, .j = 0},
  {.i = 1, .j = 0},
  {.i = 0, .j = 1},
  {.i = 1, .j = 1}
};

void sjs_set_cell_vert_ind_offsets(sjs *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->cell_vert_ind_offsets[i] = sjs_lindex(sjs, cell_vert_offsets[i]);
  }
}

ivec2 nb_cell_offsets[NUM_CELL_VERTS] = {
  {.i = -1, .j = -1},
  {.i =  0, .j = -1},
  {.i = -1, .j =  0},
  {.i =  0, .j = -1}
};

void sjs_set_nb_cell_ind_offsets(sjs *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->nb_cell_ind_offsets[i] = sjs_lindex(sjs, nb_cell_offsets[i]);
  }
}

void sjs_init(sjs *sjs, ivec2 shape, dbl h, func *s) {
  int m = shape.i, n = shape.j;
  int ncells = (m + 1)*(n + 1);
  int nnodes = (m + 2)*(n + 2);

  sjs->shape = shape;
  sjs->h = h;
  sjs->s = s;
  sjs->bicubics = malloc(ncells*sizeof(bicubic));
  sjs->jets = malloc(nnodes*sizeof(jet));
  sjs->states = malloc(nnodes*sizeof(state));
  sjs->parents = malloc(nnodes*sizeof(int));
  sjs->positions = malloc(nnodes*sizeof(int));

  sjs_set_nb_ind_offsets(sjs);
  sjs_set_tri_cell_ind_offsets(sjs);
  sjs_set_cell_vert_ind_offsets(sjs);
  sjs_set_nb_cell_ind_offsets(sjs);

  for (int l = 0; l < nnodes; ++l) {
    sjs->states[l] = FAR;
  }
}

void sjs_add_fac_pt_src(sjs *sjs, ivec2 ind0, dbl r0) {
  int m = sjs->shape.i, n = sjs->shape.j;

  int l0 = sjs_lindex(sjs, ind0);
  for (int i = 0; i < m; ++i) {
    dbl x = ((dbl) i)/((dbl) (m - 1));
    for (int j = 0; j < n; ++j) {
      dbl y = ((dbl) j)/((dbl) (n - 1));
      ivec2 ind = {.i = i, .j = j};
      int l = sjs_lindex(sjs, ind);
      sjs->parents[l] = hypot(x, y) <= r0 ? l0 : -1;
    }
  }

  jet *J = &sjs->jets[l0];
  J->f = J->fx = J->fy = J->fxy = 0;
  sjs->states[l0] = TRIAL;
  heap_insert(&sjs->heap, l0);
}

dvec2 sjs_xy(sjs *sjs, int l) {
  int mpad = sjs->shape.i + 2;
  dvec2 xy = {
    .x = sjs->h*(l/mpad - 1),
    .y = sjs->h*(l%mpad - 1)
  };
  return xy;
}

dbl sjs_s(sjs *sjs, int l) {
  return sjs->s->f(sjs_xy(sjs, l));
}

dbl sjs_T(sjs *sjs, int l) {
  return sjs->jets[l].f;
}

bicubic_variable tri_bicubic_vars[NUM_NB] = {
  MU, MU, LAMBDA, LAMBDA, MU, MU, LAMBDA, LAMBDA
};

int tri_edges[NUM_NB] = {1, 1, 0, 0, 0, 0, 1, 1};

dvec2 get_xylam(dvec2 xy0, dvec2 xy1, dbl lam) {
  dvec2 xylam = {
    .x = (1 - lam)*xy0.x + lam*xy1.x,
    .y = (1 - lam)*xy0.y + lam*xy1.y
  };
  return xylam;
}

typedef struct {
  sjs *sjs;
  bicubic_variable var;
  cubic cubic;
  dvec2 xy0, xy1;
} F_data;

dbl F(F_data *data, dbl lam) {
  dvec2 xylam = get_xylam(data->xy0, data->xy1, lam);
  dbl T = cubic_f(&data->cubic, lam);
  dbl s = data->sjs->s->f(xylam);
  dbl L = sqrt(1 + lam*lam);
  return T + data->sjs->h*s*L;
}

dbl dF_dlam(F_data *data, dbl lam) {
  dvec2 xylam = get_xylam(data->xy0, data->xy1, lam);
  dbl s = data->sjs->s->f(xylam);
  dvec2 ds = data->sjs->s->df(xylam);
  dbl ds_dlam = data->var == LAMBDA ? ds.x : ds.y;
  dbl dT_dlam = cubic_df(&data->cubic, lam);
  dbl L = sqrt(1 + lam*lam);
  dbl dL_dlam = lam/L;
  return dT_dlam + data->sjs->h*(ds_dlam*L + s*dL_dlam);
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

void sjs_tri(sjs *sjs, int l, int l0, int l1, int i0, int i1) {
  F_data data;
  data.sjs = sjs;
  bicubic *bicubic = &sjs->bicubics[l + sjs->tri_cell_ind_offsets[i0]];
  data.var = tri_bicubic_vars[i0];
  data.cubic = bicubic_restrict(bicubic, data.var, tri_edges[i0]);
  data.xy0 = sjs_xy(sjs, l0);
  data.xy1 = sjs_xy(sjs, l1);

  dbl lam, a, b, c, d, fa, fb, fc, fd, dm, df, ds, dd, tmp;

  fa = dF_dlam(&data, 0);
  if (fabs(fa) <= EPS) {
    lam = 0;
    goto found;
  }

  fb = dF_dlam(&data, 1);
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
    fd = dF_dlam(&data, d);
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
      dvec2 xylam = get_xylam(data.xy0, data.xy1, lam);
      dbl L = sqrt(1 + lam*lam);
      J->fx = sjs_s(sjs, l)*(xy.x - xylam.x)/L;
      J->fy = sjs_s(sjs, l)*(xy.y - xylam.y)/L;
    }
  }
}

void sjs_line(sjs *sjs, int l, int l0, int i0) {
  dbl s = sjs_s(sjs, l), s0 = sjs_s(sjs, l0);
  dbl T0 = sjs_T(sjs, l0);
  dbl T = T0 + sjs->h*(s + s0)/2;
  jet *J = &sjs->jets[l];
  if (T < J->f) {
    J->f = T;
    dbl dist = i0 % 2 == 0 ? SQRT2 : 1;
    J->fx = s*offsets[i0].i/dist;
    J->fy = s*offsets[i0].j/dist;
  }
}

bool sjs_valid_cell(sjs *sjs, int l) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    if (sjs->states[l + sjs->cell_vert_ind_offsets[i]] != VALID) {
      return false;
    }
  }
  return true;
}

dbl sjs_est_fxy(sjs *sjs, int ln, int lc) {
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

  int i = 0, dl = ln - lc;
  while (dl != sjs->cell_vert_ind_offsets[i]) ++i;

  dbl lam = lams[i], mu = mus[i];

  return (1 - mu)*((1 - lam)*fxy[0] + lam*fxy[1]) +
    mu*((1 - lam)*fxy[2] + lam*fxy[3]);
}

void sjs_update(sjs *sjs, int l) {
  bool updated[NUM_NB];
  memset(updated, 0x0, NUM_NB*sizeof(bool));
  for (int i = 1, l0, l1; i < 8; i += 2) {
    l0 = l + sjs->nb_ind_offsets[i];
    if (sjs->states[l0] == VALID) {
      l1 = l + sjs->nb_ind_offsets[i - 1];
      if (sjs->states[l1] == VALID) {
        sjs_tri(sjs, l, l0, l1, i, i - 1);
        updated[l0] = updated[l1] = true;
      }
      l1 = l + sjs->nb_ind_offsets[i + 1];
      if (sjs->states[l1] == VALID) {
        sjs_tri(sjs, l, l0, l1, i, i + 1);
        updated[l0] = updated[l1] = true;
      }
    }
  }
  for (int i = 0, l0; i < 8; ++i) {
    l0 = l + sjs->nb_ind_offsets[i];
    if (!updated[l0] && sjs->states[l0] == VALID) {
      sjs_line(sjs, l, l0, i);
    }
  }
}

void sjs_adjust(sjs *sjs, int l0) {
  heap_swim(&sjs->heap, sjs->positions[l0]);
}

void sjs_step(sjs *sjs) {
  int l0 = heap_front(&sjs->heap);
  heap_pop(&sjs->heap);
  sjs->states[l0] = VALID;

  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_ind_offsets[i];
    if (sjs->states[l] == FAR) {
      sjs->states[l] = TRIAL;
    }
  }

  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_ind_offsets[i];
    if (sjs->states[l] == TRIAL) {
      sjs_update(sjs, l);
      sjs_adjust(sjs, l);
    }
  }
}

void sjs_solve(sjs *sjs) {
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

int main() {
  int m = 51, n = 31;
  dbl h = 1.0/(n - 1);
  dbl rf = 0.1;

  ivec2 ind0 = {m/2, n/2};

  ivec2 shape = {.i = m, .j = n};
  func s = {.f = f, .df = df};

  sjs sjs;
  sjs_init(&sjs, shape, h, &s);
  sjs_add_fac_pt_src(&sjs, ind0, rf);
  sjs_solve(&sjs);

  free(sjs.bicubics);
  free(sjs.jets);
  free(sjs.states);
}
