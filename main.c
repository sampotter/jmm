#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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

typedef struct sjs_ sjs;

typedef struct {
  int capacity;
  int size;
  int* inds;
  sjs *sjs;
} heap;

#define NUM_NB 8

typedef struct sjs_ {
  ivec2 shape;
  dbl h;
  int nbs[NUM_NB + 1];
  int tri_cell_inds[NUM_NB];
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

void sjs_set_nb_inds(sjs *sjs) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    sjs->nbs[i] = sjs_lindex(sjs, offsets[i]);
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

void sjs_set_tri_cell_inds(sjs *sjs) {
  for (int i = 0; i < NUM_NB; ++i) {
    sjs->tri_cell_inds[i] = sjs_lindex(sjs, tri_cell_offsets[i]);
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

  sjs_set_nb_inds(sjs);
  sjs_set_tri_cell_inds(sjs);

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

void sjs_tri(sjs *sjs, int l, int l0, int l1) {

}


}

void sjs_line(sjs *sjs, int l, int l0) {
  dbl s = sjs_s(sjs, l), s0 = sjs_s(sjs, l0);
  dbl T0 = sjs_T(sjs, l0);
  dbl T = T0 + sjs->h*(s + s0)/2;
  jet *J = &sjs->jets[l];
  if (T < J->f) {
    J->f = T;
    dvec2 grad_T = sjs_est_grad_T(sjs, l, l0, NO_INDEX);
    J->fx = grad_T.x;
    J->fy = grad_T.y;
  }
}

void sjs_update(sjs *sjs, int l) {
  bool updated[NUM_NB];
  memset(updated, 0x0, NUM_NB*sizeof(bool));
  for (int i = 1, l0, l1; i < 8; i += 2) {
    l0 = l + sjs->nbs[i];
    if (sjs->states[l0] == VALID) {
      l1 = l + sjs->nbs[i - 1];
      if (sjs->states[l1] == VALID) {
        sjs_tri(sjs, l, l0, l1);
        updated[l0] = updated[l1] = true;
      }
      l1 = l + sjs->nbs[i + 1];
      if (sjs->states[l1] == VALID) {
        sjs_tri(sjs, l, l0, l1);
        updated[l0] = updated[l1] = true;
      }
    }
  }
  for (int i = 0, l0; i < 8; ++i) {
    l0 = l + sjs->nbs[i];
    if (!updated[l0] && sjs->states[l0] == VALID) {
      sjs_line(sjs, l, l0);
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
    l = l0 + sjs->nbs[i];
    if (sjs->states[l] == FAR) {
      sjs->states[l] = TRIAL;
    }
  }

  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nbs[i];
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
