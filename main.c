#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef double dbl;

typedef struct {
  dbl A[16];
} cell;

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

typedef struct sjs_ sjs;

typedef struct {
  int capacity;
  int size;
  int* inds;
  sjs *sjs;
} heap;

typedef struct sjs_ {
  ivec2 shape;
  dbl h;
  int nbs[8];
  func *s;
  cell *cells;
  jet *jets;
  state *states;
  int *parents;
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

void swap(heap *heap, int pos1, int pos2) {
  int tmp = heap->inds[pos1];
  heap->inds[pos1] = heap->inds[pos2];
  heap->inds[pos2] = tmp;
}

void swim(heap *heap, int pos) {
  int par = parent(pos);
  // TODO: this calls `value` about 2x as many times as necessary
  while (pos > 0 && value(heap, par) > value(heap, pos)) {
    swap(heap, par, pos);
    pos = par;
    par = parent(pos);
  }
}

void heap_insert(heap *heap, int ind) {
  if (heap->size == heap->capacity) {
    heap_grow(heap);
  }
  int pos = heap->size++;
  heap->inds[pos] = ind;
  swim(heap, pos);
}

int heap_front(heap *heap) {
  return heap->inds[0];
}

void sink(heap *heap, int pos) {
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
      swap(heap, pos, ch);
    }
    pos = ch;
    ch = left(pos);
    next = ch + 1;
  }
}

void heap_pop(heap *heap) {
  if (--heap->size > 0) {
    swap(heap, 0, heap->size);
    sink(heap, 0);
  }
}

int sjs_lindex(sjs *sjs, int i, int j) {
  return (sjs->shape.i + 2)*(j + 1) + i + 1;
}

void sjs_set_nb_inds(sjs *sjs) {
  static int offsets[8][2] = {
    {-1, -1}, {-1,  0}, {-1,  1},
    { 0, -1},           { 0,  1},
    { 1, -1}, { 1,  0}, { 1, -1}
  };
  for (int i = 0; i < 8; ++i) {
    sjs->nbs[i] = sjs_lindex(sjs, offsets[i][0], offsets[i][1]);
  }
}

void sjs_init(sjs *sjs, ivec2 shape, dbl h, func *s) {
  int m = shape.i, n = shape.j;
  int ncells = (m + 1)*(n + 1);
  int nnodes = (m + 2)*(n + 2);

  sjs->shape = shape;
  sjs->h = h;
  sjs->s = s;
  sjs->cells = malloc(ncells*sizeof(cell));
  sjs->jets = malloc(nnodes*sizeof(jet));
  sjs->states = malloc(nnodes*sizeof(state));
  sjs->parents = malloc(nnodes*sizeof(int));

  sjs_set_nb_inds(sjs);

  for (int l = 0; l < nnodes; ++l) {
    sjs->states[l] = FAR;
  }
}

void sjs_add_fac_pt_src(sjs *sjs, int i0, int j0, dbl r0) {
  int m = sjs->shape.i, n = sjs->shape.j;

  int l0 = sjs_lindex(sjs, i0, j0);
  for (int i = 0; i < m; ++i) {
    dbl x = ((dbl) i)/((dbl) (m - 1));
    for (int j = 0; j < n; ++j) {
      dbl y = ((dbl) j)/((dbl) (n - 1));
      int l = sjs_lindex(sjs, i, j);
      sjs->parents[l] = hypot(x, y) <= r0 ? l0 : -1;
    }
  }

  jet *J = &sjs->jets[l0];
  J->f = J->fx = J->fy = J->fxy = 0;
  sjs->states[l0] = TRIAL;
  heap_insert(&sjs->heap, l0);
}

void sjs_step(sjs *sjs) {
  int l = heap_front(&sjs->heap);
  heap_pop(&sjs->heap);
  sjs->states[l] = VALID;


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

  int i0 = m/2, j0 = n/2;

  ivec2 shape = {.i = m, .j = n};
  func s = {.f = f, .df = df};

  sjs sjs;
  sjs_init(&sjs, shape, h, &s);
  sjs_add_fac_pt_src(&sjs, i0, j0, rf);
  sjs_solve(&sjs);

  free(sjs.cells);
  free(sjs.jets);
  free(sjs.states);
}
