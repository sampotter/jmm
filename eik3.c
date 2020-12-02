#include "eik3.h"

#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "heap.h"
#include "vec.h"

struct eik3 {
  mesh3_s const *mesh;
  jet3 *jet;
  state_e *state;
  int *pos;
  heap_s *heap;
};

void eik3_alloc(eik3_s **eik) {
  *eik = malloc(sizeof(eik3_s));
}

void eik3_dealloc(eik3_s **eik) {
  free(*eik);
  *eik = NULL;
}

static dbl value(void *ptr, int l) {
  eik3_s *eik = (eik3_s *)ptr;
  assert(l >= 0);
  assert(l < mesh3_nverts(eik->mesh));
  dbl T = eik->jet[l].f;
  return T;
}

static void setpos(void *ptr, int l, int pos) {
  eik3_s *eik = (eik3_s *)ptr;
  eik->pos[l] = pos;
}

void eik3_init(eik3_s *eik, mesh3_s const *mesh) {
  eik->mesh = mesh;

  size_t nverts = mesh3_nverts(mesh);

  eik->jet = malloc(nverts*sizeof(jet3));
  for (size_t l = 0; l < nverts; ++l) {
    eik->jet[l] = (jet3_s) {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN};
  }

  eik->state = malloc(nverts*sizeof(state_e));
  for (size_t l = 0; l < nverts; ++l) {
    eik->state[l] = FAR;
  }

  eik->pos = malloc(nverts*sizeof(int));
  for (size_t l = 0; l < nverts; ++l) {
    eik->pos[l] = NO_INDEX;
  }

  /**
   * When we compute the initial heap capacity, we want to estimate
   * the number of nodes that could comprise the expanding numerical
   * front at any one time. We can't know this ahead of time, so we
   * set it to a constant multiple times (# nodes)^(1/d). In this
   * case, d=3. Even if this is an underestimate, well still reduce
   * the number of times the heap needs to be expanded at solve time.
   */
  int capacity = (int) 3*cbrt(nverts);

  heap_alloc(&eik->heap);
  heap_init(eik->heap, capacity, value, setpos, (void *)eik);


}

void eik3_deinit(eik3_s *eik) {
  free(eik->jet);
  eik->jet = NULL;

  free(eik->state);
  eik->state = NULL;

  free(eik->pos);
  eik->pos = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);
}

static void tetra(eik3_s *eik, size_t l, size_t l0, size_t l1, size_t l2) {
  // TODO: !!!
}

static void update(eik3_s *eik, size_t l, size_t l0) {
  // To find the updates incident on l and l0, first find the
  // tetrahedra that are incident on the edge (l0, l).
  int nec = mesh3_nec(eik->mesh, l0, l);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(eik->mesh, l0, l);

  // Allocate space for the l1 and l2 indices.
  size_t *l1 = malloc(nec*sizeof(size_t));
  size_t *l2 = malloc(nec*sizeof(size_t));

  int nup = 0;
  for (int i = 0; i < nec; ++i) {
    // Get the indices of cell ec[i]'s vertices
    int lv[4];
    mesh3_cv(eik->mesh, ec[i], lv);

    // Find the VALID vertices that *aren't* l0 or l
    int lnew[2];
    int k = 0;
    for (int j = 0; j < 4; ++j) {
      if (lv[j] != l0 && lv[j] != l && eik->state[lv[j]] == VALID) {
        lnew[k++] = lv[j];
      }
    }
    assert(k <= 2);

    // If we found a pair of VALID verts, add a new update, assigning
    // l1 and l2
    if (k == 2) {
      l1[nup] = lnew[0];
      l2[nup++] = lnew[1];
    }
  }

  // Don't need ec any more
  free(ec);

  // Do the updates
  for (int i = 0; i < nup; ++i) {
    tetra(eik, l, l0, l1, l2);
  }

coda:
  free(l1);
  free(l2);
}

static void adjust(eik3_s *eik, size_t l) {
  assert(eik->state[l] == TRIAL);
  assert(l >= 0);
  assert(l < mesh3_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l]);
}

void eik3_step(eik3_s *eik) {
  size_t l0 = heap_front(eik->heap);
  assert(eik->states[l0] == TRIAL);
  heap_pop(eik->heap);
  eik->states[l0] = VALID;

  // Get i0's neighboring nodes.
  int nnb = mesh3_nvv(eik->mesh, l0);
  size_t *nb = malloc(nnb*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, nb);

  // Set FAR nodes to TRIAL and insert them into the heap.
  for (int i = 0; i < nnb; ++i) {
    size_t l = nb[i];
    if (eik->states[l] == FAR) {
      eik->states[l] = TRIAL;
      heap_insert(eik->heap, l);
    }
  }

  // Update neighboring nodes.
  for (int i = 0; i < nnb; ++i) {
    size_t l = nb[i];
    if (eik->states[l] == TRIAL) {
      update(eik, l, l0);
      adjust(eik, l);
    }
  }

  free(nb);
}

void eik3_solve(eik3_s *eik) {
  while (heap_size(eik->heap) > 0) {
    eik3_step(eik);
  }
}

void eik3_add_trial(eik3_s *eik, size_t l, jet3_s jet) {
  eik->jet[l] = jet;

  // TODO: a better way to do this would be to adjust the position of
  // this node in the heap if it's already TRIAL... much better API
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = TRIAL;
  heap_insert(eik->heap, l);
}

void eik3_add_valid(eik3_s *eik, size_t l, jet3_s jet) {
  eik->jets[l] = jet;

  // TODO: see comment in eik3_add_trial... would want to make changes
  // that are compatible (we would need to adjust the API of this a
  // bit to keep things sane)
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = VALID;
}
