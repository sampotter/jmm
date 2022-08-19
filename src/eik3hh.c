#include "eik3hh.h"

#include <assert.h>
#include <stdio.h>

#include "bmesh.h"
#include "eik3.h"
#include "eik3hh_branch.h"

struct eik3hh {
  mesh3_s const *mesh;
  dbl c;
  dbl rfac;
  eik3hh_branch_s *root;
};

void eik3hh_alloc(eik3hh_s **hh) {
  *hh = malloc(sizeof(eik3hh_s));
}

static void init(eik3hh_s *hh, mesh3_s const *mesh, dbl c, dbl rfac) {
  hh->mesh = mesh;
  hh->c = c;
  hh->rfac = rfac;
  hh->root = NULL;
}

void eik3hh_init_with_pt_src(eik3hh_s *hh, mesh3_s const *mesh, dbl c,
                             dbl rfac, dbl3 const xsrc) {
  init(hh, mesh, c, rfac);

  eik3hh_branch_alloc(&hh->root);
  eik3hh_branch_init_pt_src(hh->root, hh, xsrc);
}

void eik3hh_deinit(eik3hh_s *hh) {
  hh->mesh = NULL;
  hh->c = NAN;
  hh->rfac = NAN;

  if (hh->root != NULL) {
    eik3hh_branch_deinit(hh->root, /* free_children = */ true);
    eik3hh_branch_dealloc(&hh->root);
  }
}

void eik3hh_dealloc(eik3hh_s **hh) {
  free(*hh);
  *hh = NULL;
}

mesh3_s const *eik3hh_get_mesh(eik3hh_s const *hh) {
  return hh->mesh;
}

dbl eik3hh_get_rfac(eik3hh_s const *hh) {
  return hh->rfac;
}

eik3hh_branch_s *eik3hh_get_root_branch(eik3hh_s *hh) {
  return hh->root;
}
