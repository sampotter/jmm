#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "sjs.h"

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

  sjs_s *sjs;
  sjs_alloc(&sjs);
  sjs_init(sjs, shape, h, f, df);

  int nf = sjs_add_fac_pt_src(sjs, ind, r);
  if (nf == 0) {
    printf("ERROR: nf = 0\n");
    exit(EXIT_FAILURE);
  } else {
    if (verbose) {
      printf("nf = %d\n", nf);
    }
  }

  sjs_solve(sjs);

  FILE *f = fopen(path ? path : "A.bin", "w");
  for (int i = 0; i < m; ++i) {
    for (int j = 0, l; j < n; ++j) {
      ivec2 ind = {i, j};
      fwrite(sjs_bicubic(sjs, ind)->A, sizeof(dbl), 16, f);
    }
  }
  fclose(f);

  if (path != NULL) {
    free(path);
  }

  sjs_deinit(sjs);
  sjs_dealloc(&sjs);
}
