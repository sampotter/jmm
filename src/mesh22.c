#include "mesh22.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"
#include "util.h"

struct mesh22 {
  dbl2 *verts;
  size_t nverts;
  uint3 *faces;
  size_t nfaces;

  // details:

  size_t *vf_offsets, *vf;      // vertex-to-face LUT
};

void mesh22_alloc(mesh22_s **mesh) {
  *mesh = malloc(sizeof(mesh22_s));
}

void mesh22_dealloc(mesh22_s **mesh) {
  free(*mesh);
  *mesh = NULL;
}

static void init_details(mesh22_s *mesh) {
  // build vertex-to-face LUT:

  /* first, count the number of faces incident on each vertex */
  size_t *nvf = calloc(mesh->nverts, sizeof(size_t));
  for (size_t l = 0; l < mesh->nfaces; ++l) {
    for (size_t i = 0; i < 3; ++i) {
      size_t j = mesh->faces[l][i];
      assert(j < mesh->nverts);
      ++nvf[j];
    }
  }

  /* next, use the number of incident faces to compute offsets into
   * our incidence lists */
  mesh->vf_offsets = malloc(sizeof(size_t)*(mesh->nverts + 1));
  mesh->vf_offsets[0] = 0;
  for (size_t l = 0; l < mesh->nverts; ++l)
    mesh->vf_offsets[l + 1] = mesh->vf_offsets[l] + nvf[l];

  /* we zero out nvf here and recompute it... this is a little
   * inefficient, it'd be nice to figure out how to initialize in one
   * pass... */
  memset(nvf, 0x0, mesh->nverts*sizeof(size_t));

  mesh->vf = malloc(mesh->vf_offsets[mesh->nverts]*sizeof(size_t));
  for (size_t l = 0; l < mesh->nfaces; ++l) {
    for (size_t i = 0; i < 3; ++i) {
      size_t j = mesh->faces[l][i];
      mesh->vf[mesh->vf_offsets[j] + nvf[j]++] = l;
    }
  }

  free(nvf);
}

void mesh22_init(mesh22_s *mesh, dbl2 const *verts, size_t nverts,
                 uint3 const *faces, size_t nfaces) {
  mesh->verts = malloc(nverts*sizeof(dbl2));
  memcpy(mesh->verts, verts, nverts*sizeof(dbl2));

  mesh->nverts = nverts;

  mesh->faces = malloc(nfaces*sizeof(uint3));
  memcpy(mesh->faces, faces, nfaces*sizeof(uint3));

  mesh->nfaces = nfaces;

  init_details(mesh);
}

void mesh22_init_from_binary_files(mesh22_s *mesh, char const *verts_path,
                                   char const *faces_path) {
  FILE *fp = NULL;
  size_t size = 0;

  /**
   * Read points from binary file at `verts_path`
   */

  fp = fopen(verts_path, "rb");
  if (!fp) {
    log_error("failed to open file \"%s\"", verts_path);
    exit(EXIT_FAILURE);
  }

  fseek(fp, 0, SEEK_END);
  size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  log_debug("%s size: %lu bytes", verts_path, size);
  log_debug("%llu", size/sizeof(dbl2));

  mesh->verts = malloc(size);
  mesh->nverts = size/sizeof(dbl2);
  fread(mesh->verts, sizeof(dbl2), mesh->nverts, fp);

  if (fclose(fp)) {
    perror("failed to close file");
    exit(EXIT_FAILURE);
  }

  /**
   * Read faces from binary file at `faces_path`
   */

  fp = fopen(faces_path, "rb");
  if (!fp) {
    log_error("failed to open file \"%s\"", faces_path);
    exit(EXIT_FAILURE);
  }

  fseek(fp, 0, SEEK_END);
  size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  log_debug("%s size: %lu bytes", faces_path, size);

  mesh->faces = malloc(size);
  mesh->nfaces = size/sizeof(uint3);
  fread(mesh->faces, sizeof(uint3), mesh->nfaces, fp);

  if (fclose(fp)) {
    perror("failed to close file");
    exit(EXIT_FAILURE);
  }

  init_details(mesh);
}

void mesh22_deinit(mesh22_s *mesh) {
  free(mesh->verts);
  mesh->verts = NULL;

  free(mesh->faces);
  mesh->faces = NULL;

  // details:

  free(mesh->vf_offsets);
  mesh->vf_offsets = NULL;

  free(mesh->vf);
  mesh->vf = NULL;
}

size_t mesh22_nverts(mesh22_s const *mesh) {
  return mesh->nverts;
}

size_t mesh22_nfaces(mesh22_s const *mesh) {
  return mesh->nfaces;
}

dbl2 const *mesh22_get_verts_ptr(mesh22_s const *mesh) {
  return mesh->verts;
}

uint3 const *mesh22_get_faces_ptr(mesh22_s const *mesh) {
  return mesh->faces;
}

void mesh22_get_vert(mesh22_s const *mesh, size_t l, dbl2 x) {
  dbl2_copy(mesh->verts[l], x);
}

void mesh22_fv(mesh22_s const *mesh, size_t l, uint3 fv) {
  memcpy(fv, mesh->faces[l], sizeof(uint3));
}

size_t mesh22_nvf(mesh22_s const *mesh, size_t l) {
  assert(l < mesh->nverts);
  return mesh->vf_offsets[l + 1] - mesh->vf_offsets[l];
}

void mesh22_vf(mesh22_s const *mesh, size_t l, size_t *vf) {
  assert(l < mesh->nverts);
  size_t nvf = mesh->vf_offsets[l + 1] - mesh->vf_offsets[l];
  memcpy(vf, &mesh->vf[mesh->vf_offsets[l]], nvf*sizeof(size_t));
}

size_t mesh22_nvv(mesh22_s const *mesh, size_t l) {
  assert(l < mesh->nverts);

  size_t nvf = mesh22_nvf(mesh, l);
  size_t *vf = malloc(nvf*sizeof(size_t));
  mesh22_vf(mesh, l, vf);

  /* trivial upper bound: nvv <= 2*nvf */
  size_t *vv = malloc(2*sizeof(size_t)*nvf);

  size_t nvv = 0;
  for (size_t i = 0; i < nvf; ++i) {
    size_t const *lf = mesh->faces[vf[i]];
    for (size_t j = 0; j < 3; ++j) {
      size_t m = lf[j];
      if (l == m || contains(vv, nvv, &m, sizeof(size_t)))
        continue;
      vv[nvv++] = m;
    }
  }

  free(vv);

  free(vf);

  return nvv;
}

void mesh22_vv(mesh22_s const *mesh, size_t l, size_t *vv) {
  assert(l < mesh->nverts);

  size_t nvf = mesh22_nvf(mesh, l);
  size_t *vf = malloc(nvf*sizeof(size_t));
  mesh22_vf(mesh, l, vf);

  size_t nvv = 0;
  for (size_t i = 0; i < nvf; ++i) {
    size_t const *lf = mesh->faces[vf[i]];
    for (size_t j = 0; j < 3; ++j) {
      size_t m = lf[j];
      if (l == m || contains(vv, nvv, &m, sizeof(size_t)))
        continue;
      vv[nvv++] = m;
    }
  }

  free(vf);
}
