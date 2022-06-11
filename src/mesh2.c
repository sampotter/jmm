#include "mesh2.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "def.h"
#include "log.h"
#include "vec.h"

bool mesh2_tri_equal(mesh2_tri_s const *t1, mesh2_tri_s const *t2) {
  return t1->mesh == t2->mesh && t1->l == t2->l;
}

struct mesh2 {
  dbl3 const *verts;
  size_t num_verts;
  bool owns_verts;

  size_t const (*faces)[3];
  size_t num_faces;

  dbl3 const *face_normals;
};

void mesh2_alloc(mesh2_s **mesh) {
  *mesh = malloc(sizeof(mesh2_s));
}

void mesh2_dealloc(mesh2_s **mesh) {
  free(*mesh);
  *mesh = NULL;
}

void mesh2_init(mesh2_s *mesh,
                dbl3 const *verts, size_t nverts, bool copy_verts,
                size_t const (*faces)[3], size_t nfaces,
                dbl3 const *face_normals) {
  if (copy_verts) {
    mesh->verts = malloc(nverts*sizeof(dbl3));
    memcpy((dbl3 *)mesh->verts, verts, nverts*sizeof(dbl3));
  } else {
    mesh->verts = (dbl3 *)verts;
  }
  mesh->num_verts = nverts;
  mesh->owns_verts = copy_verts;

  mesh->faces = malloc(nfaces*sizeof(size_t[3]));
  memcpy((size_t(*)[3])mesh->faces, faces, nfaces*sizeof(size_t[3]));
  mesh->num_faces = nfaces;

  mesh->face_normals = malloc(nfaces*sizeof(dbl3));
  memcpy((dbl3 *)mesh->face_normals, face_normals, nfaces*sizeof(dbl3));
}

void mesh2_init_from_binary_files(mesh2_s *mesh, char const *verts_path,
                                  char const *faces_path) {
  FILE *fp = NULL;
  size_t size = 0;

  /**
   * Read verts from binary file at `verts_path`
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

  mesh->verts = malloc(size);
  mesh->num_verts = size/sizeof(dbl3);
  mesh->owns_verts = true;

  fread((dbl3 *)mesh->verts, sizeof(dbl3), mesh->num_verts, fp);

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
  mesh->num_faces = size/sizeof(size_t[3]);
  fread((size_t(*)[3])mesh->faces, sizeof(size_t[3]), mesh->num_faces, fp);

  if (fclose(fp)) {
    perror("failed to close file");
    exit(EXIT_FAILURE);
  }
}

void mesh2_deinit(mesh2_s *mesh) {
  if (mesh->owns_verts) {
    free((dbl3 *)mesh->verts);
    mesh->num_verts = 0;
  }
  mesh->verts = NULL;
  mesh->owns_verts = false;

  free((size_t(*)[3])mesh->faces);
  mesh->faces = NULL;
  mesh->num_faces = 0;

  free((dbl3 *)mesh->face_normals);
  mesh->face_normals = NULL;
}

void mesh2_dump_verts(mesh2_s const *mesh, char const *path) {
  FILE *fp = fopen(path, "wb");
  fwrite(mesh->verts, sizeof(dbl3), mesh->num_verts, fp);
  fclose(fp);
}

void mesh2_dump_faces(mesh2_s const *mesh, char const *path) {
  FILE *fp = fopen(path, "wb");
  fwrite(mesh->faces, sizeof(size_t[3]), mesh->num_faces, fp);
  fclose(fp);
}

size_t mesh2_get_num_verts(mesh2_s const *mesh) {
  return mesh->num_verts;
}

dbl3 const *mesh2_get_verts_ptr(mesh2_s const *mesh) {
  return mesh->verts;
}

size_t mesh2_get_num_faces(mesh2_s const *mesh) {
  return mesh->num_faces;
}

uint3 const *mesh2_get_faces_ptr(mesh2_s const *mesh) {
  return mesh->faces;
}

rect3 mesh2_get_bounding_box(mesh2_s const *mesh) {
  rect3 bbox;

  bbox.min[0] = bbox.min[1] = bbox.min[2] = INFINITY;
  bbox.max[0] = bbox.max[1] = bbox.max[2] = -INFINITY;

  for (size_t i = 0; i < mesh->num_verts; ++i) {
    dbl x = mesh->verts[i][0];
    dbl y = mesh->verts[i][1];
    dbl z = mesh->verts[i][2];

    bbox.min[0] = fmin(x, bbox.min[0]);
    bbox.min[1] = fmin(y, bbox.min[1]);
    bbox.min[2] = fmin(z, bbox.min[2]);

    bbox.max[0] = fmax(x, bbox.max[0]);
    bbox.max[1] = fmax(y, bbox.max[1]);
    bbox.max[2] = fmax(z, bbox.max[2]);
  }

  return bbox;
}

void mesh2_get_centroid(mesh2_s const *mesh, size_t i, dbl *centroid) {
  assert(i < mesh->num_faces);
  size_t f[3];
  for (size_t j = 0; j < 3; ++j) {
    f[j] = mesh->faces[i][j];
  }
  for (size_t k = 0; k < 3; ++k) {
    centroid[k] = 0;
    for (size_t j = 0; j < 3; ++j) {
      centroid[k] += mesh->verts[f[j]][k];
    }
    centroid[k] /= 3;
  }
}

void mesh2_get_vertex(mesh2_s const *mesh, size_t i, size_t j, dbl *v) {
  size_t f = mesh->faces[i][j];
  memcpy((void *)v, (void *)&mesh->verts[f], sizeof(dbl3));
}

bool mesh2_tri_bbox_overlap(mesh2_s const *mesh, size_t i, rect3 const *bbox) {
  dbl boxhalfsize[3];
  rect3_get_extent(bbox, boxhalfsize);
  for (int i = 0; i < 3; ++i) {
    boxhalfsize[i] /= 2;
  }

  dbl boxcenter[3];
  for (int i = 0; i < 3; ++i) {
    boxcenter[i] = bbox->min[i] + boxhalfsize[i];
  }

  dbl triverts[3][3];
  for (size_t j = 0; j < 3; ++j) {
    memcpy(
      (void *)triverts[j],
      (void *)&mesh->verts[mesh->faces[i][j]],
      sizeof(dbl3));
  }

  return triBoxOverlap(boxcenter, boxhalfsize, triverts);
}

tri3 mesh2_get_tri(mesh2_s const *mesh, size_t i) {
  assert(i < mesh->num_faces);
  tri3 tri;
  memcpy(tri.v[0], &mesh->verts[mesh->faces[i][0]], sizeof(dbl3));
  memcpy(tri.v[1], &mesh->verts[mesh->faces[i][1]], sizeof(dbl3));
  memcpy(tri.v[2], &mesh->verts[mesh->faces[i][2]], sizeof(dbl3));
  return tri;
}

void mesh2_get_unit_surface_normal(mesh2_s const *mesh, size_t lf, dbl3 n) {
  dbl3_copy(mesh->face_normals[lf], n);
}
