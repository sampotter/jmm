#include "mesh2.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "def.h"
#include "log.h"

struct mesh2 {
  dbl *points;
  size_t num_points;
  size_t *faces;
  size_t *num_faces;
};

void mesh2_init_from_binary_files(mesh2_s *mesh, char const *verts_path,
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
  log_debug("%llu", size/(3*sizeof(dbl)));

  mesh->points = malloc(size);
  mesh->num_points = size/(3*sizeof(dbl));
  fread(mesh->points, 3*sizeof(dbl), mesh->num_points, fp);

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
  mesh->num_faces = size/(3*sizeof(size_t));
  fread(mesh->faces, 3*sizeof(size_t), mesh->num_faces, fp);

  if (fclose(fp)) {
    perror("failed to close file");
    exit(EXIT_FAILURE);
  }
}

void mesh2_deinit(mesh2_s *mesh) {
  free(mesh->points);
  mesh->points = NULL;
  mesh->num_points = 0;

  free(mesh->faces);
  mesh->faces = NULL;
  mesh->num_faces = 0;
}

size_t mesh2_get_num_points(mesh2_s const *mesh) {
  return mesh->num_points;
}

dbl *mesh2_get_points_ptr(mesh2_s const *mesh) {
  return mesh->points;
}

size_t mesh2_get_num_faces(mesh2_s const *mesh) {
  return mesh->num_faces;
}

size_t *mesh2_get_faces_ptr(mesh2_s const *mesh) {
  return mesh->faces;
}

rect3 mesh2_get_bounding_box(mesh2_s const *mesh) {
  rect3 bbox;

  bbox.min[0] = bbox.min[1] = bbox.min[2] = INFINITY;
  bbox.max[0] = bbox.max[1] = bbox.max[2] = -INFINITY;

  for (size_t i = 0; i < mesh->num_points; ++i) {
    dbl x = mesh->points[3*i];
    dbl y = mesh->points[3*i + 1];
    dbl z = mesh->points[3*i + 2];

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
    f[j] = mesh->faces[3*i + j];
  }
  for (size_t k = 0; k < 3; ++k) {
    centroid[k] = 0;
    for (size_t j = 0; j < 3; ++j) {
      centroid[k] += mesh->points[3*f[j] + k];
    }
    centroid[k] /= 3;
  }
}

void mesh2_get_vertex(mesh2_s const *mesh, size_t i, size_t j, dbl *v) {
  size_t f = mesh->faces[3*i + j];
  memcpy((void *)v, (void *)&mesh->points[3*f], 3*sizeof(dbl));
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
      (void *)triverts[j], (void *)&mesh->points[3*mesh->faces[3*i + j]],
      3*sizeof(dbl));
  }

  return triBoxOverlap(boxcenter, boxhalfsize, triverts);
}

tri3 mesh2_get_tri(mesh2_s const *mesh, size_t i) {
  tri3 tri;
  memcpy(tri.v[0], &mesh->points[3*mesh->faces[i]], sizeof(dbl[3]));
  memcpy(tri.v[1], &mesh->points[3*mesh->faces[i] + 1], sizeof(dbl[3]));
  memcpy(tri.v[2], &mesh->points[3*mesh->faces[i] + 2], sizeof(dbl[3]));
  return tri;
}
