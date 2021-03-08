#include "camera.h"

#include <assert.h>
#include <math.h>

#include "mat.h"
#include "vec.h"

void camera_init(camera_s *camera) {
  camera->type = CAMERA_TYPE_NONE;
  camera->pos[0] = camera->pos[1] = camera->pos[2] = NAN;
  camera->look[0] = camera->look[1] = camera->look[2] = NAN;
  camera->left[0] = camera->left[1] = camera->left[2] = NAN;
  camera->up[0] = camera->up[1] = camera->up[2] = NAN;
  camera->width = camera->height = NAN;
  camera->near = camera->fovy = camera->aspect = NAN;
  camera->dim[0] = camera->dim[1] = -1;
}

/**
 * Compute a camera ray with the direction determined by a pair of
 * image space coordinates:
 *
 *   |
 *  -+------->  +x[1]
 *   |
 *   |   * <- optical axis
 *   |
 *   |
 *   v
 *
 * +x[0]
 *
 * How this is done depends on, e.g., the camera type and other camera
 * parameters.
 */
ray3 camera_get_ray_for_index(camera_s const *camera, int i, int j) {
  assert(camera->type != CAMERA_TYPE_NONE);
  dbl h[2], mat[3][3] = {
    {camera->look[0], camera->left[0], camera->up[0]},
    {camera->look[1], camera->left[1], camera->up[1]},
    {camera->look[2], camera->left[2], camera->up[2]},
  };
  ray3 ray;
  if (camera->type == CAMERA_TYPE_ORTHOGRAPHIC) {
    assert(isfinite(camera->height));
    assert(isfinite(camera->width));
    h[0] = camera->height/(2*camera->dim[0]);
    h[1] = camera->width/(2*camera->dim[1]);
    dbl3_copy(camera->look, ray.dir);
    ray.org[0] = 0;
    ray.org[1] = camera->width*(0.5 - ((dbl)j)/camera->dim[1]) - h[1];
    ray.org[2] = camera->height*(0.5 - ((dbl)i)/camera->dim[0]) - h[0];
    dbl33_dbl3_mul_inplace(mat, ray.org);
    ray.org[0] += camera->pos[0];
    ray.org[1] += camera->pos[1];
    ray.org[2] += camera->pos[2];
  } else if (camera->type == CAMERA_TYPE_PERSPECTIVE) {
    assert(isfinite(camera->near));
    assert(isfinite(camera->fovy));
    assert(isfinite(camera->aspect));
    dbl height = 2*camera->near*tan(PI*(camera->fovy/2)/180);
    dbl width = camera->aspect*height;
    h[0] = height/(2*camera->dim[0]);
    h[1] = width/(2*camera->dim[1]);
    dbl3_copy(camera->pos, ray.org);
    ray.dir[0] = camera->near;
    ray.dir[1] = width/2 - j - h[1];
    ray.dir[2] = height/2 - i - h[0];
    dbl3_normalize(ray.dir);
    dbl33_dbl3_mul_inplace(mat, ray.dir);
    dbl3_copy(camera->pos, ray.org);
  } else {
    assert(false);
  }
  return ray;
}
