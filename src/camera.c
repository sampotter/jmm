#include "camera.h"

#include <assert.h>
#include <math.h>

#include "mat.h"
#include "vec.h"

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
ray3 camera_get_ray(camera_s const *camera, dbl x[2]) {
  dbl mat[3][3] = {
    {camera->look[0], camera->left[0], camera->up[0]},
    {camera->look[1], camera->left[1], camera->up[1]},
    {camera->look[2], camera->left[2], camera->up[2]},
  };
  ray3 ray;
  if (camera->type == CAMERA_TYPE_ORTHOGRAPHIC) {
    dbl h[2] = {camera->height/(2*camera->dim[0]),
                camera->width/(2*camera->dim[1])};
    dbl3_copy(camera->look, ray.dir);
    ray.org[0] = 0;
    ray.org[1] = camera->width/2 - x[1] - h[1];
    ray.org[2] = camera->height/2 - x[0] - h[0];
    // dbl33_transpose(mat);
    dbl33_dbl3_mul_inplace(mat, ray.org);
  } else if (camera->type == CAMERA_TYPE_PERSPECTIVE) {
    dbl3_copy(camera->pos, ray.org);
    ray.dir[0] = camera->near;
    ray.dir[2] = camera->near*tan(camera->fovy/2);
    ray.dir[1] = x[1] - camera->aspect*ray.dir[2];
    ray.dir[2] -= x[0];
    dbl33_dbl3_mul_inplace(mat, ray.dir);
    dbl3_normalize(ray.dir);
  } else {
    assert(false);
  }
  return ray;
}
