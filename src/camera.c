#include "camera.h"

#include <assert.h>
#include <math.h>

#include "macros.h"
#include "mat.h"
#include "vec.h"

void camera_reset(camera_s *camera) {
  camera->type = CAMERA_TYPE_UNINITIALIZED;
  camera->pos[0] = camera->pos[1] = camera->pos[2] = NAN;
  camera->look[0] = camera->look[1] = camera->look[2] = NAN;
  camera->left[0] = camera->left[1] = camera->left[2] = NAN;
  camera->up[0] = camera->up[1] = camera->up[2] = NAN;
  camera->width = camera->height = NAN;
  camera->fovy = camera->aspect = NAN;
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
  assert(camera->type != CAMERA_TYPE_UNINITIALIZED);
  dbl h[2], mat[3][3] = {
    {camera->look[0], camera->left[0], camera->up[0]},
    {camera->look[1], camera->left[1], camera->up[1]},
    {camera->look[2], camera->left[2], camera->up[2]},
  };
  ray3 ray = ray3_make_empty();
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
    assert(isfinite(camera->fovy));
    assert(isfinite(camera->aspect));

    dbl phi_max = (JMM_PI/180)*camera->fovy/2;
    dbl h = 2*tan(phi_max);
    dbl w = camera->aspect*h;
    // dbl theta_max = atan(w/2);

    // dbl dphi = 2*phi_max/(camera->dim[0] - 1);
    // dbl dtheta = 2*theta_max/(camera->dim[1] - 1);

    // TODO: wrong way to do it... gives a fisheye effect
    dbl y = h/2 - i*h/(camera->dim[0] - 1);
    dbl x = w/2 - j*w/(camera->dim[1] - 1);
    dbl phi = atan(y);
    dbl theta = atan(x);

    // dbl phi = phi_max - i*dphi;
    // dbl theta = theta_max - j*dtheta;

    dbl3_copy(camera->pos, ray.org);

    ray.dir[0] = cos(theta)*cos(phi);
    ray.dir[1] = sin(theta)*cos(phi);
    ray.dir[2] = sin(phi);

    dbl33_dbl3_mul_inplace(mat, ray.dir);
  } else {
    die();
  }
  return ray;
}
