from defs cimport dbl
from geom cimport ray3

cdef extern from "camera.h":
    cdef enum camera_type:
        CAMERA_TYPE_NONE
        CAMERA_TYPE_ORTHOGRAPHIC
        CAMERA_TYPE_PERSPECTIVE
    cdef struct camera:
        camera_type type
        dbl pos[3]
        dbl look[3]
        dbl left[3]
        dbl up[3]
        dbl width
        dbl height
        dbl near
        dbl fovy
        dbl aspect
        size_t dim[2]
    void camera_init(camera *camera)
    ray3 camera_get_ray_for_index(const camera *camera, int i, int j)
