import numpy as np

from libc.string cimport memcpy

from jmm.defs cimport dbl

cdef class Rect3:
    def __cinit__(self, rect3 rect):
        self._rect = rect

    def __repr__(self):
        return f'Rect3(xmin={self._rect.min[0]}, xmax={self._rect.max[0]}, zmin={self._rect.min[1]}, zmax={self._rect.max[1]}, zmin={self._rect.min[2]}, zmax={self._rect.max[2]})'

    def __str__(self):
        return f'[{self._rect.min[0]}, {self._rect.max[0]}] x [{self._rect.min[1]}, {self._rect.max[1]}] x [{self._rect.min[2]}, {self._rect.max[2]}]'

cdef class Ray3:
    '''A ray in 3D, consisting of an `origin` and `direction`. The
`direction` vector doesn't need to be parametrized, but its magnitude
_does_ determine the parametrization of the ray. That is, for a
parameter `t`, a point `p(t)` on the ray equals `p(t) = origin +
direction*t`.

    :param origin: A 3D vector indicating the start of the ray,
        typically an instance of :class:`numpy.ndarray`.
    :param direction: The direction of ray. See comment above about
        normalization. Also typically an :class:`numpy.ndarray`
        instance.
    '''

    def __cinit__(self, dbl[::1] origin, dbl[::1] direction):
        memcpy(&self._ray.org[0], &origin[0], 3*sizeof(dbl))
        memcpy(&self._ray.dir[0], &direction[0], 3*sizeof(dbl))

    @staticmethod
    cdef from_ptr(const ray3 *ray):
        cdef dbl[::1] org = np.array([ray.org[0], ray.org[1], ray.org[2]])
        cdef dbl[::1] dir = np.array([ray.dir[0], ray.dir[1], ray.dir[2]])
        return Ray3(org, dir)

    def __repr__(self):
        return f'Ray3(origin={self._ray.org}, direction={self._ray.dir})'

    def get_point(self, dbl t):
        x = np.empty((3,), dtype=np.float64)
        cdef dbl[::1] x_mv = x
        ray3_get_point(&self._ray, t, &x_mv[0])
        return x
