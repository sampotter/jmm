import numpy as np

from jmm.eik cimport Eik3

cdef class Utri:
    def __init__(self, *args):
        if len(args) > 0:
            raise RuntimeError('construct Utri using factory functions')

    def __cinit__(self, *args):
        if len(args) > 0:
            raise RuntimeError('construct Utri using factory functions')

    @staticmethod
    def from_eik_without_l(Eik3 eik, dbl[::1] x, size_t[::1] l):
        cdef utri_spec spec = \
            utri_spec_from_eik_without_l(eik.eik, &x[0], l[0], l[1])

        u = Utri()
        utri_alloc(&u.u)
        utri_init(u.u, &spec)
        u.h = eik.h
        return u

    def __dealloc__(self):
        utri_dealloc(&self.u)

    def solve(self):
        utri_solve(self.u)

    def approx_hess(self):
        cdef dbl[:, ::1] hess = np.empty((3, 3), dtype=np.float64)
        utri_approx_hess(self.u, self.h, <dbl (*)[3]>&hess[0, 0])
        return np.asarray(hess)
