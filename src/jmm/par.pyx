import numpy as np

from libc.string cimport memcpy

cdef class Parent3:
    def __cinit__(self, par3 par):
        self.par = par

    @property
    def size(self):
        return par3_size(&self.par)

    @property
    def l(self):
        cdef size_t[:] l = np.empty((self.size,), dtype=np.uintp)
        memcpy(&l[0], self.par.l, self.size*sizeof(size_t))
        return np.asarray(l)

    @property
    def b(self):
        cdef dbl[:] b = np.empty((self.size,), dtype=np.float64)
        memcpy(&b[0], self.par.b, self.size*sizeof(dbl))
        return np.asarray(b)
