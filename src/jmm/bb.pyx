import cython
import numpy as np

from libc.string cimport memcpy

cimport jmm.jet

@cython.auto_pickle(True)
cdef class Bb31:
    @staticmethod
    def from_1d_data(dbl[::1] f, dbl[::1] Df, dbl[::1] x):
        if f.size != 2 or f.shape[0] != 2:
            raise ValueError('`f` must be a length 2 vector')
        if Df.size != 2 or Df.shape[0] != 2:
            raise ValueError('`Df` must be a length 2 vector')
        if x.size != 2 or x.shape[0] != 2:
            raise ValueError('`x` must be a length 2 vector')
        bb = Bb31()
        bb31_init_from_1d_data(&bb.bb, &f[0], &Df[0], &x[0])
        return bb

    @property
    def c(self):
        cdef dbl[4] c = np.empty((4,), dtype=np.float64)
        memcpy(&c[0], self.bb.c, 4*sizeof(dbl))
        return np.asarray(c)

    def isfinite(self):
        return all(np.isfinite(c) for c in self.bb.c)

    @staticmethod
    cdef from_bb31(bb31 bb):
        bb_ = Bb31()
        bb_.bb = bb
        return bb_

cdef class Bb33:
    @staticmethod
    def from_3d_data(dbl[:] f, dbl[:, :] Df, dbl[:, :] x):
        if f.size != 4 or f.shape[0] != 4:
            raise ValueError('`f` must be a length 4 vector')
        if Df.size != 12 or Df.shape[0] != 4 or Df.shape[1] != 3:
            raise ValueError('`Df` must have shape (4, 3)')
        if x.size != 12 or x.shape[0] != 4 or x.shape[1] != 3:
            raise ValueError('`x` must have shape (4, 3)')
        bb = Bb33()
        bb33_init_from_3d_data(
            &bb._bb,
            &f[0],
            <const dbl (*)[3]>&Df[0, 0],
            <const dbl (*)[3]>&x[0, 0])
        return bb

    @staticmethod
    def from_jets(dbl[:, ::1] x, jets):
        if x.size != 12 or x.shape[0] != 4 or x.shape[1] != 3:
            raise ValueError('`x` must have shape (4, 3)')
        if len(jets) != 4:
            raise ValueError('`jets` must satisfy `len(jets) == 4`')
        cdef jet3[4] jets_
        cdef size_t i
        for i in range(4):
            if isinstance(jets[i], np.void):
                jets_[i] = jmm.jet.Jet3(*jets[i]).jet
            elif isinstance(jets[i], jmm.jet.Jet3):
                jets_[i] = jets[i].jet
            else:
                raise ValueError(
                    '`jets[%d]` contains unexpected type: %s' % (
                        i, type(jets[i]).__name__))
        bb = Bb33()
        bb33_init_from_jets(&bb._bb, &jets_[0], <const dbl (*)[3]>&x[0, 0])
        return bb

    @staticmethod
    cdef from_bb33(bb33 bb):
        bb_ = Bb33()
        bb_._bb = bb
        return bb_

    def f(self, dbl[:] b):
        if b.size != 4 or b.shape[0] != 4:
            raise ValueError('`b` must be a length 4 vector')
        return bb33_f(&self._bb, &b[0])

    def Df(self, dbl[:] b, dbl[:] a):
        if b.size != 4 or b.shape[0] != 4:
            raise ValueError('`b` must be a length 4 vector')
        if a.size != 4 or a.shape[0] != 4:
            raise ValueError('`a` must be a length 4 vector')
        return bb33_df(&self._bb, &b[0], &a[0])

    def D2f(self, dbl[:] b, dbl[:, :] a):
        if b.size != 4 or b.shape[0] != 4:
            raise ValueError('`b` must be a length 4 vector')
        if a.size != 12 or a.shape[0] != 3 or a.shape[4] != 4:
            raise ValueError('`a` must have shape (3, 4)')
        return bb33_d2f(&self._bb, &b[0], <const dbl (*)[4]>&a[0, 0])

    def convex_hull_brackets_value(self, dbl value):
        return bb33_convex_hull_brackets_value(&self._bb, value)
