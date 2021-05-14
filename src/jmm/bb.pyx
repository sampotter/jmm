cdef class Bb33:
    @staticmethod
    def from_3d_data(dbl[:] f, dbl[:, :] Df, dbl[:, :] x):
        if f.size != 4 or f.shape[0] != 4:
            raise Exception('`f` must be a length 4 vector')
        if Df.size != 12 or Df.shape[0] != 4 or Df.shape[1] != 3:
            raise Exception('`Df` must have shape (4, 3)')
        if x.size != 12 or x.shape[0] != 4 or x.shape[1] != 3:
            raise Exception('`x` must have shape (4, 3)')
        bb = Bb33()
        bb33_init_from_3d_data(
            &bb._bb,
            &f[0],
            <const dbl (*)[3]>&Df[0, 0],
            <const dbl (*)[3]>&x[0, 0])
        return bb

    @staticmethod
    cdef from_bb33(bb33 bb):
        bb_ = Bb33()
        bb_._bb = bb
        return bb_

    def f(self, dbl[:] b):
        if b.size != 4 or b.shape[0] != 4:
            raise Exception('`b` must be a length 4 vector')
        return bb33_f(&self._bb, &b[0])

    def Df(self, dbl[:] b, dbl[:] a):
        if b.size != 4 or b.shape[0] != 4:
            raise Exception('`b` must be a length 4 vector')
        if a.size != 4 or a.shape[0] != 4:
            raise Exception('`a` must be a length 4 vector')
        return bb33_df(&self._bb, &b[0], &a[0])

    def D2f(self, dbl[:] b, dbl[:, :] a):
        if b.size != 4 or b.shape[0] != 4:
            raise Exception('`b` must be a length 4 vector')
        if a.size != 12 or a.shape[0] != 3 or a.shape[4] != 4:
            raise Exception('`a` must have shape (3, 4)')
        return bb33_d2f(&self._bb, &b[0], <const dbl (*)[4]>&a[0, 0])

    def convex_hull_brackets_value(self, dbl value):
        return bb33_convex_hull_brackets_value(&self._bb, value)
