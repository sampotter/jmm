cdef class Jet3:
    def __cinit__(self, dbl f, dbl fx, dbl fy, dbl fz):
        self.jet.f = f
        self.jet.fx = fx
        self.jet.fy = fy
        self.jet.fz = fz

    @staticmethod
    def make_point_source(dbl tau=0):
        cdef jet3 jet = jet3_make_point_source(tau)
        return Jet3(jet.f, jet.fx, jet.fy, jet.fz)

    def __repr__(self):
        return 'Jet3(f = %f, fx = %f, fy = %f, fz = %f)' % (
            self.jet.f, self.jet.fx, self.jet.fy, self.jet.fz)

    @property
    def f(self):
        return self.jet.f

    @property
    def fx(self):
        return self.jet.fx

    @property
    def fy(self):
        return self.jet.fy

    @property
    def fz(self):
        return self.jet.fz

    def is_finite(self):
        return jet3_is_finite(&self.jet)
