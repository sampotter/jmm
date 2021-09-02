cdef class Jet2:
    def __cinit__(self, dbl f, dbl fx, dbl fy, dbl fxy):
        self.jet.f = f
        self.jet.fx = fx
        self.jet.fy = fy
        self.jet.fxy = fxy

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
    def fxy(self):
        return self.jet.fxy

    def __repr__(self):
        return 'Jet2(f = %f, fx = %f, fy = %f, fxy = %f)' % (
            self.jet.f, self.jet.fx, self.jet.fy, self.jet.fxy)

cdef class Jet22t:
    def __cinit__(self, dbl f, dbl fx, dbl fy, dbl fxx, dbl fyx, dbl fxy, dbl fyy):
        self.jet.f = f
        self.jet.fx = fx
        self.jet.fy = fy
        self.jet.fxx = fxx
        self.jet.fyx = fyx
        self.jet.fxy = fxy
        self.jet.fyy = fyy

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
    def fxx(self):
        return self.jet.fxx

    @property
    def fyx(self):
        return self.jet.fyx

    @property
    def fxy(self):
        return self.jet.fxy

    @property
    def fyy(self):
        return self.jet.fyy

    def __repr__(self):
        return 'Jet2(f = %f, fx = %f, fy = %f, fxx = %f, fyx = %f, fxy = %f, fyy = %f)' % (
            self.jet.f,
            self.jet.fx, self.jet.fy,
            self.jet.fxx, self.jet.fyx, self.jet.fxy, self.jet.fyy)

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
