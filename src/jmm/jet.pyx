cdef class Jet2:
    def __cinit__(self, dbl f, dbl fx, dbl fy, dbl fxy):
        self.jet.f = f
        self.jet.Df[0] = fx
        self.jet.Df[1] = fy
        self.jet.fxy = fxy

    @property
    def f(self):
        return self.jet.f

    @property
    def fx(self):
        return self.jet.Df[0]

    @property
    def fy(self):
        return self.jet.Df[1]

    @property
    def fxy(self):
        return self.jet.fxy

    def __repr__(self):
        return 'Jet2(f = %f, fx = %f, fy = %f, fxy = %f)' % (
            self.jet.f, self.jet.Df[0], self.jet.Df[1], self.jet.fxy)

cdef class Jet22t:
    def __cinit__(self, dbl f, dbl fx, dbl fy, dbl fxx, dbl fyx, dbl fxy, dbl fyy):
        self.jet.f = f
        self.jet.Df[0] = fx
        self.jet.Df[1] = fy
        self.jet.D2f[0][0] = fxx
        self.jet.D2f[0][1] = fyx
        self.jet.D2f[1][0] = fxy
        self.jet.D2f[1][1] = fyy

    @property
    def f(self):
        return self.jet.f

    @property
    def fx(self):
        return self.jet.Df[0]

    @property
    def fy(self):
        return self.jet.Df[1]

    @property
    def fxx(self):
        return self.jet.D2f[0][0]

    @property
    def fyx(self):
        return self.jet.D2f[0][1]

    @property
    def fxy(self):
        return self.jet.D2f[1][0]

    @property
    def fyy(self):
        return self.jet.D2f[1][1]

    def __repr__(self):
        return 'Jet2(f = %f, fx = %f, fy = %f, fxx = %f, fyx = %f, fxy = %f, fyy = %f)' % (
            self.jet.f,
            self.jet.Df[0], self.jet.Df[1],
            self.jet.D2f[0][0], self.jet.D2f[0][1], self.jet.D2f[1][0], self.jet.D2f[1][1])

cdef class Jet3:
    def __cinit__(self, dbl f, dbl fx, dbl fy, dbl fz):
        self.jet.f = f
        self.jet.Df[0] = fx
        self.jet.Df[1] = fy
        self.jet.Df[2] = fz

    @staticmethod
    def make_point_source(dbl tau=0):
        cdef jet3 jet = jet3_make_point_source(tau)
        return Jet3(jet.f, jet.Df[0], jet.Df[1], jet.Df[2])

    def __repr__(self):
        return 'Jet3(f = %f, fx = %f, fy = %f, fz = %f)' % (
            self.jet.f, self.jet.Df[0], self.jet.Df[1], self.jet.Df[2])

    @property
    def f(self):
        return self.jet.f

    @property
    def fx(self):
        return self.jet.Df[0]

    @property
    def fy(self):
        return self.jet.Df[1]

    @property
    def fz(self):
        return self.jet.Df[2]

    def is_finite(self):
        return jet3_is_finite(&self.jet)
