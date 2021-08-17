cdef class Bicubic:
    @staticmethod
    cdef from_bicubic(bicubic bicubic):
        bicubic_ = Bicubic()
        bicubic_.bicubic = bicubic
        return bicubic_

    def f(self, dbl x, dbl y):
        cdef dbl[2] cc = [x, y]
        return bicubic_f(&self.bicubic, cc)

    def fx(self, dbl x, dbl y):
        cdef dbl[2] cc = [x, y]
        return bicubic_fx(&self.bicubic, cc)

    def fy(self, dbl x, dbl y):
        cdef dbl[2] cc = [x, y]
        return bicubic_fy(&self.bicubic, cc)

    def fxx(self, dbl x, dbl y):
        cdef dbl[2] cc = [x, y]
        return bicubic_fxx(&self.bicubic, cc)

    def fxy(self, dbl x, dbl y):
        cdef dbl[2] cc = [x, y]
        return bicubic_fxy(&self.bicubic, cc)

    def fyy(self, dbl x, dbl y):
        cdef dbl[2] cc = [x, y]
        return bicubic_fyy(&self.bicubic, cc)
