from jmm.defs cimport dbl, dbl2, dbl22, dbl44

cdef extern from "bicubic.h":
    cdef struct bicubic:
        dbl44 A

    dbl bicubic_f(const bicubic *bicubic, dbl2 cc)
    dbl bicubic_fx(const bicubic *bicubic, dbl2 cc)
    dbl bicubic_fy(const bicubic *bicubic, dbl2 cc)
    dbl bicubic_fxx(const bicubic *bicubic, dbl2 cc)
    dbl bicubic_fxy(const bicubic *bicubic, dbl2 cc)
    dbl bicubic_fyy(const bicubic *bicubic, dbl2 cc)

cdef class Bicubic:
    cdef bicubic bicubic

    @staticmethod
    cdef from_bicubic(bicubic bicubic)
