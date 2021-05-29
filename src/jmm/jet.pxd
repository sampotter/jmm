from jmm.defs cimport bool, dbl

cdef extern from "jet.h":
    ctypedef struct jet3:
        dbl f
        dbl fx
        dbl fy
        dbl fz
    jet3 jet3_make_point_source(dbl tau)
    bool jet3_is_finite(const jet3 *jet)

cdef class Jet3:
    cdef jet3 jet
