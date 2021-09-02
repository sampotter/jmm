from jmm.defs cimport bool, dbl

cdef extern from "jet.h":
    ctypedef struct jet2:
        dbl f
        dbl fx
        dbl fy
        dbl fxy

    ctypedef struct jet22t:
        dbl f
        dbl fx
        dbl fy
        dbl fxx
        dbl fyx
        dbl fxy
        dbl fyy

    ctypedef struct jet3:
        dbl f
        dbl fx
        dbl fy
        dbl fz
    jet3 jet3_make_point_source(dbl tau)
    bool jet3_is_finite(const jet3 *jet)

cdef class Jet2:
    cdef jet2 jet

cdef class Jet22t:
    cdef jet22t jet

cdef class Jet3:
    cdef jet3 jet
