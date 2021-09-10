from jmm.defs cimport bool, dbl, dbl2, dbl22, dbl3

cdef extern from "jet.h":
    ctypedef struct jet2:
        dbl f
        dbl2 Df
        dbl fxy

    ctypedef struct jet22t:
        dbl f
        dbl2 Df
        dbl22 D2f

    ctypedef struct jet3:
        dbl f
        dbl3 Df

    jet3 jet3_make_point_source(dbl tau)
    bool jet3_is_finite(const jet3 *jet)

cdef class Jet2:
    cdef jet2 jet

cdef class Jet22t:
    cdef jet22t jet

cdef class Jet3:
    cdef jet3 jet
