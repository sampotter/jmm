from defs cimport dbl

cdef extern from "jet.h":
    ctypedef struct jet3:
        dbl f
        dbl fx
        dbl fy
        dbl fz
    jet3 jet3_make_point_source(dbl tau)
