from defs cimport dbl

cdef extern from "jet.h":
    ctypedef struct jet3:
        dbl f
        dbl fx
        dbl fy
        dbl fz
