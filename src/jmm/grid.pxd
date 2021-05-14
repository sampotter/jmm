from jmm.defs cimport dbl

cdef extern from "grid3.h":
    cdef struct grid3:
        int dim[3]
        dbl min[3]
        dbl h

cdef class Grid3:
    cdef grid3 _grid
