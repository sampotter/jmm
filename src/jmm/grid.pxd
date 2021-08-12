from jmm.defs cimport dbl, dbl2, int2, order, ORDER_ROW_MAJOR

cdef extern from "grid2.h":
    cdef struct grid2:
        int2 shape
        dbl2 xymin
        dbl h
        order order

cdef class Grid2:
    cdef grid2 _grid

cdef extern from "grid3.h":
    cdef struct grid3:
        int dim[3]
        dbl min[3]
        dbl h

cdef class Grid3:
    cdef grid3 _grid
