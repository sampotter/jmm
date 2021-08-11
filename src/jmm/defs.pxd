cimport cython

ctypedef bint bool

cdef extern from "def.h":
    ctypedef double dbl

    ctypedef dbl[2] dbl2

    ctypedef dbl[3] dbl3
    ctypedef dbl3[3] dbl33

    ctypedef int[2] int2

    cdef enum state:
        FAR
        TRIAL
        VALID
        BOUNDARY
        ADJACENT_TO_BOUNDARY
        NEW_VALID
        SHADOW

    cdef enum ftype:
        FTYPE_POINT_SOURCE
        FTYPE_REFLECTION
        FTYPE_EDGE_DIFFRACTION

    cdef enum stype:
        CONSTANT
        NUM_STYPE

    cdef enum error:
        SUCCESS
        BAD_ARGUMENT

ctypedef double complex dblz

ctypedef fused dbl_or_dblz:
    dbl
    dblz
