ctypedef bint bool

cdef extern from "def.h":
    ctypedef double dbl

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
