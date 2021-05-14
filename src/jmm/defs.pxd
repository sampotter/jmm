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
    cdef enum stype:
        CONSTANT
        NUM_STYPE
    cdef enum error:
        SUCCESS
        BAD_ARGUMENT
