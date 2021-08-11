from jmm.defs cimport dbl, dbl2

cdef extern from "field.h":
    cdef struct field2:
        dbl (*f)(dbl, dbl, void *)
        void (*grad_f)(dbl, dbl, void *, dbl2)
        void *context

cdef class SlownessFunc2:
    cdef field2 field

cdef class LinearSpeedFunc2(SlownessFunc2):
    cdef:
        dbl c[3]
        dbl x0
        dbl y0
