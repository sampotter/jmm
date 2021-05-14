from jmm.defs cimport dbl

cdef extern from "geom.h":
    ctypedef struct rect3:
        dbl min[3]
        dbl max[3]

    ctypedef struct tri3:
        dbl v[3][3]

    ctypedef struct ray3:
        dbl org[3]
        dbl dir[3]

    void ray3_get_point(const ray3 *ray, dbl t, dbl x[3])

cdef class Rect3:
    cdef rect3 _rect

cdef class Ray3:
    cdef ray3 _ray

    @staticmethod
    cdef from_ptr(const ray3 *ray)
