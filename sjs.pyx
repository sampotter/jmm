# distutils: language=c99
# cython: embedsignature=True
# cython: language_level=3


cdef extern from "def.h":
    ctypedef double dbl


cdef extern from "vec.h":
    cdef struct dvec2:
        dbl x
        dbl y
    cdef struct ivec2:
        int i
        int j

cdef extern from "sjs.h":
    ctypedef dbl (*sfield)(dvec2)
    ctypedef dvec2 (*vfield)(dvec2)

    cdef struct sjs:
        pass

    void sjs_alloc(sjs **sjs)
    void sjs_init(sjs *sjs, ivec2 shape, dvec2 xymin, dbl h, sfield s,
                  vfield grad_s)
    void sjs_dealloc(sjs **sjs)


cdef class Sjs:
    cdef:
        sjs *sjs

    def __cinit__(self, shape, xymin, dbl h, sfield s, vfield grad_s):
        sjs_alloc(&self.sjs)
        cdef ivec2 shape_ = {shape[0], shape[1]}
        cdef dvec2 xymin_ = {xymin[0], xymin[1]}
        sjs_init(self.sjs, shape_, xymin_, h, s, grad_s)
