from jmm.defs cimport bool, dbl, state
from jmm.eik cimport eik3
from jmm.jet cimport jet3

cdef extern from "utri.h":
    cdef struct utri_spec:
        const eik3 *eik
        size_t lhat
        size_t l[2]
        state state[2]
        dbl xhat[3]
        dbl x[2][3]
        jet3 jet
        size_t orig_index

    utri_spec utri_spec_from_eik_without_l(const eik3 *eik, const dbl *x,
                                           size_t l0, size_t l1)

    cdef struct utri:
        pass

    void utri_alloc(utri **u)
    void utri_dealloc(utri **u)
    bool utri_init(utri *u, const utri_spec *spec)
    void utri_solve(utri *u)
    bool utri_approx_hess(const utri *u, dbl h, dbl hess[3][3])

cdef class Utri:
    cdef:
        utri *u
        dbl h