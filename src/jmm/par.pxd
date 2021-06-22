from jmm.defs cimport dbl

cdef extern from "par.h":
    cdef struct par3:
        size_t l[3]
        dbl b[3]
    size_t par3_size(const par3 *par)
    size_t par3_num_active(const par3 *par)
    void par3_get_active(const par3 *par, size_t *l, dbl *b)

cdef class Parent3:
    cdef par3 par
