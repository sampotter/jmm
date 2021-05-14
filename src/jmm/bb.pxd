from jmm.defs cimport bool, dbl

cdef extern from "bb.h":
    ctypedef struct bb33:
        dbl c[20]
    void bb33_init_from_3d_data(bb33 *bb, const dbl f[4], const dbl Df[4][3], const dbl x[4][3])
    dbl bb33_f(const bb33 *bb, const dbl b[4])
    dbl bb33_df(const bb33 *bb, const dbl b[4], const dbl a[4])
    dbl bb33_d2f(const bb33 *bb, const dbl b[4], const dbl a[2][4])
    bool bb33_convex_hull_brackets_value(const bb33 *bb, dbl value)

cdef class Bb33:
    cdef bb33 _bb

    @staticmethod
    cdef from_bb33(bb33 bb)
