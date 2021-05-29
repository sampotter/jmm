from jmm.defs cimport dbl

cdef extern from "slerp.h":
    void slerp2(const dbl p0[3], const dbl p1[3], const dbl w[2], dbl q[3])
    void slerp3(const dbl p[3][3], const dbl w[3], const dbl q[3], dbl tol)
    void slerp4(const dbl p[4][3], const dbl w[4], const dbl q[3], dbl tol)
