from defs cimport dbl

cdef extern from "bb.h":
    void bb3tet_interp3(const dbl f[4], const dbl Df[4][3], const dbl x[4][3], dbl c[20])
    dbl bb3tet(const dbl c[20], const dbl b[4])
    dbl dbb3tet(const dbl c[20], const dbl b[4], const dbl a[4])
    dbl d2bb3tet(const dbl c[20], const dbl b[4], const dbl a[2][4])
