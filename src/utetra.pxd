from defs cimport bool, dbl
from eik3 cimport *

cdef extern from "utetra.h":
    cdef struct utetra:
        pass
    void utetra_alloc(utetra **cf)
    void utetra_dealloc(utetra **cf)
    void utetra_init_from_eik3(utetra *cf, const eik3 *eik,
                              size_t l, size_t l0, size_t l1, size_t l2)
    void utetra_init(utetra *cf, const dbl x[3], const dbl Xt[3][3],
                     const jet3 jet[3])
    void utetra_deinit(utetra *cf)
    bool utetra_is_degenerate(const utetra *cf)
    bool utetra_is_causal(const utetra *cf)
    void utetra_solve(utetra *cf, const dbl *lam)
    dbl utetra_get_value(const utetra *cf)
    void utetra_get_jet(const utetra *cf, jet3 *jet)
