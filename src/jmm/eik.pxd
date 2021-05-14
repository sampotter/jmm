from jmm.array_view cimport ArrayView
from jmm.defs cimport bool, dbl, state
from jmm.jet cimport jet3
from jmm.par cimport par3
from jmm.mesh cimport mesh3

cdef extern from "eik3.h":
    cdef struct cutedge:
        dbl t
        dbl n[3]
        jet3 jet

    cdef struct eik3:
        pass

    void eik3_alloc(eik3 **eik)
    void eik3_dealloc(eik3 **eik)
    void eik3_init(eik3 *eik, const mesh3 *mesh)
    void eik3_deinit(eik3 *eik)
    size_t eik3_peek(const eik3 *eik)
    size_t eik3_step(eik3 *eik)
    void eik3_solve(eik3 *eik)
    void eik3_add_trial(eik3 *eik, size_t ind, jet3 jet)
    const mesh3 *eik3_get_mesh(const eik3 *eik)
    bool eik3_is_far(const eik3 *eik, size_t ind)
    bool eik3_is_trial(const eik3 *eik, size_t ind)
    bool eik3_is_valid(const eik3 *eik, size_t ind)
    jet3 *eik3_get_jet_ptr(const eik3 *eik)
    state *eik3_get_state_ptr(const eik3 *eik)
    par3 eik3_get_par(const eik3 *eik, size_t l)
    dbl *eik3_get_t_in_ptr(const eik3 *eik)
    dbl *eik3_get_t_out_ptr(const eik3 *eik)
    void eik3_add_valid_bdf(eik3 *eik, const size_t lf[3], const jet3 jet[3])
    void eik3_add_valid_bde(eik3 *eik, const size_t le[2], const jet3 jet[2])

cdef class Eik3:
    cdef:
        eik3 *eik
        ArrayView jet_view
        ArrayView state_view
        ArrayView t_in_view
        ArrayView t_out_view
