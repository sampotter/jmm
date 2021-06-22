from jmm.array_view cimport ArrayView
from jmm.defs cimport bool, dbl, ftype, state
from jmm.jet cimport jet3
from jmm.par cimport par3
from jmm.mesh cimport mesh3

cdef extern from "eik3.h":
    cdef struct eik3:
        pass

    void eik3_alloc(eik3 **eik)
    void eik3_dealloc(eik3 **eik)
    void eik3_init(eik3 *eik, const mesh3 *mesh, ftype ftype)
    void eik3_deinit(eik3 *eik)
    size_t eik3_peek(const eik3 *eik)
    size_t eik3_step(eik3 *eik)
    void eik3_solve(eik3 *eik)
    const mesh3 *eik3_get_mesh(const eik3 *eik)
    bool eik3_is_far(const eik3 *eik, size_t ind)
    bool eik3_is_trial(const eik3 *eik, size_t ind)
    bool eik3_is_valid(const eik3 *eik, size_t ind)
    jet3 *eik3_get_jet_ptr(const eik3 *eik)
    state *eik3_get_state_ptr(const eik3 *eik)
    par3 eik3_get_par(const eik3 *eik, size_t l)
    dbl *eik3_get_t_in_ptr(const eik3 *eik)
    dbl *eik3_get_t_out_ptr(const eik3 *eik)
    void eik3_add_pt_src_BCs(eik3 *eik, size_t l, jet3 jet)
    void eik3_add_refl_BCs(eik3 *eik, const size_t lf[3], const jet3 jet[3],
                           const dbl t_in[3][3])
    void eik3_add_diff_edge_BCs(eik3 *eik, const size_t le[2], const jet3 jet[2])
    ftype eik3_get_ftype(const eik3 *eik)
    dbl eik3_get_slerp_tol(const eik3 *eik)
    bool eik3_has_BCs(const eik3 *eik, size_t l)
    void eik3_transport_dbl(const eik3 *eik, dbl *values, bool skip_filled)
    void eik3_transport_dblz(const eik3 *eik, dblz *values, bool skip_filled)
    void eik3_transport_curvature(const eik3 *eik, dbl *kappa, bool skip_filled)
    dbl eik3_get_h(const eik3* eik)

cdef class Eik3:
    cdef:
        eik3 *eik
        ArrayView jet_view
        ArrayView T_view
        ArrayView grad_T_view
        ArrayView state_view
        ArrayView t_in_view
        ArrayView t_out_view
