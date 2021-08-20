from jmm.array_view cimport ArrayView
from jmm.bb cimport bb31
from jmm.bicubic cimport bicubic
from jmm.defs cimport bool, dbl, dblz, dbl_or_dblz, dbl2, dbl3, dbl33, ftype, \
    int2, state
from jmm.field cimport field2
from jmm.grid cimport grid2
from jmm.jet cimport jet2, jet3
from jmm.par cimport par2, par3
from jmm.mesh cimport mesh3

cdef extern from "eik.h":
    cdef struct eik:
        pass

    void eik_alloc(eik **eik)
    void eik_dealloc(eik **eik)
    void eik_init(eik *eik, const field2 *slow, const grid2 *grid)
    void eik_deinit(eik *eik)
    size_t eik_peek(const eik *eik)
    void eik_step(eik *eik)
    void eik_solve(eik *eik)
    void eik_add_trial(eik *eik, int2 ind, jet2 jet)
    void eik_add_valid(eik *eik, int2 ind, jet2 jet)
    void eik_get_shape(const eik *eik, int2 shape)
    jet2 *eik_get_jets_ptr(const eik *eik)
    state *eik_get_states_ptr(const eik *eik)
    dbl eik_T(const eik *eik, dbl2 xy)
    dbl eik_Tx(const eik *eik, dbl2 xy)
    dbl eik_Ty(const eik *eik, dbl2 xy)
    dbl eik_Txx(const eik *eik, dbl2 xy)
    dbl eik_Txy(const eik *eik, dbl2 xy)
    dbl eik_Tyy(const eik *eik, dbl2 xy)
    void eik_build_cells(eik *eik)
    bicubic eik_get_bicubic(const eik *eik, int2 indc)
    par2 eik_get_par(const eik *eik, int2 ind)
    bool eik_has_par(const eik *eik, int2 ind)
    const size_t *eik_get_accepted_ptr(const eik *eik)

cdef class Eik:
    cdef:
        eik *eik
        ArrayView state_view
        ArrayView T_view
        ArrayView Tx_view
        ArrayView Ty_view
        ArrayView Txy_view
        ArrayView accepted_view

cdef extern from "eik3.h":
    cdef struct eik3:
        pass

    void eik3_alloc(eik3 **eik)
    void eik3_dealloc(eik3 **eik)
    void eik3_init(eik3 *eik, const mesh3 *mesh, ftype ftype, const eik3 *orig)
    void eik3_deinit(eik3 *eik)
    size_t eik3_peek(const eik3 *eik)
    size_t eik3_step(eik3 *eik)
    void eik3_solve(eik3 *eik)
    bool eik3_is_solved(const eik3 *eik)
    const mesh3 *eik3_get_mesh(const eik3 *eik)
    bool eik3_is_far(const eik3 *eik, size_t ind)
    bool eik3_is_trial(const eik3 *eik, size_t ind)
    bool eik3_is_valid(const eik3 *eik, size_t ind)
    jet3 *eik3_get_jet_ptr(const eik3 *eik)
    dbl33 *eik3_get_hess_ptr(const eik3 *eik)
    state *eik3_get_state_ptr(const eik3 *eik)
    par3 eik3_get_par(const eik3 *eik, size_t l)
    bool eik3_has_par(const eik3 *eik, size_t l)
    dbl *eik3_get_t_in_ptr(const eik3 *eik)
    dbl *eik3_get_t_out_ptr(const eik3 *eik)
    void eik3_add_trial_w_data(eik3 *eik, size_t l, jet3 jet, dbl hess[3][3],
                               const dbl t_in[3], const dbl t_out[3])
    void eik3_add_pt_src_BCs(eik3 *eik, size_t l, jet3 jet)
    void eik3_add_refl_BCs(eik3 *eik, const size_t lf[3], const jet3 jet[3],
                           const dbl hess[3][3][3], const dbl t_in[3][3])
    void eik3_add_diff_edge_BCs(eik3 *eik, const size_t le[2],
                                const bb31 *T, const dbl rho1[2],
                                const dbl3 t_in[2])
    ftype eik3_get_ftype(const eik3 *eik)
    dbl eik3_get_slerp_tol(const eik3 *eik)
    bool eik3_has_BCs(const eik3 *eik, size_t l)
    void eik3_transport_dbl(const eik3 *eik, dbl *values, bool skip_filled)
    void eik3_transport_dblz(const eik3 *eik, dblz *values, bool skip_filled)
    void eik3_transport_curvature(const eik3 *eik, dbl *kappa, bool skip_filled)
    dbl eik3_get_h(const eik3* eik)
    const size_t *eik3_get_accepted_ptr(const eik3* eik)

cdef class Eik3:
    cdef:
        eik3 *eik
        ArrayView jet_view
        ArrayView T_view
        ArrayView grad_T_view
        ArrayView hess_view
        ArrayView state_view
        ArrayView t_in_view
        ArrayView t_out_view
        ArrayView accepted_view
