from jmm.array_view cimport ArrayView
from jmm.defs cimport bool, dbl, error, state, stype

cdef extern from "dial.h":
    cdef struct dial3:
        pass

    void dial3_alloc(dial3 **dial)
    error dial3_init(dial3 *dial, stype stype, const int *shape, dbl h)
    void dial3_deinit(dial3 *dial)
    void dial3_dealloc(dial3 **dial)
    void dial3_add_point_source(dial3 *dial, const int *ind0, dbl T)
    void dial3_add_boundary_points(dial3 *dial, const int *inds, size_t n)
    bool dial3_step(dial3 *dial)
    void dial3_solve(dial3 *dial)
    dbl dial3_get_T(const dial3 *dial, int l)
    void dial3_get_grad_T(const dial3 *dial, int l, dbl *grad_T)
    dbl *dial3_get_Toff_ptr(const dial3 *dial)
    dbl *dial3_get_xsrc_ptr(const dial3 *dial)
    state *dial3_get_state_ptr(const dial3 *dial)

cdef class _Dial3:
    cdef:
        dial3 *dial
        Py_ssize_t shape[3]
        ArrayView Toff_view
        ArrayView xsrc_view
        ArrayView state_view
