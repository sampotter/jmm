import numpy as np

from jmm.defs import Ftype

cimport jmm.defs

from jmm.bb cimport Bb31, Bb33
from jmm.defs cimport int2
from jmm.field cimport SlownessFunc2
from jmm.grid cimport Grid2, Grid3
from jmm.jet cimport jet2, Jet2, Jet3
from jmm.mesh cimport mesh3, Mesh3, mesh3_nverts
from jmm.par cimport par3, Parent3
from jmm.xfer cimport xfer

cdef class Eik:
    def __init__(self, *args):
        if len(args) > 0:
            raise RuntimeError('construct Eik using factory functions')

    def __cinit__(self, *args):
        if len(args) > 0:
            raise RuntimeError('construct Eik using factory functions')

    @staticmethod
    def from_s_and_grid(SlownessFunc2 s, Grid2 grid):
        eik = Eik()
        eik_alloc(&eik.eik)
        eik_init(eik.eik, &s.field, &grid._grid)
        eik._init_views()
        return eik

    def _init_views(self):
        self.state_view = ArrayView(2)
        self.state_view.readonly = True
        self.state_view.ptr = <void *>eik_get_states_ptr(self.eik)
        self.state_view.shape[0] = self.shape[0]
        self.state_view.shape[1] = self.shape[1]
        self.state_view.strides[0] = self.shape[1]*sizeof(state)
        self.state_view.strides[1] = 1*sizeof(state)
        self.state_view.format = 'i'
        self.state_view.itemsize = sizeof(state)

        self.T_view = ArrayView(2)
        self.T_view.readonly = True
        self.T_view.ptr = <void *>&eik_get_jets_ptr(self.eik).f
        self.T_view.shape[0] = self.shape[0]
        self.T_view.shape[1] = self.shape[1]
        self.T_view.strides[0] = self.shape[1]*sizeof(jet2)
        self.T_view.strides[1] = sizeof(jet2)
        self.T_view.format = 'd'
        self.T_view.itemsize = sizeof(dbl)

        self.Tx_view = ArrayView(2)
        self.Tx_view.readonly = True
        self.Tx_view.ptr = <void *>&eik_get_jets_ptr(self.eik).fx
        self.Tx_view.shape[0] = self.shape[0]
        self.Tx_view.shape[1] = self.shape[1]
        self.Tx_view.strides[0] = self.shape[1]*sizeof(jet2)
        self.Tx_view.strides[1] = sizeof(jet2)
        self.Tx_view.format = 'd'
        self.Tx_view.itemsize = sizeof(dbl)

        self.Ty_view = ArrayView(2)
        self.Ty_view.readonly = True
        self.Ty_view.ptr = <void *>&eik_get_jets_ptr(self.eik).fy
        self.Ty_view.shape[0] = self.shape[0]
        self.Ty_view.shape[1] = self.shape[1]
        self.Ty_view.strides[0] = self.shape[1]*sizeof(jet2)
        self.Ty_view.strides[1] = sizeof(jet2)
        self.Ty_view.format = 'd'
        self.Ty_view.itemsize = sizeof(dbl)

        self.Txy_view = ArrayView(2)
        self.Txy_view.readonly = True
        self.Txy_view.ptr = <void *>&eik_get_jets_ptr(self.eik).fxy
        self.Txy_view.shape[0] = self.shape[0]
        self.Txy_view.shape[1] = self.shape[1]
        self.Txy_view.strides[0] = self.shape[1]*sizeof(jet2)
        self.Txy_view.strides[1] = sizeof(jet2)
        self.Txy_view.format = 'd'
        self.Txy_view.itemsize = sizeof(dbl)

    @property
    def shape(self):
        cdef int[::1] shape = np.empty((2,), dtype=np.intc)
        eik_get_shape(self.eik, &shape[0])
        return np.asarray(shape)

    def add_valid(self, int[::1] ind, Jet2 jet):
        eik_add_valid(self.eik, &ind[0], jet.jet)

    def add_trial(self, int[::1] ind, Jet2 jet):
        eik_add_trial(self.eik, &ind[0], jet.jet)

    def build_cells(self):
        eik_build_cells(self.eik)

    def step(self):
        eik_step(self.eik)

    def solve(self):
        eik_solve(self.eik)

    @property
    def T(self):
        return np.asarray(self.T_view)

    @property
    def Tx(self):
        return np.asarray(self.Tx_view)

    @property
    def Ty(self):
        return np.asarray(self.Ty_view)

    @property
    def Txy(self):
        return np.asarray(self.Txy_view)

cdef class Eik3:
    def __init__(self, *args):
        if len(args) > 0:
            raise RuntimeError('construct Eik3 using factory functions')

    def __cinit__(self, *args):
        if len(args) > 0:
            raise RuntimeError('construct Eik3 using factory functions')

    @staticmethod
    def from_mesh_and_ftype(Mesh3 mesh, ftype, Eik3 orig=None):
        if not isinstance(ftype, Ftype):
            raise ValueError('ftype argument not an jmm.defs.Ftype instance')

        cdef const eik3 *orig_ptr = NULL
        if orig is not None:
            orig_ptr = orig.eik

        eik = Eik3()
        eik3_alloc(&eik.eik)
        eik3_init(eik.eik, mesh.mesh, ftype.value, orig_ptr)
        eik._init_views()
        return eik

    def _init_views(self):
        self.jet_view = ArrayView(1)
        self.jet_view.readonly = True
        self.jet_view.ptr = <void *>eik3_get_jet_ptr(self.eik)
        self.jet_view.shape[0] = self.size
        self.jet_view.strides[0] = 4*sizeof(dbl)
        self.jet_view.format = 'dddd'
        self.jet_view.itemsize = 4*sizeof(dbl)

        self.T_view = ArrayView(1)
        self.T_view.readonly = True
        self.T_view.ptr = <void *>&eik3_get_jet_ptr(self.eik).f
        self.T_view.shape[0] = self.size
        self.T_view.strides[0] = 4*sizeof(dbl)
        self.T_view.format = 'd'
        self.T_view.itemsize = 1*sizeof(dbl)

        self.grad_T_view = ArrayView(2)
        self.grad_T_view.readonly = True
        self.grad_T_view.ptr = <void *>&eik3_get_jet_ptr(self.eik).fx
        self.grad_T_view.shape[0] = self.size
        self.grad_T_view.shape[1] = 3
        self.grad_T_view.strides[0] = 4*sizeof(dbl)
        self.grad_T_view.strides[1] = 1*sizeof(dbl)
        self.grad_T_view.format = 'd'
        self.grad_T_view.itemsize = 1*sizeof(dbl)

        self.hess_view = ArrayView(3)
        self.hess_view.readonly = True
        self.hess_view.ptr = <void *>eik3_get_hess_ptr(self.eik)
        self.hess_view.shape[0] = self.size
        self.hess_view.shape[1] = 3
        self.hess_view.shape[2] = 3
        self.hess_view.strides[0] = 9*sizeof(dbl)
        self.hess_view.strides[1] = 3*sizeof(dbl)
        self.hess_view.strides[2] = 1*sizeof(dbl)
        self.hess_view.format = 'd'
        self.hess_view.itemsize = 1*sizeof(dbl)

        self.state_view = ArrayView(1)
        self.state_view.readonly = True
        self.state_view.ptr = <void *>eik3_get_state_ptr(self.eik)
        self.state_view.shape[0] = self.size
        self.state_view.strides[0] = sizeof(state)
        self.state_view.format = 'i'
        self.state_view.itemsize = sizeof(state)

        if self.ftype != Ftype.PointSource:
            self.t_in_view = ArrayView(2)
            self.t_in_view.readonly = True
            self.t_in_view.ptr = <void *>eik3_get_t_in_ptr(self.eik)
            self.t_in_view.shape[0] = self.size
            self.t_in_view.shape[1] = 3
            self.t_in_view.strides[0] = 3*sizeof(dbl)
            self.t_in_view.strides[1] = 1*sizeof(dbl)
            self.t_in_view.format = 'd'
            self.t_in_view.itemsize = sizeof(dbl)

        self.t_out_view = ArrayView(2)
        self.t_out_view.readonly = True
        self.t_out_view.ptr = <void *>eik3_get_t_out_ptr(self.eik)
        self.t_out_view.shape[0] = self.size
        self.t_out_view.shape[1] = 3
        self.t_out_view.strides[0] = 3*sizeof(dbl)
        self.t_out_view.strides[1] = 1*sizeof(dbl)
        self.t_out_view.format = 'd'
        self.t_out_view.itemsize = sizeof(dbl)

        # TODO: the buffer format character 'Q' below should be 'N'
        # according to:
        #
        #   https://docs.python.org/3/library/struct.html#module-struct
        #
        # but this fails at run-time with:
        #
        #   "ValueError: 'N' is not a valid PEP 3118 buffer format string"
        self.accepted_view = ArrayView(1)
        self.accepted_view.readonly = True
        self.accepted_view.ptr = <void *>eik3_get_accepted_ptr(self.eik)
        self.accepted_view.shape[0] = self.size
        self.accepted_view.strides[0] = sizeof(size_t)
        self.accepted_view.format = 'Q'
        self.accepted_view.itemsize = sizeof(size_t)

    def __dealloc__(self):
        eik3_deinit(self.eik)
        eik3_dealloc(&self.eik)

    def peek(self):
        return eik3_peek(self.eik)

    def step(self):
        return eik3_step(self.eik)

    def solve(self):
        eik3_solve(self.eik)

    @property
    def is_solved(self):
        return eik3_is_solved(self.eik)

    def is_far(self, size_t ind):
        return eik3_is_far(self.eik, ind)

    def is_trial(self, size_t ind):
        return eik3_is_trial(self.eik, ind)

    def is_valid(self, size_t ind):
        return eik3_is_valid(self.eik, ind)

    def transfer_solution_to_grid(self, Grid3 grid):
        cdef dbl[::1] y = np.empty((grid.size,), dtype=np.float64)
        xfer(eik3_get_mesh(self.eik), eik3_get_jet_ptr(self.eik),
             &grid._grid, &y[0])
        return np.asarray(y).reshape(grid.shape)

    @property
    def front(self):
        return self.peek()

    @property
    def size(self):
        cdef const mesh3 *mesh = eik3_get_mesh(self.eik)
        return mesh3_nverts(mesh)

    @property
    def jet(self):
        return np.asarray(self.jet_view)

    @property
    def T(self):
        return np.asarray(self.T_view)

    @property
    def grad_T(self):
        return np.asarray(self.grad_T_view)

    @property
    def hess(self):
        return np.asarray(self.hess_view)

    @property
    def state(self):
        return np.asarray(self.state_view)

    def has_par(self, ind):
        return eik3_has_par(self.eik, ind)

    def get_par(self, ind):
        if not eik3_has_par(self.eik, ind):
            raise ValueError('node %d has no parent' % ind)
        return Parent3(eik3_get_par(self.eik, ind))

    @property
    def mesh(self):
        cdef mesh3 *mesh = eik3_get_mesh(self.eik)
        return Mesh3.from_ptr(mesh)

    def get_bezier_tetra(self, size_t cell_ind):
        vert_inds = self.mesh.cells[cell_ind]
        jets = self.jet[vert_inds]
        cdef dbl[::1] f = np.array([_[0] for _ in jets])
        cdef dbl[:, ::1] Df = np.array([(_[1], _[2], _[3]) for _ in jets])
        cdef dbl[:, ::1] x = self.mesh.verts[vert_inds]
        return Bb33.from_3d_data(f, Df, x)

    @property
    def t_in(self):
        if self.ftype == Ftype.PointSource:
            raise RuntimeError('t_in vector only defined for scattered fields')
        return np.asarray(self.t_in_view)

    @property
    def t_out(self):
        return np.asarray(self.t_out_view)

    @property
    def ftype(self):
        ftype = Ftype(eik3_get_ftype(self.eik))
        return ftype

    def add_trial(self, size_t l, Jet3 jet, dbl[:, ::1] hess,
                  dbl[::1] t_in, dbl[::1] t_out):
        eik3_add_trial_w_data(self.eik, l, jet.jet,
                              <dbl(*)[3]>&hess[0, 0],
                              &t_in[0],
                              &t_out[0])

    def add_pt_src_BCs(self, size_t l, Jet3 jet):
        if self.ftype != Ftype.PointSource:
            raise RuntimeError(
                'tried to add point source BCs to %s field' % str(self.ftype))
        eik3_add_pt_src_BCs(self.eik, l, jet.jet)

    def add_refl_BCs(self, size_t l0, size_t l1, size_t l2,
                     Jet3 jet0, Jet3 jet1, Jet3 jet2,
                     dbl[:, :, ::1] hess, dbl[:, ::1] t_in):
        if self.ftype != Ftype.Reflection:
            raise RuntimeError(
                'tried to add reflection BCs to %s field' % str(self.ftype))
        if t_in.shape[0] != 3 or t_in.shape[1] != 3:
            raise ValueError('t_in should have shape (3, 3)')
        cdef size_t[3] lf = [l0, l1, l2]
        cdef jet3[3] jet = [jet0.jet, jet1.jet, jet2.jet]
        eik3_add_refl_BCs(self.eik, lf, jet,
                          <const dbl33 *>&hess[0, 0, 0],
                          <const dbl(*)[3]>&t_in[0, 0])

    def add_diff_edge_BCs(self, size_t[::1] le, Bb31 T,
                          dbl[::1] rho1, dbl[:, ::1] t_in):
        if self.ftype != Ftype.EdgeDiffraction:
            raise RuntimeError(
                'tried to add diff. edge BCs to a %s' % str(self.ftype))
        if not T.isfinite():
            raise ValueError(
                'tried to add diff. edge BCs with bad edge cubic for T')
        eik3_add_diff_edge_BCs(
            self.eik, &le[0], &T.bb, &rho1[0], <const dbl(*)[3]>&t_in[0, 0])

    @property
    def slerp_tol(self):
        return eik3_get_slerp_tol(self.eik)

    def has_BCs(self, size_t l):
        return eik3_has_BCs(self.eik, l)

    def transport_scalars(self, dbl_or_dblz[::1] values, bool skip_filled):
        if values.size != mesh3_nverts(eik3_get_mesh(self.eik)):
            raise ValueError('must pass as many values as domain points')
        if dbl_or_dblz is dbl:
            eik3_transport_dbl(self.eik, &values[0], skip_filled)
        if dbl_or_dblz is dblz:
            eik3_transport_dblz(self.eik, &values[0], skip_filled)
        return np.asarray(values)

    def transport_curvature(self, dbl[::1] values, bool skip_filled):
        if values.size != mesh3_nverts(eik3_get_mesh(self.eik)):
            raise ValueError('must pass as many values as domain points')
        eik3_transport_curvature(self.eik, &values[0], skip_filled)
        return np.asarray(values)

    @property
    def h(self):
        return eik3_get_h(self.eik)

    @property
    def accepted(self):
        return np.asarray(self.accepted_view)
