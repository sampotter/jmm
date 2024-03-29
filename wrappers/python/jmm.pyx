# cython: language_level=3

from cpython cimport bool

from libc.stdio cimport printf
from libc.stdlib cimport malloc

import numpy as np

from enum import Enum

NAN = np.nan

cdef size_t NO_INDEX = -1

cdef extern from "jmm/common.h":
    struct eik3:
        pass
    struct eik3hh_branch:
        pass
    struct eik3hh:
        pass
    struct mesh3:
        pass

cdef extern from "jmm/def.h":
    ctypedef double dbl
    ctypedef double[3] dbl3
    ctypedef double[3][3] dbl33

    cdef enum error:
        SUCCESS
        BAD_ARGUMENT

cdef extern from "jmm/jet.h":
    struct jet31t:
        dbl f
        dbl3 Df

cdef extern from "jmm/array.h":
    # NOTE: unlike the other structs in this file, we use the typdef
    # name here ("array_s") instead of the struct name ("array"),
    # since it collides with another struct named "array" in the C
    # code generated by Cython...
    #
    # TODO: fix this by namespacing the structs in the C library
    struct array_s:
        pass

    void array_dealloc(array_s **arr)
    void array_deinit(array_s *arr)
    size_t array_size(const array_s *arr)
    void array_get(const array_s *arr, size_t i, void *elt)

cdef extern from "jmm/bmesh.h":
    struct bmesh33:
        pass

    void bmesh33_alloc(bmesh33 **bmesh)
    void bmesh33_init_from_mesh3_and_jets(bmesh33 *bmesh, const mesh3 *mesh, const jet31t *jet)
    dbl bmesh33_f(const bmesh33 *bmesh, const dbl3 x)

cdef class Bmesh33:
    cdef bmesh33 *bmesh

    @staticmethod
    cdef from_ptr(bmesh33 *bmesh):
        cdef Bmesh33 _ = Bmesh33.__new__(Bmesh33)
        _.bmesh = bmesh
        return _

    def __call__(self, *args):
        _ = args[0] if len(args) == 1 else args
        cdef dbl3 x = [_[0], _[1], _[2]]
        return bmesh33_f(self.bmesh, x)

cdef extern from "jmm/mesh3.h":
    struct mesh3_data:
        pass

    void mesh3_data_init_from_bin(mesh3_data *data, const char *verts_path, const char *cells_path)
    void mesh3_data_init_from_off_file(mesh3_data *data, const char *path, dbl maxvol, bool verbose)
    void mesh3_data_deinit(mesh3_data *data)
    error mesh3_data_insert_vert(mesh3_data *data, const dbl3 x, dbl eps)

    void mesh3_alloc(mesh3 **mesh)
    void mesh3_dealloc(mesh3 **mesh)
    void mesh3_init(mesh3 *mesh, const mesh3_data *data, bool compute_bd_info, const dbl *eps)
    const size_t *mesh3_get_cells_ptr(const mesh3 *mesh)
    const dbl *mesh3_get_verts_ptr(const mesh3 *mesh)
    size_t mesh3_ncells(const mesh3 *mesh)
    size_t mesh3_nverts(const mesh3 *mesh)
    size_t mesh3_find_cell_containing_point(const mesh3 *mesh, const dbl x[3], size_t lc)

cdef class Mesh3Data:
    cdef mesh3_data data

    @staticmethod
    def from_bin(str verts_path, str cells_path):
        cdef bytes verts_path_bytes = verts_path.encode()
        cdef const char *verts_path_c_str = verts_path_bytes

        cdef bytes cells_path_bytes = cells_path.encode()
        cdef const char *cells_path_c_str = cells_path_bytes

        mesh_data = Mesh3Data()
        mesh3_data_init_from_bin(&mesh_data.data, verts_path_c_str, cells_path_c_str)

        return mesh_data

    @staticmethod
    def from_off(str path, dbl maxvol, bool verbose=False):
        '''Create a new Mesh3Data instance from an OFF file, running
        TetGen to generate the tetrahedron mesh.

        Args:
            path (str): the path to the OFF file
            maxvol (float): maximum volumne constraint for TetGen
            verbose (bool): whether to allow verbose output from TetGen

        Returns:
            The new Mesh3Data instance.

        '''
        cdef bytes path_bytes = path.encode()
        cdef const char *path_c_str = path_bytes

        mesh_data = Mesh3Data()
        mesh3_data_init_from_off_file(&mesh_data.data, path_c_str, maxvol, verbose)

        return mesh_data

    def insert_vert(self, const dbl[:] x, dbl eps):
        if len(x) != 3:
            raise ValueError('must have len(x) == 3')
        cdef dbl3 x_ = [x[0], x[1], x[2]]
        mesh3_data_insert_vert(&self.data, &x_[0], eps)

cdef class Mesh3:
    cdef mesh3 *mesh

    def __cinit__(self):
        mesh3_alloc(&self.mesh)

    def __init__(self, Mesh3Data mesh_data, bool compute_bd_info=True, eps=None):
        cdef dbl eps_ = np.nan if eps is None else eps
        mesh3_init(self.mesh, &mesh_data.data, compute_bd_info, &eps_)

    @staticmethod
    cdef from_ptr(mesh3 *mesh):
        cdef Mesh3 _ = Mesh3.__new__(Mesh3)
        _.mesh = mesh
        return _

    @property
    def cells(self):
        cdef size_t ncells = mesh3_ncells(self.mesh)
        return np.asarray(<const size_t[:ncells, :4]> mesh3_get_cells_ptr(self.mesh))

    @property
    def verts(self):
        cdef size_t nverts = mesh3_nverts(self.mesh)
        return np.asarray(<const dbl[:nverts, :3]> mesh3_get_verts_ptr(self.mesh))

    @property
    def ncells(self):
        return mesh3_ncells(self.mesh)

    @property
    def nverts(self):
        return mesh3_nverts(self.mesh)

    def find_cell_containing_point(self, *args):
        _ = args[0] if len(args) == 1 else args
        cdef dbl[3] x = [_[0], _[1], _[2]]
        return mesh3_find_cell_containing_point(self.mesh, x, NO_INDEX)

cdef extern from "jmm/slow.h":
    cdef enum stype:
        STYPE_CONSTANT = 0
        STYPE_FUNC_PTR = 1
        STYPE_JET31T = 2
        STYPE_NUM_STYPE = 3

    ctypedef dbl (*sfunc_s)(dbl3)
    ctypedef void (*sfunc_Ds)(dbl3, dbl3)
    ctypedef void (*sfunc_D2s)(dbl3, dbl33)

    cdef struct sfunc_funcs:
        sfunc_s s
        sfunc_Ds Ds
        sfunc_D2s D2s

    cdef struct sfunc:
        stype stype
        sfunc_funcs funcs
        jet31t *data_jet31t

class Stype(Enum):
    Constant = STYPE_CONSTANT
    Func = STYPE_FUNC_PTR
    Jet31t = STYPE_JET31T
    NumStype = STYPE_NUM_STYPE

cdef sfunc make_constant_sfunc():
    cdef sfunc sfunc
    sfunc.stype = STYPE_CONSTANT
    sfunc.funcs.s = NULL
    sfunc.funcs.Ds = NULL
    sfunc.funcs.D2s = NULL
    return sfunc

cdef class Sfunc:
    cdef sfunc sfunc

    Constant = Sfunc(Stype.Constant)

    def __init__(self, stype):
        if stype == Stype.Constant:
            self.sfunc = make_constant_sfunc()
        else:
            raise NotImplementedError(f'stype == {stype}')

cdef extern from "jmm/eik3.h":
    void eik3_alloc(eik3 **eik)
    void eik3_dealloc(eik3 **eik)
    void eik3_init(eik3 *eik, const mesh3 *mesh, const sfunc *sfunc)
    const mesh3 *eik3_get_mesh(const eik3 *eik)
    jet31t *eik3_get_jet_ptr(const eik3 *eik)

cdef class Eik3:
    cdef eik3 *eik

    def __cinit__(self):
        eik3_alloc(&self.eik)

    def __init__(self, Mesh3 mesh, Sfunc sfunc):
        eik3_init(self.eik, mesh.mesh, &sfunc.sfunc)

    @staticmethod
    cdef Eik3 from_ptr(eik3 *eik):
        cdef Eik3 _ = Eik3.__new__(Eik3)
        _.eik = eik
        return _

    def build_T_bmesh(self):
        cdef mesh3 *mesh = eik3_get_mesh(self.eik)
        cdef jet31t *jet = eik3_get_jet_ptr(self.eik)
        cdef bmesh33 *bmesh
        bmesh33_alloc(&bmesh)
        bmesh33_init_from_mesh3_and_jets(bmesh, mesh, jet)
        return Bmesh33.from_ptr(bmesh)

cdef extern from "jmm/eik3hh_branch.h":
    cdef enum eik3hh_branch_type:
        EIK3HH_BRANCH_TYPE_UNINITIALIZED
        EIK3HH_BRANCH_TYPE_PT_SRC
        EIK3HH_BRANCH_TYPE_REFL

    void eik3hh_branch_alloc(eik3hh_branch **branch)
    void eik3hh_branch_dealloc(eik3hh_branch **hh)
    void eik3hh_branch_init_pt_src(eik3hh_branch *branch, const eik3hh *hh, const dbl3 xsrc)
    void eik3hh_branch_solve(eik3hh_branch *branch, bool verbose)
    eik3 *eik3hh_branch_get_eik(const eik3hh_branch *branch)
    const dbl *eik3hh_branch_get_spread(const eik3hh_branch *branch)
    const dbl *eik3hh_branch_get_org(const eik3hh_branch *branch)
    array_s *eik3hh_branch_get_visible_refls(const eik3hh_branch *branch)
    eik3hh_branch *eik3hh_branch_add_refl(const eik3hh_branch *branch, size_t refl_index)

cdef class Eik3hhBranch:
    cdef eik3hh_branch *branch

    def __cinit__(self):
        eik3hh_branch_alloc(&self.branch)

    def solve(self, bool verbose=False):
        eik3hh_branch_solve(self.branch, verbose)

    def get_eik(self):
        cdef eik3 *eik = eik3hh_branch_get_eik(self.branch)
        return Eik3.from_ptr(eik)

    def get_visible_refls(self):
        cdef array_s *refl_inds_array = eik3hh_branch_get_visible_refls(self.branch)
        cdef size_t num_refls = array_size(refl_inds_array)
        refl_inds_lst = []
        cdef size_t i
        cdef size_t refl_ind
        for i in range(num_refls):
            array_get(refl_inds_array, i, &refl_ind)
            refl_inds_lst.append(refl_ind)
        return refl_inds_lst

    def add_refl(self, refl_index):
        refl = Eik3hhBranch()
        refl.branch = eik3hh_branch_add_refl(self.branch, refl_index)
        return refl

    @property
    def mesh(self):
        return Mesh3.from_ptr(eik3_get_mesh(eik3hh_branch_get_eik(self.branch)))

    @property
    def spread(self):
        return np.asarray(<const dbl[:self.mesh.nverts]> eik3hh_branch_get_spread(self.branch))

    @property
    def org(self):
        return np.asarray(<const dbl[:self.mesh.nverts]> eik3hh_branch_get_org(self.branch))

cdef extern from "jmm/eik3hh.h":
    void eik3hh_alloc(eik3hh **hh);
    void eik3hh_dealloc(eik3hh **hh);
    void eik3hh_init_with_pt_src(eik3hh *hh, const mesh3 *mesh, dbl c, dbl rfac, const dbl3 xsrc)
    eik3hh_branch *eik3hh_get_root_branch(eik3hh *hh)

cdef class Eik3hh:
    cdef eik3hh *hh

    def __cinit__(self):
        eik3hh_alloc(&self.hh)

    @staticmethod
    def new_with_pt_src(Mesh3 mesh, dbl c, dbl rfac, dbl[:] xsrc):
        cdef dbl3 xsrc_ = [xsrc[0], xsrc[1], xsrc[2]]

        hh = Eik3hh()
        eik3hh_init_with_pt_src(hh.hh, mesh.mesh, c, rfac, xsrc_)

        return hh

    def get_root_branch(self):
        branch = Eik3hhBranch()
        branch.branch = eik3hh_get_root_branch(self.hh)
        return branch
