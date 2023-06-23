# cython: language_level=3

from cpython cimport bool

from libc.stdio cimport printf
from libc.stdlib cimport malloc

import numpy as np

from enum import Enum

NAN = np.nan

cdef extern from "jmm/common.h":
    struct eik3:
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

    def __dealloc__(self):
        mesh3_dealloc(&self.mesh)

    def __init__(self, Mesh3Data mesh_data, bool compute_bd_info=True, eps=None):
        cdef dbl eps_ = np.nan if eps is None else eps
        mesh3_init(self.mesh, &mesh_data.data, compute_bd_info, &eps_)

    @property
    def cells(self):
        cdef size_t ncells = mesh3_ncells(self.mesh)
        return np.asarray(<const size_t[:ncells, :4]> mesh3_get_cells_ptr(self.mesh))

    @property
    def verts(self):
        cdef size_t nverts = mesh3_nverts(self.mesh)
        return np.asarray(<const dbl[:nverts, :3]> mesh3_get_verts_ptr(self.mesh))

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

cdef class Eik3:
    cdef eik3 *eik

    def __cinit__(self):
        eik3_alloc(&self.eik)

    def __dealloc__(self):
        eik3_dealloc(&self.eik)

    def __init__(self, Mesh3 mesh, Sfunc sfunc):
        eik3_init(self.eik, mesh.mesh, &sfunc.sfunc)
