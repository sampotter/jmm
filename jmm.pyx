# cython: embedsignature=True
# cython: language_level=3

import numpy as np

import array

from cython cimport Py_buffer

from libc.stdlib cimport free, malloc

from enum import Enum

ctypedef bint bool

cdef extern from "def.h":
    ctypedef double dbl
    cdef enum state:
        FAR
        TRIAL
        VALID
        BOUNDARY
        NEW_VALID
    cdef enum stype:
        CONSTANT
        NUM_STYPE
    cdef enum error:
        SUCCESS
        BAD_ARGUMENT

cdef extern from "immintrin.h":
    ctypedef double __m256d

cdef extern from "vec.h":
    ctypedef struct dvec3:
        dbl data[4]
        __m256d packed

cdef extern from "bb.h":
    void bb3tri_interp3(dbl *f, dvec3 *Df, dvec3 *x, dbl *c)
    dbl bb3tri(dbl *c, dbl *b)
    dbl dbb3tri(dbl *c, dbl *b, dbl *a)
    dbl d2bb3tri(dbl *c, dbl *b, dbl *a1, dbl *a2)

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

cdef class ArrayView:
    cdef:
        bool readonly
        int ndim
        void *ptr
        Py_ssize_t *shape
        Py_ssize_t *strides
        char *format
        size_t itemsize

    def __cinit__(self, int ndim):
        self.ndim = ndim
        self.shape = <Py_ssize_t *>malloc(sizeof(Py_ssize_t)*ndim)
        self.strides = <Py_ssize_t *> malloc(sizeof(Py_ssize_t)*ndim)

    def __dealloc__(self):
        free(self.shape)
        free(self.strides)

    def __getbuffer__(self, Py_buffer *buf, int flags):
        buf.buf = <char *>self.ptr
        buf.format = self.format
        buf.internal = NULL
        buf.itemsize = self.itemsize
        buf.len = self.size
        buf.ndim = self.ndim
        buf.obj = self
        buf.readonly = self.readonly
        buf.shape = self.shape
        buf.strides = self.strides
        buf.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buf):
        pass

    @property
    def size(self):
        cdef Py_ssize_t size = 1
        for i in range(self.ndim):
            size *= self.shape[i]
        return size


cdef class Bb3Tri:
    cdef:
        dbl[::1] c

    def __cinit__(self, dbl[::1] f, dbl[:, ::1] Df, dbl[:, ::1] x):
        self.c = np.empty((10,), dtype=np.float64)

        cdef int i, j

        cdef dvec3 Df_[3]
        for i in range(3):
            for j in range(3):
                Df_[i].data[j] = Df[i, j]

        cdef dvec3 x_[3]
        for i in range(3):
            for j in range(3):
                x_[i].data[j] = x[i, j]

        bb3tri_interp3(&f[0], Df_, x_, &self.c[0])

    def f(self, dbl[::1] b):
        return bb3tri(&self.c[0], &b[0])

    def Df(self, dbl[::1] b, dbl[::1] a):
        return dbb3tri(&self.c[0], &b[0], &a[0])

    def D2f(self, dbl[::1] b, dbl[::1] a1, dbl[::1] a2):
        return d2bb3tri(&self.c[0], &b[0], &a1[0], &a2[0])

cdef class _Dial3:
    cdef:
        dial3 *dial
        Py_ssize_t shape[3]
        ArrayView Toff_view
        ArrayView xsrc_view
        ArrayView state_view

    def __cinit__(self, stype stype, int[:] shape, dbl h):
        dial3_alloc(&self.dial)
        dial3_init(self.dial, stype, &shape[0], h)

        self.shape[0] = shape[0]
        self.shape[1] = shape[1]
        self.shape[2] = shape[2]

        # Strides that haven't been scaled by the size of the
        # underlying type
        cdef Py_ssize_t base_strides[3]
        base_strides[2] = 1
        base_strides[1] = self.shape[2]
        base_strides[0] = self.shape[2]*self.shape[1]

        self.Toff_view = ArrayView(3)
        self.Toff_view.readonly = False
        self.Toff_view.ptr = <void *>dial3_get_Toff_ptr(self.dial)
        self.Toff_view.shape[0] = self.shape[0]
        self.Toff_view.shape[1] = self.shape[1]
        self.Toff_view.shape[2] = self.shape[2]
        self.Toff_view.strides[0] = sizeof(dbl)*base_strides[0]
        self.Toff_view.strides[1] = sizeof(dbl)*base_strides[1]
        self.Toff_view.strides[2] = sizeof(dbl)*base_strides[2]
        self.Toff_view.format = 'd'
        self.Toff_view.itemsize = sizeof(dbl)

        self.xsrc_view = ArrayView(4)
        self.xsrc_view.readonly = False
        self.xsrc_view.ptr = <void *>dial3_get_xsrc_ptr(self.dial)
        self.xsrc_view.shape[0] = self.shape[0]
        self.xsrc_view.shape[1] = self.shape[1]
        self.xsrc_view.shape[2] = self.shape[2]
        self.xsrc_view.shape[3] = 3
        self.xsrc_view.strides[0] = 4*sizeof(dbl)*base_strides[0]
        self.xsrc_view.strides[1] = 4*sizeof(dbl)*base_strides[1]
        self.xsrc_view.strides[2] = 4*sizeof(dbl)*base_strides[2]
        self.xsrc_view.strides[3] = sizeof(dbl)
        self.xsrc_view.format = 'd'
        self.xsrc_view.itemsize = sizeof(dbl)

        self.state_view = ArrayView(3)
        self.state_view.readonly = False
        self.state_view.ptr = <void *>dial3_get_state_ptr(self.dial)
        self.state_view.shape[0] = self.shape[0]
        self.state_view.shape[1] = self.shape[1]
        self.state_view.shape[2] = self.shape[2]
        self.state_view.strides[0] = sizeof(state)*base_strides[0]
        self.state_view.strides[1] = sizeof(state)*base_strides[1]
        self.state_view.strides[2] = sizeof(state)*base_strides[2]
        self.state_view.format = 'i'
        self.state_view.itemsize = sizeof(state)

    def __dealloc__(self):
        dial3_deinit(self.dial)
        dial3_dealloc(&self.dial)

    def add_point_source(self, int[:] ind0, dbl Toff):
        dial3_add_point_source(self.dial, &ind0[0], Toff)

    def add_boundary_points(self, int[::1, :] inds):
        # TODO: handle the case where inds is in a weird format
        if inds.shape[0] != 3:
            raise Exception('inds must be an 3xN array')
        dial3_add_boundary_points(self.dial, &inds[0, 0], inds.shape[1])

    def step(self):
        dial3_step(self.dial)

    def solve(self):
        dial3_solve(self.dial)

    @property
    def Toff(self):
        return self.Toff_view

    @property
    def xsrc(self):
        return self.xsrc_view

    @property
    def state(self):
        return self.state_view

class Stype(Enum):
    Constant = 0

class State(Enum):
    Far = 0
    Trial = 1
    Valid = 2
    Boundary = 3
    AdjacentToBoundary = 4
    NewValid = 5

class Dial:

    def __init__(self, stype, shape, h):
        self.shape = shape
        self.h = h
        if len(self.shape) == 3:
            self._dial = _Dial3(stype.value, array.array('i', [*shape]), h)
        else:
            raise Exception('len(shape) == %d not supported yet' % len(shape))

    def add_point_source(self, ind0, Toff):
        self._dial.add_point_source(array.array('i', [*ind0]), Toff)

    def add_boundary_points(self, inds):
        self._dial.add_boundary_points(inds)

    def step(self):
        self._dial.step()

    def solve(self):
        self._dial.solve()

    @property
    def _x(self):
        x = np.linspace(0, self.h*self.shape[0], self.shape[0])
        return x.reshape(self.shape[0], 1, 1)

    @property
    def _y(self):
        y = np.linspace(0, self.h*self.shape[1], self.shape[1])
        return y.reshape(1, self.shape[1], 1)

    @property
    def _z(self):
        z = np.linspace(0, self.h*self.shape[2], self.shape[2])
        return z.reshape(1, 1, self.shape[2])

    @property
    def T(self):
        dx = self._x - self.xsrc[:, :, :, 0]
        dy = self._y - self.xsrc[:, :, :, 1]
        dz = self._z - self.xsrc[:, :, :, 2]
        return self.Toff + np.sqrt(dx**2 + dy**2 + dz**2)

    @property
    def Toff(self):
        return np.asarray(self._dial.Toff)

    @property
    def xsrc(self):
        return np.asarray(self._dial.xsrc)

    @property
    def state(self):
        return np.asarray(self._dial.state)


cdef extern from "index.h":
    cdef struct ind2:
        size_t data[2]


cdef extern from "mesh3.h":
    struct mesh3:
        pass
    void mesh3_alloc(mesh3 **mesh)
    void mesh3_dealloc(mesh3 **mesh)
    void mesh3_init(mesh3 *mesh,
                    dbl *verts, size_t nverts,
                    size_t *cells, size_t ncells)
    void mesh3_deinit(mesh3 *mesh)
    int mesh3_nvv(mesh3 *mesh, size_t i)
    void mesh3_vv(mesh3 *mesh, size_t i, size_t *vv)
    int mesh3_nve(mesh3 *mesh, size_t i)
    void mesh3_ve(mesh3 *mesh, size_t i, size_t *ve)
    int mesh3_nvf(mesh3 *mesh, size_t i)
    void mesh3_vf(mesh3 *mesh, size_t i, size_t *vf)
    int mesh3_nvc(mesh3 *mesh, size_t i)
    void mesh3_vc(mesh3 *mesh, size_t i, size_t *vc)
    int mesh3_ncc(mesh3 *mesh, size_t i)
    void mesh3_cc(mesh3 *mesh, size_t i, size_t *cc)
    bool mesh3_bdc(mesh3 *mesh, size_t i)


cdef class Mesh3:
    cdef:
        mesh3 *mesh

    def __cinit__(self, dbl[:, ::1] verts, size_t[:, ::1] cells):
        mesh3_alloc(&self.mesh)
        cdef size_t nverts = verts.shape[0]
        cdef size_t ncells = cells.shape[0]
        mesh3_init(self.mesh, &verts[0, 0], nverts, &cells[0, 0], ncells)

    def __dealloc__(self):
        mesh3_deinit(self.mesh)
        mesh3_dealloc(&self.mesh)

    def vv(self, size_t i):
        cdef int nvv = mesh3_nvv(self.mesh, i)
        cdef size_t[::1] vv = np.empty((nvv,), dtype=np.uintp)
        mesh3_vv(self.mesh, i, &vv[0])
        return np.asarray(vv)

    def ve(self, size_t i):
        cdef int nve = mesh3_nve(self.mesh, i)
        cdef size_t[:, ::1] ve = np.empty((nve, 2), dtype=np.uintp)
        mesh3_ve(self.mesh, i, &ve[0, 0])
        return np.asarray(ve)

    def vf(self, size_t i):
        cdef int nvf = mesh3_nvf(self.mesh, i)
        cdef size_t[:, ::1] vf = np.empty((nvf, 3), dtype=np.uintp)
        mesh3_vf(self.mesh, i, &vf[0, 0])
        return np.asarray(vf)

    def vc(self, size_t i):
        cdef int nvc = mesh3_nvc(self.mesh, i)
        cdef size_t[::1] vc = np.empty((nvc,), dtype=np.uintp)
        mesh3_vc(self.mesh, i, &vc[0])
        return np.asarray(vc)

    def cc(self, size_t i):
        cdef int ncc = mesh3_ncc(self.mesh, i)
        cdef size_t[::1] cc = np.empty((ncc,), dtype=np.uintp)
        mesh3_cc(self.mesh, i, &cc[0])
        return np.asarray(cc)

    def bdc(self, size_t i):
        return mesh3_bdc(self.mesh, i)
