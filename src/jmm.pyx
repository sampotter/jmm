import numpy as np

import array

from enum import Enum

from cython cimport Py_buffer

from libc.math cimport isfinite
from libc.stdlib cimport free, malloc, qsort
from libc.string cimport memcpy

from defs cimport bool, error, state, stype

from bb cimport *
from bmesh cimport *
from dial cimport *
from edge cimport *
from edgemap cimport *
from grid3 cimport *
from jet cimport *
from mesh2 cimport *
from mesh3 cimport *
from par cimport *
from rtree cimport *
from utetra cimport *
from xfer cimport *

class Stype(Enum):
    Constant = 0

class State(Enum):
    Far = 0
    Trial = 1
    Valid = 2
    Boundary = 3
    AdjacentToBoundary = 4
    NewValid = 5
    Shadow = 6

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

cdef class Bb33:
    cdef bb33 _bb

    @staticmethod
    def from_3d_data(dbl[:] f, dbl[:, :] Df, dbl[:, :] x):
        if f.size != 4 or f.shape[0] != 4:
            raise Exception('`f` must be a length 4 vector')
        if Df.size != 12 or Df.shape[0] != 4 or Df.shape[1] != 3:
            raise Exception('`Df` must have shape (4, 3)')
        if x.size != 12 or x.shape[0] != 4 or x.shape[1] != 3:
            raise Exception('`x` must have shape (4, 3)')
        bb = Bb33()
        bb33_init_from_3d_data(
            &bb._bb,
            &f[0],
            <const dbl (*)[3]>&Df[0, 0],
            <const dbl (*)[3]>&x[0, 0])
        return bb

    @staticmethod
    cdef from_bb33(bb33 bb):
        bb_ = Bb33()
        bb_._bb = bb
        return bb_

    def f(self, dbl[:] b):
        if b.size != 4 or b.shape[0] != 4:
            raise Exception('`b` must be a length 4 vector')
        return bb33_f(&self._bb, &b[0])

    def Df(self, dbl[:] b, dbl[:] a):
        if b.size != 4 or b.shape[0] != 4:
            raise Exception('`b` must be a length 4 vector')
        if a.size != 4 or a.shape[0] != 4:
            raise Exception('`a` must be a length 4 vector')
        return bb33_df(&self._bb, &b[0], &a[0])

    def D2f(self, dbl[:] b, dbl[:, :] a):
        if b.size != 4 or b.shape[0] != 4:
            raise Exception('`b` must be a length 4 vector')
        if a.size != 12 or a.shape[0] != 3 or a.shape[4] != 4:
            raise Exception('`a` must have shape (3, 4)')
        return bb33_d2f(&self._bb, &b[0], <const dbl (*)[4]>&a[0, 0])

    def convex_hull_brackets_value(self, dbl value):
        return bb33_convex_hull_brackets_value(&self._bb, value)

cdef class Bmesh33Cell:
    cdef bmesh33_cell cell

    @staticmethod
    cdef from_robj_ptr(const robj *obj):
        cell = Bmesh33Cell()
        cell.cell = (<bmesh33_cell *>robj_get_data(obj))[0]
        return cell

    @staticmethod
    cdef from_cell(bmesh33_cell cell):
        cell_ = Bmesh33Cell()
        cell_.cell = cell
        return cell_

    @property
    def index(self):
        return self.cell.l

cdef class Bmesh33:
    cdef:
        bool ptr_owner
        bmesh33 *bmesh

    def __dealloc__(self):
        if self.ptr_owner:
            bmesh33_deinit(self.bmesh)
            bmesh33_dealloc(&self.bmesh)

    @staticmethod
    def from_eik3(Eik3 eik):
        bmesh = Bmesh33()
        bmesh.ptr_owner = True
        bmesh33_alloc(&bmesh.bmesh)
        cdef const mesh3 *mesh = eik3_get_mesh(eik.eik)
        cdef const jet3 *jet = eik3_get_jet_ptr(eik.eik)
        bmesh33_init_from_mesh3_and_jets(bmesh.bmesh, mesh, jet)
        return bmesh

    @staticmethod
    def from_mesh_and_jets(Mesh3 mesh, dbl[::1] f, dbl[:, ::1] Df):
        cdef jet3 *jet = <jet3 *>malloc(f.size*sizeof(jet3))
        cdef int i
        for i in range(f.size):
            jet[i].f = f[i]
            jet[i].fx = Df[i, 0]
            jet[i].fy = Df[i, 1]
            jet[i].fz = Df[i, 2]
        bmesh = Bmesh33()
        bmesh.ptr_owner = True
        bmesh33_alloc(&bmesh.bmesh)
        bmesh33_init_from_mesh3_and_jets(bmesh.bmesh, mesh.mesh, jet)
        free(jet)
        return bmesh

    @staticmethod
    cdef from_ptr(bmesh33 *bmesh_ptr, ptr_owner=False):
        bmesh = Bmesh33()
        bmesh.ptr_owner = ptr_owner
        bmesh.bmesh = bmesh_ptr
        return bmesh

    @property
    def num_cells(self):
        return bmesh33_num_cells(self.bmesh)

    @property
    def mesh(self):
        return Mesh3.from_ptr(<mesh3 *>bmesh33_get_mesh_ptr(self.bmesh))

    def restrict_to_level(self, dbl level):
        return Bmesh33.from_ptr(
            bmesh33_restrict_to_level(self.bmesh, level), ptr_owner=True)

    def get_cell(self, size_t l):
        return Bmesh33Cell.from_cell(bmesh33_get_cell(self.bmesh, l))

cdef class Grid3:
    cdef grid3 _grid

    def __cinit__(self, int[:] dim, dbl[:] xmin, dbl h):
        self._grid.dim[0] = dim[0]
        self._grid.dim[1] = dim[1]
        self._grid.dim[2] = dim[2]
        self._grid.min[0] = xmin[0]
        self._grid.min[1] = xmin[1]
        self._grid.min[2] = xmin[2]
        self._grid.h = h

    @property
    def size(self):
        return self._grid.dim[0]*self._grid.dim[1]*self._grid.dim[2]

    @property
    def shape(self):
        return np.array([self._grid.dim[0],
                         self._grid.dim[1],
                         self._grid.dim[2]])

cdef class Rect3:
    cdef rect3 _rect

    def __cinit__(self, rect3 rect):
        self._rect = rect

    def __repr__(self):
        return f'Rect3(xmin={self._rect.min[0]}, xmax={self._rect.max[0]}, zmin={self._rect.min[1]}, zmax={self._rect.max[1]}, zmin={self._rect.min[2]}, zmax={self._rect.max[2]})'

    def __str__(self):
        return f'[{self._rect.min[0]}, {self._rect.max[0]}] x [{self._rect.min[1]}, {self._rect.max[1]}] x [{self._rect.min[2]}, {self._rect.max[2]}]'

cdef class Ray3:
    '''A ray in 3D, consisting of an `origin` and `direction`. The
`direction` vector doesn't need to be parametrized, but its magnitude
_does_ determine the parametrization of the ray. That is, for a
parameter `t`, a point `p(t)` on the ray equals `p(t) = origin +
direction*t`.

    :param origin: A 3D vector indicating the start of the ray,
        typically an instance of :class:`numpy.ndarray`.
    :param direction: The direction of ray. See comment above about
        normalization. Also typically an :class:`numpy.ndarray`
        instance.
    '''

    cdef ray3 _ray

    def __cinit__(self, dbl[::1] origin, dbl[::1] direction):
        memcpy(&self._ray.org[0], &origin[0], 3*sizeof(dbl))
        memcpy(&self._ray.dir[0], &direction[0], 3*sizeof(dbl))

    @staticmethod
    cdef from_ptr(const ray3 *ray):
        cdef dbl[::1] org = np.array([ray.org[0], ray.org[1], ray.org[2]])
        cdef dbl[::1] dir = np.array([ray.dir[0], ray.dir[1], ray.dir[2]])
        return Ray3(org, dir)

    def __repr__(self):
        return f'Ray3(origin={self._ray.org}, direction={self._ray.dir})'

cdef class Mesh2Tri:
    cdef mesh2_tri *_tri

    @staticmethod
    cdef from_robj_ptr(const robj *obj):
        tri = Mesh2Tri()
        tri._tri = <mesh2_tri *>robj_get_data(obj)
        return tri

    @property
    def mesh(self):
        return Mesh2.from_ptr(self._tri.mesh)

    @property
    def index(self):
        return self._tri.l

cdef class Mesh2:
    cdef:
        bool ptr_owner
        mesh2 *mesh
        ArrayView verts_view
        ArrayView faces_view

    def __dealloc__(self):
        if self.ptr_owner:
            mesh2_deinit(self.mesh)
            mesh2_dealloc(&self.mesh)

    @staticmethod
    def from_verts_and_faces(dbl[:, ::1] verts, size_t[:, ::1] faces):
        mesh = Mesh2()
        mesh.ptr_owner = True
        mesh2_alloc(&mesh.mesh)
        cdef size_t nverts = verts.shape[0]
        cdef size_t nfaces = faces.shape[0]
        mesh2_init(mesh.mesh, &verts[0, 0], nverts, &faces[0, 0], nfaces)
        mesh._set_views()
        return mesh

    @staticmethod
    cdef from_ptr(const mesh2 *mesh_ptr, ptr_owner=False):
        mesh = Mesh2()
        mesh.ptr_owner = ptr_owner
        mesh.mesh = <mesh2 *>mesh_ptr
        mesh._set_views()
        return mesh

    cdef _set_views(self):
        self.verts_view = ArrayView(2)
        self.verts_view.readonly = True
        self.verts_view.ptr = <void *>mesh2_get_points_ptr(self.mesh)
        self.verts_view.shape[0] = self.num_verts
        self.verts_view.shape[1] = 3
        self.verts_view.strides[0] = 3*sizeof(dbl)
        self.verts_view.strides[1] = 1*sizeof(dbl)
        self.verts_view.format = 'd'
        self.verts_view.itemsize = sizeof(dbl)

        self.faces_view = ArrayView(2)
        self.faces_view.readonly = True
        self.faces_view.ptr = <void *>mesh2_get_faces_ptr(self.mesh)
        self.faces_view.shape[0] = self.num_faces
        self.faces_view.shape[1] = 3
        self.faces_view.strides[0] = 3*sizeof(size_t)
        self.faces_view.strides[1] = 1*sizeof(size_t)
        self.faces_view.format = 'L'
        self.faces_view.itemsize = sizeof(size_t)

    @property
    def num_verts(self):
        return mesh2_get_num_points(self.mesh)

    @property
    def num_faces(self):
        return mesh2_get_num_faces(self.mesh)

    @property
    def verts(self):
        return np.asarray(self.verts_view)

    @property
    def faces(self):
        return np.asarray(self.faces_view)

    @property
    def bounding_box(self):
        return Rect3(mesh2_get_bounding_box(self.mesh))

class RobjType(Enum):
    Mesh2Tri = 0
    Mesh3Tetra = 1
    Tri3 = 2
    Tetra3 = 3

class RtreeSplitStrategy(Enum):
    SurfaceArea = 0

cdef class Robj:
    cdef const robj *_obj

    @staticmethod
    cdef from_ptr(const robj *obj):
        robj = Robj()
        robj._obj = obj
        return robj

    @property
    def type(self):
        return RobjType(robj_get_type(self._obj))

    @property
    def centroid(self):
        p = np.empty((3,), dtype=np.float64)
        cdef dbl[::1] p_mv = p
        robj_get_centroid(self._obj, &p_mv[0])
        return p

    def overlaps_box(self, Rect3 box):
        return robj_isects_bbox(self._obj, &box._rect)

    def intersect(self, Ray3 ray):
        cdef dbl t
        cdef bool hit = robj_intersect(self._obj, &ray._ray, &t)
        return hit, t

    def astype(self, Class):
        Classes = {Bmesh33Cell, Mesh2Tri, Mesh3Tetra}
        if isinstance(Class, type):
            if Class is Bmesh33Cell:
                return Bmesh33Cell.from_robj_ptr(self._obj)
            elif Class is Mesh2Tri:
                return Mesh2Tri.from_robj_ptr(self._obj)
            elif Class is Mesh3Tetra:
                return Mesh3Tetra.from_robj_ptr(self._obj)
            else:
                raise Exception(f'`Class` should be one of: {Classes}')
        else:
            raise Exception(f'`Class` should be an instance of type')

cdef class Isect:
    cdef isect _isect

    def __repr__(self):
        return f'Isect(t = {self.t}, obj = {self.obj})'

    @property
    def hit(self):
        return isfinite(self._isect.t)

    @property
    def t(self):
        return self._isect.t

    @property
    def obj(self):
        if np.isfinite(self.t):
            return Robj.from_ptr(self._isect.obj)

cdef class Rtree:
    '''An R-tree data structure, intended to be used for speeding up
raytracing and other basic geometric queries. Right now, it can only
be constructed from :class:`Mesh2` instances, but this limitation will
be removed in the near future.

    :param mesh: A triangle mesh stored as a :class:`Mesh2`
        instance. An R-tree will be built for this mesh.
    '''
    cdef:
        bool ptr_owner
        rtree *_rtree

    def __cinit__(self, leaf_thresh=32,
                  split_strategy=RtreeSplitStrategy.SurfaceArea,
                  should_alloc_and_init=True):
        if should_alloc_and_init:
            self.ptr_owner = True
            rtree_alloc(&self._rtree)
            rtree_init(self._rtree, leaf_thresh, split_strategy.value)

    def __dealloc__(self):
        if self.ptr_owner:
            rtree_deinit(self._rtree)
            rtree_dealloc(&self._rtree)

    @staticmethod
    cdef from_ptr(rtree *rtree_ptr, ptr_owner=False):
        rtree = Rtree(should_alloc_and_init=False)
        rtree.ptr_owner = ptr_owner
        rtree._rtree = rtree_ptr
        return rtree

    @staticmethod
    def from_mesh2(Mesh2 mesh, leaf_thresh=32,
                  split_strategy=RtreeSplitStrategy.SurfaceArea):
        rtree = Rtree(leaf_thresh, split_strategy)
        rtree_insert_mesh2(rtree._rtree, mesh.mesh)
        rtree_build(rtree._rtree)
        return rtree

    @property
    def bounding_box(self):
        return Rect3(rtree_get_bbox(self._rtree))

    @property
    def num_leaf_nodes(self):
        return rtree_get_num_leaf_nodes(self._rtree)

    def copy(self):
        cdef rtree *rtree = rtree_copy(self._rtree)
        return Rtree.from_ptr(rtree, ptr_owner=True)

    cdef insert_bmesh33(self, Bmesh33 bmesh):
        rtree_insert_bmesh33(self._rtree, bmesh.bmesh)

    cdef insert_mesh2(self, Mesh2 mesh):
        rtree_insert_mesh2(self._rtree, mesh.mesh)

    cdef insert_mesh3(self, Mesh3 mesh):
        rtree_insert_mesh3(self._rtree, mesh.mesh)

    def insert(self, obj):
        # TODO: this is a gross way to do this, but it works for
        # now. Would be nice to come up with a simpler way to dispatch
        # on the type of obj. One thing that would probably work is
        # implementing an `insert_into_rtree` method for each of the
        # classes listed below.
        Classes = {'Bmesh33', 'Mesh2', 'Mesh3'}
        if isinstance(obj, Bmesh33):
            self.insert_bmesh33(obj)
        elif isinstance(obj, Mesh2):
            self.insert_mesh2(obj)
        elif isinstance(obj, Mesh3):
            self.insert_mesh3(obj)
        else:
            raise Exception(f'obj should be an instance of: {Classes}')

    def build(self):
        rtree_build(self._rtree)

    def intersect(self, dbl[::1] org, dbl[::1] dir):
        cdef ray3 ray
        memcpy(&ray.org[0], &org[0], 3*sizeof(dbl))
        memcpy(&ray.dir[0], &dir[0], 3*sizeof(dbl))
        isect = Isect()
        rtree_intersect(self._rtree, &ray, &isect._isect)
        return isect

cdef class UpdateTetra:
    cdef utetra *_utetra

    def __cinit__(self, dbl[:] x, dbl[:, :] Xt, dbl[:] T, dbl[:, :] DT):
        cdef jet3 jet[3]
        cdef int i
        cdef int j
        for i in range(3):
            jet[i].f = T[i]
            jet[i].fx = DT[i, 0]
            jet[i].fy = DT[i, 1]
            jet[i].fz = DT[i, 2]
        cdef dbl Xt_[3][3]
        for i in range(3):
            for j in range(3):
                Xt_[i][j] = Xt[i, j]
        utetra_alloc(&self._utetra)
        utetra_init(self._utetra, &x[0], Xt_, jet)

    def __dealloc__(self):
        utetra_dealloc(&self._utetra)

    def is_degenerate(self):
        return utetra_is_degenerate(self._utetra)

    def reset(self):
        utetra_reset(self._utetra)

    def solve(self):
        utetra_solve(self._utetra)

    def get_lambda(self):
        cdef dbl[:] lam = np.empty((2,), dtype=np.float64)
        utetra_get_lambda(self._utetra, &lam[0])
        return np.asarray(lam)

    def set_lambda(self, dbl[:] lam):
        utetra_set_lambda(self._utetra, &lam[0])

    def get_value(self):
        return utetra_get_value(self._utetra)

    def get_gradient(self):
        cdef dbl[:] g = np.empty((2,), dtype=np.float64)
        utetra_get_gradient(self._utetra, &g[0])
        return np.asarray(g)

    def get_jet(self):
        cdef jet3 jet
        utetra_get_jet(self._utetra, &jet)
        return Jet3(jet.f, jet.fx, jet.fy, jet.fz)

    def get_lag_mults(self):
        cdef dbl[:] alpha = np.empty((3,), dtype=np.float64)
        utetra_get_lag_mults(self._utetra, &alpha[0])
        return np.asarray(alpha)

    def get_num_iter(self):
        return utetra_get_num_iter(self._utetra)

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

cdef class Mesh3Tetra:
    cdef mesh3_tetra *_tetra

    @staticmethod
    cdef from_robj_ptr(const robj *obj):
        tetra = Mesh3Tetra()
        tetra._tetra = <mesh3_tetra *>robj_get_data(obj)
        return tetra

    @property
    def mesh(self):
        return Mesh3.from_ptr(<mesh3 *>self._tetra.mesh)

    @property
    def index(self):
        return self._tetra.l

cdef class Mesh3:
    cdef:
        bool ptr_owner
        mesh3 *mesh
        ArrayView verts_view
        ArrayView cells_view

    def __dealloc__(self):
        if self.ptr_owner:
            mesh3_deinit(self.mesh)
            mesh3_dealloc(&self.mesh)

    @staticmethod
    def from_verts_and_cells(dbl[:, ::1] verts, size_t[:, ::1] cells,
                             compute_bd_info=True):
        mesh = Mesh3()
        mesh.ptr_owner = True
        mesh3_alloc(&mesh.mesh)
        cdef size_t nverts = verts.shape[0]
        cdef size_t ncells = cells.shape[0]
        mesh3_init(mesh.mesh, &verts[0, 0], nverts, &cells[0, 0], ncells,
                   compute_bd_info)
        mesh._set_views()
        return mesh

    @staticmethod
    cdef from_ptr(mesh3 *mesh_ptr):
        mesh = Mesh3()
        mesh.ptr_owner = False
        mesh.mesh = mesh_ptr
        mesh._set_views()
        return mesh

    cdef _set_views(self):
        self.verts_view = ArrayView(2)
        self.verts_view.readonly = True
        self.verts_view.ptr = <void *>mesh3_get_verts_ptr(self.mesh)
        self.verts_view.shape[0] = self.num_verts
        self.verts_view.shape[1] = 3
        self.verts_view.strides[0] = 4*sizeof(dbl)
        self.verts_view.strides[1] = 1*sizeof(dbl)
        self.verts_view.format = 'd'
        self.verts_view.itemsize = sizeof(dbl)

        self.cells_view = ArrayView(2)
        self.cells_view.readonly = True
        self.cells_view.ptr = <void *>mesh3_get_cells_ptr(self.mesh)
        self.cells_view.shape[0] = self.num_cells
        self.cells_view.shape[1] = 4
        self.cells_view.strides[0] = 4*sizeof(size_t)
        self.cells_view.strides[1] = 1*sizeof(size_t)
        self.cells_view.format = 'L'
        self.cells_view.itemsize = sizeof(size_t)

    @property
    def num_verts(self):
        return mesh3_nverts(self.mesh)

    @property
    def num_cells(self):
        return mesh3_ncells(self.mesh)

    @property
    def verts(self):
        return np.asarray(self.verts_view)

    @property
    def cells(self):
        return np.asarray(self.cells_view)

    def get_bbox(self):
        cdef rect3 bbox
        mesh3_get_bbox(self.mesh, &bbox)
        return ((bbox.min[0], bbox.min[1], bbox.min[2]),
                (bbox.max[0], bbox.max[1], bbox.max[2]))

    def vc(self, size_t i):
        cdef int nvc = mesh3_nvc(self.mesh, i)
        cdef size_t[::1] vc = np.empty((nvc,), dtype=np.uintp)
        mesh3_vc(self.mesh, i, &vc[0])
        return np.asarray(vc)

    def vv(self, size_t i):
        cdef int nvv = mesh3_nvv(self.mesh, i)
        cdef size_t[::1] vv = np.empty((nvv,), dtype=np.uintp)
        mesh3_vv(self.mesh, i, &vv[0])
        return np.asarray(vv)

    def cc(self, size_t i):
        cdef int ncc = mesh3_ncc(self.mesh, i)
        cdef size_t[::1] cc = np.empty((ncc,), dtype=np.uintp)
        mesh3_cc(self.mesh, i, &cc[0])
        return np.asarray(cc)

    def cv(self, size_t i):
        cdef size_t[::1] cv = np.empty((4,), dtype=np.uintp)
        mesh3_cv(self.mesh, i, &cv[0])
        return np.asarray(cv)

    def ec(self, size_t i, size_t j):
        cdef int nec = mesh3_nec(self.mesh, i, j)
        cdef size_t[::1] ec = np.empty((nec,), dtype=np.uintp)
        mesh3_ec(self.mesh, i, j, &ec[0])
        return np.asarray(ec)

    @property
    def has_bd_info(self):
        return mesh3_has_bd_info(self.mesh)

    def bdc(self, size_t i):
        if not self.has_bd_info:
            raise Exception("mesh wasn't built with boundary info");
        return mesh3_bdc(self.mesh, i)

    def bdv(self, size_t i):
        if not self.has_bd_info:
            raise Exception("mesh wasn't built with boundary info");
        return mesh3_bdv(self.mesh, i)

    def bde(self, size_t i, size_t j):
        if not self.has_bd_info:
            raise Exception("mesh wasn't built with boundary info");
        cdef size_t l[2]
        l[0] = i
        l[1] = j
        return mesh3_bde(self.mesh, l)

    def bdf(self, size_t i, size_t j, size_t k):
        if not self.has_bd_info:
            raise Exception("mesh wasn't built with boundary info");
        cdef size_t l[3]
        l[0] = i
        l[1] = j
        l[2] = k
        return mesh3_bdf(self.mesh, l)

    def is_diff_edge(self, size_t i, size_t j):
        if not self.has_bd_info:
            raise Exception("mesh wasn't built with boundary info");
        cdef size_t l[2]
        l[0] = i
        l[1] = j
        return mesh3_is_diff_edge(self.mesh, l)

    @property
    def min_tetra_alt(self):
        '''The minimum tetrahedron altitude, taken over all cells in this
instance of `jmm.Mesh3`.

        '''
        return mesh3_get_min_tetra_alt(self.mesh)

    def get_surface_mesh(self):
        if not self.has_bd_info:
            raise Exception("mesh wasn't built with boundary info");
        cdef mesh2 *surface_mesh = mesh3_get_surface_mesh(self.mesh)
        return Mesh2.from_ptr(surface_mesh, ptr_owner=True)

cdef class Jet3:
    cdef jet3 jet

    def __cinit__(self, dbl f, dbl fx, dbl fy, dbl fz):
        self.jet.f = f
        self.jet.fx = fx
        self.jet.fy = fy
        self.jet.fz = fz

    @property
    def f(self):
        return self.jet.f

    @property
    def fx(self):
        return self.jet.fx

    @property
    def fy(self):
        return self.jet.fy

    @property
    def fz(self):
        return self.jet.fz

cdef class Parent3:
    cdef par3 par

    def __cinit__(self, par3 par):
        self.par = par

    @property
    def size(self):
        return par3_size(&self.par)

    @property
    def l(self):
        cdef size_t[:] l = np.empty((self.size,), dtype=np.uintp)
        memcpy(&l[0], self.par.l, self.size*sizeof(size_t))
        return np.asarray(l)

    @property
    def b(self):
        cdef dbl[:] b = np.empty((self.size,), dtype=np.float64)
        memcpy(&b[0], self.par.b, self.size*sizeof(dbl))
        return np.asarray(b)

cdef class Cutedge:
    cdef:
        dbl t
        dbl[::1] n

    def __cinit__(self, dbl t, dbl[::1] n):
        self.t = t
        self.n = np.empty((3,), dtype=np.float64)
        self.n[:] = n[:]

    @property
    def t(self):
        return self.t

    @property
    def n(self):
        return np.asarray(self.n)

cdef class Eik3:
    cdef:
        eik3 *eik
        ArrayView jet_view
        ArrayView state_view

    def __cinit__(self, Mesh3 mesh):
        eik3_alloc(&self.eik)
        eik3_init(self.eik, mesh.mesh)

        self.jet_view = ArrayView(1)
        self.jet_view.readonly = True
        self.jet_view.ptr = <void *>eik3_get_jet_ptr(self.eik)
        self.jet_view.shape[0] = self.size
        self.jet_view.strides[0] = 4*sizeof(dbl)
        self.jet_view.format = 'dddd'
        self.jet_view.itemsize = 4*sizeof(dbl)

        self.state_view = ArrayView(1)
        self.state_view.readonly = True
        self.state_view.ptr = <void *>eik3_get_state_ptr(self.eik)
        self.state_view.shape[0] = self.size
        self.state_view.strides[0] = sizeof(state)
        self.state_view.format = 'i'
        self.state_view.itemsize = sizeof(state)

    def __dealloc__(self):
        eik3_deinit(self.eik)
        eik3_dealloc(&self.eik)

    def peek(self):
        return eik3_peek(self.eik)

    def step(self):
        return eik3_step(self.eik)

    def solve(self):
        eik3_solve(self.eik)

    def add_trial(self, size_t ind, Jet3 jet):
        eik3_add_trial(self.eik, ind, jet.jet)

    def add_valid(self, size_t ind, Jet3 jet):
        eik3_add_valid(self.eik, ind, jet.jet)

    def is_far(self, size_t ind):
        return eik3_is_far(self.eik, ind)

    def is_trial(self, size_t ind):
        return eik3_is_trial(self.eik, ind)

    def is_valid(self, size_t ind):
        return eik3_is_valid(self.eik, ind)

    def is_shadow(self, size_t ind):
        return eik3_is_shadow(self.eik, ind)

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
    def state(self):
        return np.asarray(self.state_view)

    def get_parent(self, ind):
        cdef par3 p = eik3_get_par(self.eik, ind)
        return Parent3(p)

    @property
    def shadow_cutset(self):
        cdef:
            edgemap_iter *it
            edge e
            cutedge c
        edgemap_iter_alloc(&it)
        edgemap_iter_init(it, eik3_get_cutset(self.eik))
        cutset = dict()
        while edgemap_iter_next(it, &e, &c):
            cutset[e.l[0], e.l[1]] = Cutedge(c.t, <double[:3]>c.n)
        edgemap_iter_dealloc(&it)
        return cutset

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
        return Bb33(f, Df, x)
