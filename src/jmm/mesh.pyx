import numpy as np

from libc.string cimport memcpy

from jmm.geom cimport Rect3
from jmm.rtree cimport robj, robj_get_data

cdef class Mesh2Tri:
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

cdef class Mesh3Tetra:
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
    def __dealloc__(self):
        if self.ptr_owner:
            mesh3_deinit(self.mesh)
            mesh3_dealloc(&self.mesh)

    def __init__(self, dbl[:, ::1] verts, size_t[:, ::1] cells,
                 bool compute_bd_info=True):
        self.ptr_owner = True
        mesh3_alloc(&self.mesh)
        cdef size_t nverts = verts.shape[0]
        cdef size_t ncells = cells.shape[0]
        mesh3_init(self.mesh, &verts[0, 0], nverts, &cells[0, 0], ncells,
                   compute_bd_info)
        self._set_views()

    @staticmethod
    def from_verts_and_cells(dbl[:, ::1] verts, size_t[:, ::1] cells,
                             compute_bd_info=True):
        return Mesh3(verts, cells, compute_bd_info)

    @staticmethod
    cdef from_ptr(mesh3 *mesh_ptr):
        cdef Mesh3 mesh = Mesh3.__new__(Mesh3)
        mesh.ptr_owner = False
        mesh.mesh = mesh_ptr
        mesh._set_views()
        return mesh

    def __reduce__(self):
        # TODO: each time we do this, we rebuild the mesh and
        # recompute all the boundary information. We should be able to
        # do this faster by reaching into mesh3 and pulling out that
        # info. This is just a quick fix.

        cdef dbl[:, ::1] verts = np.empty((self.num_verts, 3), dtype=np.float64)
        memcpy(&verts[0, 0], mesh3_get_verts_ptr(self.mesh),
               3*self.num_verts*sizeof(dbl))

        cdef size_t[:, ::1] cells = np.empty((self.num_cells, 4), dtype=np.uintp)
        memcpy(&cells[0, 0], mesh3_get_cells_ptr(self.mesh),
               4*self.num_cells*sizeof(size_t))

        cdef bool compute_bd_info = mesh3_has_bd_info(self.mesh)

        return (
            self.__class__,
            (np.asarray(verts), np.asarray(cells), compute_bd_info)
        )

    cdef _set_views(self):
        self.verts_view = ArrayView(2)
        self.verts_view.readonly = True
        self.verts_view.ptr = <void *>mesh3_get_verts_ptr(self.mesh)
        self.verts_view.shape[0] = self.num_verts
        self.verts_view.shape[1] = 3
        self.verts_view.strides[0] = 3*sizeof(dbl)
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

    def get_diff_edges(self):
        cdef size_t le[2]
        for l in range(mesh3_nbde(self.mesh)):
            mesh3_get_bde_inds(self.mesh, l, le)
            if mesh3_is_diff_edge(self.mesh, le):
                yield np.array([le[0], le[1]], dtype=np.uintp)

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

    @property
    def min_edge_length(self):
        '''The minimum edge length of the tetrahedron mesh, taken over all
cell edges of `jmm.Mesh3`.

        '''
        return mesh3_get_min_edge_length(self.mesh)

    def get_surface_mesh(self):
        if not self.has_bd_info:
            raise Exception("mesh wasn't built with boundary info");
        cdef mesh2 *surface_mesh = mesh3_get_surface_mesh(self.mesh)
        return Mesh2.from_ptr(surface_mesh, ptr_owner=True)

    @property
    def num_reflectors(self):
        '''The number of distinct reflecting components that partition the
        surface of the mesh.

        '''
        return mesh3_get_num_reflectors(self.mesh)

    def get_reflector(self, size_t i):
        '''Get the faces that make up the `i`th reflector.'''
        cdef size_t num_faces = mesh3_get_reflector_size(self.mesh, i)
        cdef size_t[:, ::1] lf = np.empty((num_faces, 3), dtype=np.uintp)
        mesh3_get_reflector(self.mesh, i, <size_t[3]*>&lf[0, 0])
        return np.asarray(lf)

    def get_active_reflectors(self, shadow_mask):
        '''Return the reflectors which contain any points that don't lie in
        the shadow indicated by `shadow_mask`.

        '''
        for i in range(self.num_reflectors):
            faces = self.get_reflector(i)
            refl_shadow_mask = np.logical_not(shadow_mask[faces]).any(1)
            if refl_shadow_mask.any():
                yield i, faces

    @property
    def num_diffractors(self):
        '''The number of distinct diffracting components. For a `Mesh3`, these
        are maximal straight line components embedded in the boundary.

        '''
        return mesh3_get_num_diffractors(self.mesh)

    def get_diffractor(self, size_t i):
        '''Get the edges that make up the `i`th diffractor.'''
        cdef size_t num_edges = mesh3_get_diffractor_size(self.mesh, i)
        cdef size_t[:, ::1] le = np.empty((num_edges, 2), dtype=np.uintp)
        mesh3_get_diffractor(self.mesh, i, <size_t[2]*>&le[0, 0])
        return np.asarray(le)

    def get_active_diffractors(self, shadow_mask):
        '''Return the diffractors that contain any points that don't lie in
        the shadow specified by `shadow_mask`.

        '''
        for i in range(self.num_diffractors):
            edges = self.get_diffractor(i)
            restricted_mask = np.logical_not(shadow_mask[edges]).any(1)
            if restricted_mask.any():
                yield i, edges

    def set_boundary_edge(self, size_t i, size_t j, bool diff):
        cdef size_t le[2]
        le[0] = i
        le[1] = j
        mesh3_set_bde(self.mesh, le, diff)

    @property
    def eps(self):
        return mesh3_get_eps(self.mesh)

    def get_face_normal(self, size_t l0, size_t l1, size_t l2):
        if not self.bdf(l0, l1, l2):
            raise KeyError('(%d, %d, %d) is not a face' % sorted([l0, l1, l2]))
        cdef size_t[3] lf = [l0, l1, l2]
        cdef dbl[::1] normal = np.empty((3,), dtype=np.float64)
        mesh3_get_face_normal(self.mesh, lf, &normal[0])
        return np.asarray(normal)
