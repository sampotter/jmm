from libc.stdlib cimport free, malloc

from jmm.eik cimport Eik3, eik3_get_jet_ptr, eik3_get_mesh
from jmm.mesh cimport Mesh3
from jmm.rtree cimport robj, robj_get_data

cdef class Bmesh33Cell:
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
