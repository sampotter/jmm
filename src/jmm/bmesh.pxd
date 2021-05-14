from jmm.bb cimport bb33
from jmm.defs cimport bool, dbl
from jmm.geom cimport ray3
from jmm.jet cimport jet3
from jmm.mesh cimport mesh3
from jmm.rtree cimport robj

cdef extern from "bmesh.h":
    cdef struct bmesh33_cell:
        const bb33 *bb
        const mesh3 *mesh
        dbl level
        size_t l

    bool bmesh33_cell_intersect(const bmesh33_cell *cell, const ray3 *ray, dbl *t)

    cdef struct bmesh33:
        pass

    void bmesh33_alloc(bmesh33 **bmesh)
    void bmesh33_dealloc(bmesh33 **bmesh)
    void bmesh33_init_from_mesh3_and_jets(bmesh33 *bmesh, const mesh3 *mesh,
                                          const jet3 *jet)
    void bmesh33_deinit(bmesh33 *bmesh)
    size_t bmesh33_num_cells(const bmesh33 *bmesh)
    const mesh3 *bmesh33_get_mesh_ptr(const bmesh33 *bmesh)
    bmesh33 *bmesh33_restrict_to_level(const bmesh33 *bmesh, dbl level)
    bmesh33_cell bmesh33_get_cell(const bmesh33 *bmesh, size_t l)

cdef class Bmesh33Cell:
    cdef bmesh33_cell cell

    @staticmethod
    cdef from_robj_ptr(const robj *obj)

    @staticmethod
    cdef from_cell(bmesh33_cell cell)

cdef class Bmesh33:
    cdef:
        bool ptr_owner
        bmesh33 *bmesh

    @staticmethod
    cdef from_ptr(bmesh33 *bmesh_ptr, ptr_owner=?)
