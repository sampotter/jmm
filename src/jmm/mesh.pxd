from jmm.array_view cimport ArrayView
from jmm.defs cimport bool, dbl, dbl2, dbl3, uint3
from jmm.geom cimport rect3, tri3
from jmm.rtree cimport robj

cdef extern from "mesh2.h":
    cdef struct mesh2:
        pass

    cdef struct mesh2_tri:
        mesh2 *mesh
        size_t l

    void mesh2_alloc(mesh2 **mesh)
    void mesh2_dealloc(mesh2 **mesh)
    void mesh2_init(mesh2 *mesh, const dbl *verts, size_t nverts,
                    const size_t *faces, size_t nfaces)
    void mesh2_deinit(mesh2 *mesh)
    size_t mesh2_get_num_points(const mesh2 *mesh)
    dbl *mesh2_get_points_ptr(const mesh2 *mesh)
    size_t mesh2_get_num_faces(const mesh2 *mesh)
    size_t *mesh2_get_faces_ptr(const mesh2 *mesh)
    rect3 mesh2_get_bounding_box(const mesh2 *mesh)
    void mesh2_get_centroid(const mesh2 *mesh, size_t i, dbl *centroid)
    void mesh2_get_vertex(const mesh2 *mesh, size_t i, size_t j, dbl *v)
    bool mesh2_tri_bbox_overlap(const mesh2 *mesh, size_t i, const rect3 *bbox)
    tri3 mesh2_get_tri(const mesh2 *mesh, size_t i)

cdef class Mesh2Tri:
    cdef mesh2_tri *_tri

    @staticmethod
    cdef from_robj_ptr(const robj *obj)

cdef class Mesh2:
    cdef:
        bool ptr_owner
        mesh2 *mesh
        ArrayView verts_view
        ArrayView faces_view

    @staticmethod
    cdef from_ptr(const mesh2 *mesh_ptr, ptr_owner=?)

    cdef _set_views(self)

cdef extern from "mesh22.h":
    cdef struct mesh22:
        pass

    void mesh22_alloc(mesh22 **mesh)
    void mesh22_dealloc(mesh22 **mesh)
    void mesh22_init(mesh22 *mesh, const dbl2 *verts, size_t nverts,
                     const uint3 *faces, size_t nfaces)
    void mesh22_deinit(mesh22 *mesh)
    const dbl2 *mesh22_get_verts_ptr(const mesh22 *mesh)
    const uint3 *mesh22_get_faces_ptr(const mesh22 *mesh)
    size_t mesh22_nverts(const mesh22 *mesh)
    size_t mesh22_nfaces(const mesh22 *mesh)
    void mesh22_get_vert(const mesh22 *mesh, size_t l, dbl2 x)
    size_t mesh22_nvf(const mesh22 *mesh, size_t l)
    void mesh22_vf(const mesh22 *mesh, size_t l, size_t *vf)
    size_t mesh22_nvv(const mesh22 *mesh, size_t l)
    void mesh22_vv(const mesh22 *mesh, size_t l, size_t *vv)

cdef class Mesh22:
    cdef:
        bool ptr_owner
        mesh22 *mesh
        ArrayView verts_view
        ArrayView faces_view

    @staticmethod
    cdef from_ptr(mesh22 *mesh_ptr, ptr_owner=?)

    cdef _set_views(self)

cdef extern from "mesh3.h":
    cdef struct mesh3:
        pass

    cdef struct mesh3_tetra:
        mesh3 *mesh
        size_t l

    void mesh3_alloc(mesh3 **mesh)
    void mesh3_dealloc(mesh3 **mesh)
    void mesh3_init(mesh3 *mesh,
                    dbl *verts, size_t nverts,
                    size_t *cells, size_t ncells,
                    bool compute_bd_info, const dbl *eps)
    void mesh3_deinit(mesh3 *mesh)
    void mesh3_get_bbox(const mesh3 *mesh, rect3 *bbox)
    const size_t *mesh3_get_cells_ptr(const mesh3 *mesh)
    const dbl *mesh3_get_verts_ptr(const mesh3 *mesh)
    size_t mesh3_ncells(const mesh3 *mesh)
    size_t mesh3_nverts(const mesh3 *mesh)
    int mesh3_nvc(const mesh3 *mesh, size_t i)
    void mesh3_vc(const mesh3 *mesh, size_t i, size_t *vc)
    int mesh3_nvv(const mesh3 *mesh, size_t i)
    void mesh3_vv(const mesh3 *mesh, size_t i, size_t *vv)
    int mesh3_ncc(const mesh3 *mesh, size_t i)
    void mesh3_cc(const mesh3 *mesh, size_t i, size_t *cc)
    void mesh3_cv(const mesh3 *mesh, size_t i, size_t *cv)
    int mesh3_nec(const mesh3 *mesh, size_t i, size_t j)
    void mesh3_ec(const mesh3 *mesh, size_t i, size_t j, size_t *ec)
    bool mesh3_has_bd_info(const mesh3 *mesh)
    bool mesh3_bdc(const mesh3 *mesh, size_t i)
    bool mesh3_bdv(const mesh3 *mesh, size_t i)
    size_t mesh3_nbde(const mesh3 *mesh)
    bool mesh3_bde(const mesh3 *mesh, const size_t l[2])
    size_t mesh3_nbdf(const mesh3 *mesh)
    void mesh3_get_bdf_inds(const mesh3 *mesh, size_t l, size_t lf[3])
    void mesh3_set_bdf(mesh3 *mesh, const size_t lf[3])
    bool mesh3_is_bdf(const mesh3 *mesh, const size_t l[3])
    bool mesh3_is_diff_edge(const mesh3 *mesh, const size_t l[2])
    dbl mesh3_get_min_tetra_alt(const mesh3 *mesh)
    dbl mesh3_get_min_edge_length(const mesh3 *mesh)
    mesh2 *mesh3_get_surface_mesh(const mesh3 *mesh)
    size_t mesh3_get_num_inc_diff_edges(const mesh3 *mesh, size_t l)
    void mesh3_get_inc_diff_edges(const mesh3 *mesh, size_t l, size_t (*le)[2])
    size_t mesh3_get_num_reflectors(const mesh3 *mesh)
    size_t mesh3_get_reflector_size(const mesh3 *mesh, size_t i)
    void mesh3_get_reflector(const mesh3 *mesh, size_t i, size_t (*lf)[3])
    size_t mesh3_get_num_diffractors(const mesh3 *mesh)
    size_t mesh3_get_diffractor_size(const mesh3 *mesh, size_t i)
    void mesh3_get_diffractor(const mesh3 *mesh, size_t i, size_t (*le)[2])
    void mesh3_get_bde_inds(const mesh3 *mesh, size_t l, size_t le[2])
    void mesh3_set_bde(mesh3 *mesh, const size_t le[2], bool diff)
    dbl mesh3_get_eps(const mesh3 *mesh)
    void mesh3_get_face_normal(const mesh3 *mesh, const size_t lf[3], dbl normal[3])
    dbl mesh3_get_edge_ext_angle(const mesh3 *mesh, const size_t le[2])

cdef class Mesh3Tetra:
    cdef mesh3_tetra *_tetra

    @staticmethod
    cdef from_robj_ptr(const robj *obj)

cdef class Mesh3:
    cdef:
        bool ptr_owner
        mesh3 *mesh
        ArrayView verts_view
        ArrayView cells_view

    @staticmethod
    cdef from_ptr(mesh3 *mesh_ptr)

    cdef _set_views(self)
