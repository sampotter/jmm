from mesh2 cimport mesh2_alloc

cdef class Mesh2:
    def __cinit__(self, const char *verts_path, const char *faces_path):
        mesh2_alloc(&self._mesh)
        mesh2_init_from_binary_files(self._mesh, verts_path, faces_path)

    def __dealloc__(self):
        mesh2_deinit(self._mesh)
        mesh2_dealloc(&self._mesh)

    @property
    def num_points(self):
        return mesh2_get_num_points(self._mesh)

    @property
    def num_faces(self):
        return mesh2_get_num_faces(self._mesh)

    @property
    def bounding_box(self):
        return Rect3(mesh2_get_bounding_box(self._mesh))

    cdef get_mesh_ptr(self):
        return self._mesh
