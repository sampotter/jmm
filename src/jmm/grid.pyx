import numpy as np

cdef class Grid2:
    def __cinit__(self, int[:] shape, dbl[:] xymin, dbl h):
        self._grid.shape = [shape[0], shape[1]]
        self._grid.xymin = [xymin[0], xymin[1]]
        self._grid.h = h
        self._grid.order = ORDER_ROW_MAJOR

cdef class Grid3:
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
