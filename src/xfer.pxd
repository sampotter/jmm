from defs cimport dbl
from grid3 cimport *
from jet cimport *
from mesh3 cimport *

cdef extern from "xfer.h":
    void xfer(const mesh3 *mesh, const jet3 *jet, const grid3 *grid, dbl *y)
