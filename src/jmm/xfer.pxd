from jmm.defs cimport dbl
from jmm.grid cimport grid3
from jmm.jet cimport jet3
from jmm.mesh cimport mesh3

cdef extern from "xfer.h":
    void xfer(const mesh3 *mesh, const jet3 *jet, const grid3 *grid, dbl *y)
