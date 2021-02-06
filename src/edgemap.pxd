from defs cimport bool

from edge cimport *

cdef extern from "edgemap.h":
    cdef struct edgemap_iter:
        pass
    cdef struct edgemap:
        pass
    void edgemap_iter_alloc(edgemap_iter **iter)
    void edgemap_iter_dealloc(edgemap_iter **iter)
    void edgemap_iter_init(edgemap_iter *iter, const edgemap *edgemap)
    bool edgemap_iter_next(edgemap_iter *iter, edge *edge, void *elt)
