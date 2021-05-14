cdef extern from "edge.h":
    cdef struct edge:
        size_t l[2]
    edge make_edge(size_t l0, size_t l1)
