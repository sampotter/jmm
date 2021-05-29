import numpy as np

def slerp(dbl[:, ::1] p, dbl[::1] w, dbl tol):
    if p.shape[0] <= 1:
        raise ValueError('slerp requires at least two input vectors')
    if p.shape[0] != w.shape[0]:
        raise ValueError('p and w should have compatible shapes')
    cdef dbl[::1] q = np.empty(3, dtype=np.float64)
    if p.shape[0] == 2:
        slerp2(&p[0][0], &p[1][0], &w[0], &q[0])
    elif p.shape[0] == 3:
        slerp3(<const dbl (*)[3]>&p[0][0], &w[0], &q[0], tol)
    elif p.shape[0] == 4:
        slerp4(<const dbl (*)[3]>&p[0][0], &w[0], &q[0], tol)
    else:
        raise RuntimeError('slerp only takes 2, 3, or 4 input vectors for now')
    return np.asarray(q)
