# cython: embedsignature=True
# cython: language_level=3

import numpy as np

import array

from cython cimport Py_buffer

from enum import Enum

ctypedef bint bool

cdef extern from "def.h":
    ctypedef double dbl
    cdef enum stype:
        CONSTANT
        NUM_STYPE
    cdef enum error:
        SUCCESS
        BAD_ARGUMENT

cdef extern from "dial.h":
    cdef struct dial3:
        pass
    void dial3_alloc(dial3 **dial)
    error dial3_init(dial3 *dial, stype stype, const int *shape, dbl h)
    void dial3_deinit(dial3 *dial)
    void dial3_dealloc(dial3 **dial)
    void dial3_add_trial(dial3 *dial, const int *ind, dbl T, const dbl *grad_T)
    void dial3_add_point_source_with_trial_nbs(dial3 *dial, const int *ind0, dbl T0)
    bool dial3_step(dial3 *dial)
    void dial3_solve(dial3 *dial)
    dbl dial3_get_T(const dial3 *dial, int l)
    dbl *dial3_get_T_ptr(const dial3 *dial)
    void dial3_get_grad_T(const dial3 *dial, int l, dbl *grad_T)

cdef class _Dial3:
    cdef:
        dial3 *dial
        Py_ssize_t shape[3]
        Py_ssize_t strides[3]

    def __cinit__(self, stype stype, int[:] shape, dbl h):
        dial3_alloc(&self.dial)
        dial3_init(self.dial, stype, &shape[0], h)
        self.shape[0] = shape[0]
        self.shape[1] = shape[1]
        self.shape[2] = shape[2]
        self.strides[2] = sizeof(dbl)
        self.strides[1] = sizeof(dbl)*self.shape[2]
        self.strides[0] = sizeof(dbl)*self.shape[2]*self.shape[1]

    def __dealloc__(self):
        dial3_deinit(self.dial)
        dial3_dealloc(&self.dial)

    def add_trial(self, int[:] ind, dbl T, dbl[:] grad_T):
        dial3_add_trial(self.dial, &ind[0], T, &grad_T[0])

    def add_point_source_with_trial_nbs(self, int[:] ind0, dbl T0):
        dial3_add_point_source_with_trial_nbs(self.dial, &ind0[0], T0)

    def step(self):
        dial3_step(self.dial)

    def solve(self):
        dial3_solve(self.dial)

    def __getbuffer__(self, Py_buffer *buf, int flags):
        buf.buf = <char *>dial3_get_T_ptr(self.dial)
        buf.format = 'd'
        buf.internal = NULL
        buf.itemsize = sizeof(dbl)
        buf.len = self.shape[0]*self.shape[1]*self.shape[2]
        buf.ndim = 3
        buf.obj = self
        buf.readonly = 1
        buf.shape = self.shape
        buf.strides = self.strides
        buf.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buf):
        pass

class Stype(Enum):
    Constant = 0

class Dial:

    def __init__(self, stype, shape, h):
        self.shape = shape
        self._dial = _Dial3(stype.value, array.array('i', [*shape]), h)

    def add_point_source_with_trial_nbs(self, ind0, T0):
        self._dial.add_point_source_with_trial_nbs(array.array('i', [*ind0]), T0)

    def step(self):
        self._dial.step()

    def solve(self):
        self._dial.solve()

    @property
    def T(self):
        return np.asarray(self._dial)
