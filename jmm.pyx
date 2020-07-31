# cython: embedsignature=True
# cython: language_level=3

import numpy as np

import array

from cython cimport Py_buffer

from libc.stdlib cimport free, malloc

from enum import Enum

ctypedef bint bool

cdef extern from "def.h":
    ctypedef double dbl
    cdef enum state:
        FAR
        TRIAL
        VALID
        BOUNDARY
        NEW_VALID
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
    state *dial3_get_state_ptr(const dial3 *dial)

cdef class ArrayView:
    cdef:
        bool readonly
        int ndim
        void *ptr
        Py_ssize_t *shape
        Py_ssize_t *strides
        char *format
        size_t itemsize

    def __cinit__(self, int ndim):
        self.ndim = ndim
        self.shape = <Py_ssize_t *>malloc(sizeof(Py_ssize_t)*ndim)
        self.strides = <Py_ssize_t *> malloc(sizeof(Py_ssize_t)*ndim)

    def __dealloc__(self):
        free(self.shape)
        free(self.strides)

    def __getbuffer__(self, Py_buffer *buf, int flags):
        buf.buf = <char *>self.ptr
        buf.format = self.format
        buf.internal = NULL
        buf.itemsize = self.itemsize
        buf.len = self.size
        buf.ndim = self.ndim
        buf.obj = self
        buf.readonly = self.readonly
        buf.shape = self.shape
        buf.strides = self.strides
        buf.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buf):
        pass

    @property
    def size(self):
        cdef Py_ssize_t size = 1
        for i in range(self.ndim):
            size *= self.shape[i]
        return size

cdef class _Dial3:
    cdef:
        dial3 *dial
        Py_ssize_t shape[3]
        ArrayView T_view
        ArrayView state_view

    def __cinit__(self, stype stype, int[:] shape, dbl h):
        dial3_alloc(&self.dial)
        dial3_init(self.dial, stype, &shape[0], h)

        self.shape[0] = shape[0]
        self.shape[1] = shape[1]
        self.shape[2] = shape[2]

        # Strides that haven't been scaled by the size of the
        # underlying type
        cdef Py_ssize_t base_strides[3]
        base_strides[2] = 1
        base_strides[1] = self.shape[2]
        base_strides[0] = self.shape[2]*self.shape[1]

        self.T_view = ArrayView(3)
        self.T_view.readonly = False
        self.T_view.ptr = <void *>dial3_get_T_ptr(self.dial)
        self.T_view.shape[0] = self.shape[0]
        self.T_view.shape[1] = self.shape[1]
        self.T_view.shape[2] = self.shape[2]
        self.T_view.strides[0] = sizeof(dbl)*base_strides[0]
        self.T_view.strides[1] = sizeof(dbl)*base_strides[1]
        self.T_view.strides[2] = sizeof(dbl)*base_strides[2]
        self.T_view.format = 'd'
        self.T_view.itemsize = sizeof(dbl)

        self.state_view = ArrayView(3)
        self.state_view.readonly = False
        self.state_view.ptr = <void *>dial3_get_state_ptr(self.dial)
        self.state_view.shape[0] = self.shape[0]
        self.state_view.shape[1] = self.shape[1]
        self.state_view.shape[2] = self.shape[2]
        self.state_view.strides[0] = sizeof(state)*base_strides[0]
        self.state_view.strides[1] = sizeof(state)*base_strides[1]
        self.state_view.strides[2] = sizeof(state)*base_strides[2]
        self.state_view.format = 'i'
        self.state_view.itemsize = sizeof(state)

        print(self.state_view.itemsize)


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

    @property
    def T(self):
        return self.T_view

    @property
    def state(self):
        return self.state_view

class Stype(Enum):
    Constant = 0

class State(Enum):
    Far = 0
    Trial = 1
    Valid = 2
    Boundary = 3
    NewValid = 4

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
        return np.asarray(self._dial.T)

    @property
    def state(self):
        return np.asarray(self._dial.state)
