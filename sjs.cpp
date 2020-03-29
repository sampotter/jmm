#include <pybind11/pybind11.h>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "bicubic.h"
#include "eik.h"
#include "eik_F3.h"
#include "eik_F4.h"
#include "field.h"
#include "heap.h"
#include "hybrid.h"
#include "index.h"
#include "jet.h"
#include "vec.h"

static dbl value_wrapper(void * vp, int l);
static void setpos_wrapper(void * vp, int l, int pos);

struct heap_wrapper
{
  heap * ptr {nullptr};
  bool should_call_dtor {false};

  std::optional<std::function<dbl(int)>> value;
  std::optional<std::function<void(int, int)>> setpos;

  heap_wrapper(heap * ptr): ptr {ptr} {}

  heap_wrapper(int capacity, decltype(value) value, decltype(setpos) setpos):
    should_call_dtor {true},
    value {value},
    setpos {setpos}
  {
    heap_alloc(&ptr);
    heap_init(
      ptr,
      capacity,
      value_wrapper,
      setpos_wrapper,
      (void *) this
    );
  }

  ~heap_wrapper() {
    if (should_call_dtor) {
      heap_deinit(ptr);
      heap_dealloc(&ptr);
    }
  }
};

static dbl value_wrapper(void * vp, int l) {
  heap_wrapper * hwp = (heap_wrapper *) vp;
  if (!hwp->value) {
    throw std::runtime_error {"ERROR: No value function for heap!"};
  }
  return (*hwp->value)(l);
}

static void setpos_wrapper(void * vp, int l, int pos) {
  heap_wrapper * hwp = (heap_wrapper *) vp;
  if (!hwp->setpos) {
    throw std::runtime_error {"ERROR: No setpos function for heap!"};
  }
  (*hwp->setpos)(l, pos);
}

static dbl f_wrapper(dbl t, void * context) {
  return (*(std::function<dbl(dbl)> *) context)(t);
}

dbl field_f_wrapper(dbl x, dbl y, void *context);
dvec2 field_grad_f_wrapper(dbl x, dbl y, void *context);

struct field2_wrapper {
  field2 field;

  std::function<dbl(dbl, dbl)> f;
  std::function<std::array<dbl, 2>(dbl, dbl)> grad_f;

  field2_wrapper(decltype(f) const & f, decltype(grad_f) const & grad_f):
    field {field_f_wrapper, field_grad_f_wrapper, (void *)this},
    f {f},
    grad_f {grad_f}
  {}
};

dbl field_f_wrapper(dbl x, dbl y, void *context) {
  field2_wrapper * wrap = (field2_wrapper *) context;
  return wrap->f(x, y);
}

dvec2 field_grad_f_wrapper(dbl x, dbl y, void *context) {
  field2_wrapper * wrap = (field2_wrapper *) context;
  dvec2 tmp;
  std::tie(tmp.x, tmp.y) = wrap->grad_f(x, y);
  return tmp;
}

struct eik_wrapper {
  field2_wrapper slow;
  eik * ptr {nullptr};

  eik_wrapper(field2_wrapper const & slow,
              std::array<int, 2> const & shape,
              std::array<dbl, 2> const & xymin,
              dbl h):
    slow {slow}
  {
    eik_alloc(&ptr);
    eik_init(
      ptr,
      &slow.field,
      ivec2 {shape[0], shape[1]},
      dvec2 {xymin[0], xymin[1]},
      h
    );
  }

  ~eik_wrapper() {
    eik_deinit(ptr);
    eik_dealloc(&ptr);
  }
};

PYBIND11_MODULE (_sjs, m) {
  m.doc() = R"pbdoc(
_sjs
----

TODO!
)pbdoc";

  // bicubic.h

  py::enum_<bicubic_variable>(m, "BicubicVariable")
    .value("Lambda", bicubic_variable::LAMBDA)
    .value("Mu", bicubic_variable::MU)
    ;

  py::class_<bicubic>(m, "Bicubic")
    .def(py::init(
           [] (std::array<std::array<dbl, 4>, 4> const & data) {
             auto ptr = std::make_unique<bicubic>();
             dmat44 data_;
             for (int i = 0; i < 4; ++i) {
               for (int j = 0; j < 4; ++j) {
                 data_.rows[i].data[j] = data[i][j];
               }
             }
             bicubic_set_data(ptr.get(), data_);
             return ptr;
           }
         ))
    .def_property_readonly(
      "A",
      [] (bicubic const & B) {
        std::array<std::array<dbl, 4>, 4> A;
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            A[i][j] = B.A.rows[i].data[j];
          }
        }
        return A;
      }
    )
    .def(
      "set_data",
      [] (bicubic & B, std::array<std::array<dbl, 4>, 4> const & data) {
        dmat44 data_;
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            data_.rows[i].data[j] = data[i][j];
          }
        }
        bicubic_set_data(&B, data_);
      }
    )
    .def(
      "get_f_on_edge",
      [] (bicubic const & B, bicubic_variable var, int edge) {
        return bicubic_get_f_on_edge(&B, var, edge);
      }
    )
    .def(
      "get_fx_on_edge",
      [] (bicubic const & B, bicubic_variable var, int edge) {
        return bicubic_get_fx_on_edge(&B, var, edge);
      }
    )
    .def(
      "get_fy_on_edge",
      [] (bicubic const & B, bicubic_variable var, int edge) {
        return bicubic_get_fy_on_edge(&B, var, edge);
      }
    )
    .def(
      "f",
      [] (bicubic const & B, dbl lambda, dbl mu) {
        return bicubic_f(&B, dvec2 {lambda, mu});
      }
    )
    .def(
      "fx",
      [] (bicubic const & B, dbl lambda, dbl mu) {
        return bicubic_fx(&B, dvec2 {lambda, mu});
      }
    )
    .def(
      "fy",
      [] (bicubic const & B, dbl lambda, dbl mu) {
        return bicubic_fy(&B, dvec2 {lambda, mu});
      }
    )
    .def(
      "fxy",
      [] (bicubic const & B, dbl lambda, dbl mu) {
        return bicubic_fxy(&B, dvec2 {lambda, mu});
      }
    )
    ;

  m.def(
    "interpolate_fxy_at_verts",
    [] (std::array<dbl, 4> const & fx, std::array<dbl, 4> const & fy, dbl h) {
      dvec4 fx_ = {.data = {fx[0], fx[1], fx[2], fx[3]}};
      dvec4 fy_ = {.data = {fy[0], fy[1], fy[2], fy[3]}};
      return interpolate_fxy_at_verts(fx_, fy_, h);
    }
  );

  // cubic.h

  py::class_<cubic>(m, "Cubic")
    .def(py::init(
           [] (std::array<dbl, 4> const & data) {
             auto ptr = std::make_unique<cubic>();
             cubic_set_data_from_ptr(ptr.get(), &data[0]);
             return ptr;
           }
         ))
    .def(
      "set_data",
      [] (cubic & C, std::array<dbl, 4> const & data) {
        cubic_set_data_from_ptr(&C, &data[0]);
      }
    )
    .def(
      "reverse_on_unit_interval",
      [] (cubic & C) {
        cubic_reverse_on_unit_interval(&C);
      }
    )
    .def_property_readonly(
      "a",
      [] (cubic const & C) {
        return C.a;
      }
    )
    .def(
      "f",
      [] (cubic const & C, dbl lam) { return cubic_f(&C, lam); }
    )
    .def(
      "df",
      [] (cubic const & C, dbl lam) { return cubic_df(&C, lam); }
    )
    ;

  // def.h

  py::enum_<state>(m, "State")
    .value("Far", state::FAR)
    .value("Trial", state::TRIAL)
    .value("Valid", state::VALID)
    .value("Boundary", state::BOUNDARY)
    ;

  // eik.h

  py::class_<eik_wrapper>(m, "Eik")
    .def(py::init<
           field2_wrapper const &,
           std::array<int, 2> const &,
           std::array<dbl, 2> const &,
           dbl
         >())
    .def(
      "step",
      [] (eik_wrapper const & w) { eik_step(w.ptr); }
    )
    .def(
      "solve",
      [] (eik_wrapper const & w) { eik_solve(w.ptr); }
    )
    .def(
      "add_trial",
      [] (eik_wrapper const & w, int i, int j, jet_s jet) {
        eik_add_trial(w.ptr, ivec2 {i, j}, jet);
      }
    )
    .def(
      "add_valid",
      [] (eik_wrapper const & w, int i, int j, jet_s jet) {
        eik_add_valid(w.ptr, ivec2 {i, j}, jet);
      }
    )
    .def(
      "make_bd",
      [] (eik_wrapper const & w, int i, int j) {
        eik_make_bd(w.ptr, ivec2 {i, j});
      }
    )
    .def_readonly("slow", &eik_wrapper::slow)
    .def_property_readonly(
      "shape",
      [] (eik_wrapper const &w) { return eik_get_shape(w.ptr); }
    )
    .def(
      "get_jet",
      [] (eik_wrapper const & w, int i, int j) {
        return eik_get_jet(w.ptr, {i, j});
      }
    )
    .def(
      "get_state",
      [] (eik_wrapper const & w, int i, int j) {
        return eik_get_state(w.ptr, {i, j});
      }
    )
    .def_property_readonly(
      "states",
      [] (eik_wrapper const & w) {
        ivec2 shape = eik_get_shape(w.ptr);
        return py::array {
          {shape.i, shape.j},
          {sizeof(state)*shape.j, sizeof(state)},
          eik_get_states_ptr(w.ptr)
        };
      }
    )
    .def(
      "T",
      [] (eik_wrapper const & w, dbl x, dbl y) {
        return eik_T(w.ptr, dvec2 {x, y});
      }
    )
    .def(
      "Tx",
      [] (eik_wrapper const & w, dbl x, dbl y) {
        return eik_Tx(w.ptr, dvec2 {x, y});
      }
    )
    .def(
      "Ty",
      [] (eik_wrapper const & w, dbl x, dbl y) {
        return eik_Ty(w.ptr, dvec2 {x, y});
      }
    )
    .def(
      "Txy",
      [] (eik_wrapper const & w, dbl x, dbl y) {
        return eik_Txy(w.ptr, dvec2 {x, y});
      }
    )
    .def(
      "can_build_cell",
      [] (eik_wrapper const & w, int i, int j) {
        return eik_can_build_cell(w.ptr, ivec2 {i, j});
      }
    )
    .def(
      "build_cells",
      [] (eik_wrapper const & w) { eik_build_cells(w.ptr); }
    )
    .def(
      "get_bicubic",
      [] (eik_wrapper const & w, int i, int j) {
        return eik_get_bicubic(w.ptr, ivec2 {i, j});
      }
    )
    .def_property_readonly(
      "bicubics",
      [] (eik_wrapper const & w) {
        ivec2 shape = eik_get_shape(w.ptr);
        return py::array {
          {shape.i - 1, shape.j - 1},
          {sizeof(bicubic)*(shape.j - 1), sizeof(bicubic)},
          eik_get_bicubics_ptr(w.ptr)
        };
      }
    )
    .def_property_readonly(
      "heap",
      [] (eik_wrapper const & w) {
        return heap_wrapper {eik_get_heap(w.ptr)};
      }
    )
    ;

  // field.h

  py::class_<field2_wrapper>(m, "Field2")
    .def(py::init<
           std::function<dbl(dbl, dbl)> const &,
           std::function<std::array<dbl, 2>(dbl, dbl)> const &
         >())
    .def(
      "s",
      [] (field2_wrapper const & wrap, dbl x, dbl y) {
        return field2_f(&wrap.field, {x, y});
      }
    )
    .def(
      "grad_s",
      [] (field2_wrapper const & wrap, dbl x, dbl y) {
        return field2_grad_f(&wrap.field, {x, y});
      }
    )
    ;

  // heap.h

  py::class_<heap_wrapper>(m, "Heap")
    .def(py::init<int, std::function<dbl(int)>, std::function<void(int,int)>>())
    .def(
      "insert",
      [] (heap_wrapper & w, int ind) { heap_insert(w.ptr, ind); }
    )
    .def(
      "swim",
      [] (heap_wrapper & w, int ind) { heap_swim(w.ptr, ind); }
    )
    .def_property_readonly(
      "front",
      [] (heap_wrapper const & w) {
        std::optional<int> l0;
        if (heap_size(w.ptr) > 0) {
          *l0 = heap_front(w.ptr);
        }
        return l0;
      }
    )
    .def(
      "pop",
      [] (heap_wrapper & w) { heap_pop(w.ptr); }
    )
    .def_property_readonly(
      "size",
      [] (heap_wrapper const & w) { return heap_size(w.ptr); }
    )
    ;

  // hybrid.h

  m.def(
    "hybrid",
    [] (std::function<dbl(dbl)> const & f, dbl a, dbl b) {
      void * context = (void *) &f;
      return hybrid(f_wrapper, a, b, context);
    }
  );

  // index.h

  m.def(
    "_ind2l",
    [] (std::array<int, 2> shape, std::array<int, 2> ind) {
      return ind2l({shape[0], shape[1]}, {ind[0], ind[1]});
    }
  );

  m.def(
    "_ind2lc",
    [] (std::array<int, 2> shape, std::array<int, 2> ind) {
      return ind2lc({shape[0], shape[1]}, {ind[0], ind[1]});
    }
  );

  m.def(
    "_indc2l",
    [] (std::array<int, 2> shape, std::array<int, 2> indc) {
      return indc2l({shape[0], shape[1]}, {indc[0], indc[1]});
    }
  );

  m.def(
    "_indc2lc",
    [] (std::array<int, 2> shape, std::array<int, 2> indc) {
      return indc2lc({shape[0], shape[1]}, {indc[0], indc[1]});
    }
  );

  m.def(
    "_l2ind",
    [] (std::array<int, 2> shape, int l) {
      return l2ind({shape[0], shape[1]}, l);
    }
  );

  m.def(
    "_l2indc",
    [] (std::array<int, 2> shape, int l) {
      return l2indc({shape[0], shape[1]}, l);
    }
  );

  m.def(
    "_lc2ind",
    [] (std::array<int, 2> shape, int lc) {
      return lc2ind({shape[0], shape[1]}, lc);
    }
  );

  m.def(
    "_lc2indc",
    [] (std::array<int, 2> shape, int lc) {
      return lc2indc({shape[0], shape[1]}, lc);
    }
  );

  m.def(
    "_l2lc",
    [] (std::array<int, 2> shape, int l) {
      return l2lc({shape[0], shape[1]}, l);
    }
  );

  m.def(
    "_lc2l",
    [] (std::array<int, 2> shape, int lc) {
      return lc2l({shape[0], shape[1]}, lc);
    }
  );

  m.def(
    "_xy_to_lc_and_cc",
    [] (std::array<int, 2> shape, std::array<dbl, 2> xymin, dbl h,
        std::array<dbl, 2> xy) {
      dvec2 cc;
      int lc = xy_to_lc_and_cc(
        {shape[0], shape[1]},
        {xymin[0], xymin[1]},
        h,
        {xy[0], xy[1]},
        &cc
      );
      return std::make_pair(lc, cc);
    }
  );

  // jet.h

  py::class_<jet>(m, "Jet")
    .def(py::init<dbl, dbl, dbl, dbl>())
    .def_readwrite("f", &jet::f)
    .def_readwrite("fx", &jet::fx)
    .def_readwrite("fy", &jet::fy)
    .def_readwrite("fxy", &jet::fxy)
    ;

  // mat.h

  py::class_<dmat22>(m, "Dmat22")
    .def(py::init(
           [] (std::array<std::array<dbl, 2>, 2> const & data) {
             auto ptr = std::make_unique<dmat22>();
             for (int i = 0; i < 2; ++i) {
               for (int j = 0; j < 2; ++j) {
                 ptr->data[i][j] = data[i][j];
               }
             }
             return ptr;
           }
         ))
    .def(
      "__getitem__",
      [] (dmat22 const & A, std::pair<int, int> ij) {
        return A.data[ij.first][ij.second];
      }
    )
    .def(
      "__matmul__",
      [] (dmat22 const & A, dvec2 const & x) {
        return dmat22_dvec2_mul(A, x);
      }
    )
    ;

  py::class_<dmat44>(m, "Dmat44")
    .def(py::init(
           [] (std::array<std::array<dbl, 4>, 4> const & data) {
             auto ptr = std::make_unique<dmat44>();
             for (int i = 0; i < 4; ++i) {
               for (int j = 0; j < 4; ++j) {
                 ptr->rows[i].data[j] = data[i][j];
               }
             }
             return ptr;
           }
         ))
    .def(py::init(
           [] (std::array<dvec4, 4> const & data) {
             auto ptr = std::make_unique<dmat44>();
             for (int i = 0; i < 4; ++i) {
               ptr->rows[i] = data[i];
             }
             return ptr;
           }
         ))
    .def(
      "__getitem__",
      [] (dmat44 const & A, std::pair<int, int> ij) {
        return A.data[ij.first][ij.second];
      }
    )
    .def(
      "__mul__",
      [] (dmat44 const & A, dvec4 const & x) { return dmat44_dvec4_mul(A, x); }
    )
    .def(
      "__mul__",
      [] (dmat44 const & A, dmat44 const & B) { return dmat44_dmat44_mul(A, B); }
    )
    ;

  m.def(
    "col",
    [] (dmat44 const & A, int j) { return dmat44_col(A, j); }
  );

  // eik_F3.h

  py::class_<F3_context>(m, "F3Context")
    .def(py::init(
           [] (cubic const & T_cubic, dvec2 const & xy, dvec2 const & xy0,
               dvec2 const & xy1, field2_wrapper const & slow) {
             auto ptr = std::make_unique<F3_context>();
             ptr->T_cubic = T_cubic;
             ptr->xy = xy;
             ptr->xy0 = xy0;
             ptr->xy1 = xy1;
             ptr->slow = &slow.field;
             ptr->F3 = NAN;
             ptr->F3_eta = NAN;
             return ptr;
           }
         ))
    .def_readwrite("T_cubic", &F3_context::T_cubic)
    .def_readwrite("xy", &F3_context::xy)
    .def_readwrite("xy0", &F3_context::xy0)
    .def_readwrite("xy1", &F3_context::xy1)
    .def(
      "compute",
      [] (F3_context & context, dbl eta) { F3_compute(eta, &context); }
    )
    .def_readonly("F3", &F3_context::F3)
    .def_readonly("F3_eta", &F3_context::F3_eta)
    ;

  // eik_F4.h

  py::class_<F4_context>(m, "F4Context")
    .def(py::init(
           [] (cubic const & T_cubic, cubic const & Tx_cubic,
               cubic const & Ty_cubic, dvec2 const & xy, dvec2 const & xy0,
               dvec2 const & xy1, field2_wrapper const & slow) {
             auto ptr = std::make_unique<F4_context>();
             ptr->T_cubic = T_cubic;
             ptr->Tx_cubic = Tx_cubic;
             ptr->Ty_cubic = Ty_cubic;
             ptr->xy = xy;
             ptr->xy0 = xy0;
             ptr->xy1 = xy1;
             ptr->slow = &slow.field;
             ptr->F4 = NAN;
             ptr->F4_eta = NAN;
             ptr->F4_th = NAN;
             return ptr;
           }
         ))
    .def_readwrite("T_cubic", &F4_context::T_cubic)
    .def_readwrite("Tx_cubic", &F4_context::Tx_cubic)
    .def_readwrite("Ty_cubic", &F4_context::Ty_cubic)
    .def_readwrite("xy", &F4_context::xy)
    .def_readwrite("xy0", &F4_context::xy0)
    .def_readwrite("xy1", &F4_context::xy1)
    .def(
      "compute",
      [] (F4_context & context, dbl eta, dbl th) {
        F4_compute(eta, th, &context);
      }
    )
    .def(
      "hess_fd",
      [] (F4_context & context, dbl eta, dbl th, dbl eps = 1e-7) {
        return F4_hess_fd(eta, th, eps, &context);
      }
    )
    .def(
      "bfgs_init",
      [] (F4_context & context, dbl eta, dbl th) {
        dvec2 x0, g0;
        dmat22 H0;
        F4_bfgs_init(eta, th, &x0, &g0, &H0, &context);
        return std::make_tuple(x0, g0, H0);
      }
    )
    .def(
      "bfgs_step",
      [] (F4_context & context, dvec2 xk, dvec2 gk, dmat22 Hk) {
        dvec2 xk1, gk1;
        dmat22 Hk1;
        return std::make_tuple(
          F4_bfgs_step(xk, gk, Hk, &xk1, &gk1, &Hk1, &context),
          xk1, gk1, Hk1);
      }
    )
    .def_readonly("F4", &F4_context::F4)
    .def_readonly("F4_eta", &F4_context::F4_eta)
    .def_readonly("F4_th", &F4_context::F4_th)
    ;

  // vec.h

  py::class_<dvec2>(m, "Dvec2")
    .def(py::init<dbl, dbl>())
    .def_readwrite("x", &dvec2::x)
    .def_readwrite("y", &dvec2::y)
    .def(
      "__repr__",
      [] (dvec2 const & u) {
        std::stringstream ss;
        ss << "Dvec2(" << u.x << ", " << u.y << ")";
        return ss.str();
      }
    )
    .def(
      "__getitem__",
      [] (dvec2 const & v, int i) {
        if (i < 0 || 2 <= i) {
          throw std::runtime_error {
            "index for Dvec2 must be in the interval [0, 2)"
          };
        }
        return i == 0 ? v.x : v.y;
      }
    )
    .def(
      "__len__",
      [] (dvec2 const & v) {
        (void) v;
        return 2;
      }
    )
    .def(
      "__mul__",
      [] (dvec2 const & u, dvec2 const & v) {
        return dvec2_dot(u, v);
      }
    )
    .def(
      "__add__",
      [] (dvec2 const & u, dvec2 const & v) {
        return dvec2_add(u, v);
      }
    )
    .def(
      "__sub__",
      [] (dvec2 const & u, dvec2 const & v) {
        return dvec2_sub(u, v);
      }
    )
    .def(
      "__truediv__",
      [] (dvec2 const & v, dbl a) {
        return dvec2_dbl_div(v, a);
      }
    )
    .def(
      "__mul__",
      [] (dvec2 const & v, dbl a) {
        return dvec2_dbl_mul(v, a);
      }
    )
    .def(
      "floor",
      [] (dvec2 const & v) {
        return dvec2_floor(v);
      }
    )
    .def(
      "norm",
      [] (dvec2 const & v) {
        return dvec2_norm(v);
      }
    )
    .def(
      "normalize",
      [] (dvec2 & v) {
        dvec2_normalize(&v);
      }
    )
    ;

  py::class_<dvec4>(m, "Dvec4")
    .def(py::init(
           [] (std::array<dbl, 4> const & data) {
             auto ptr = std::make_unique<dvec4>();
             for (int i = 0; i < 4; ++i) {
               ptr->data[i] = data[i];
             }
             return ptr;
           }
         ))
    .def(
      "__repr__",
      [] (dvec4 const & u) {
        std::stringstream ss;
        ss << "Dvec4("
           << u.data[0] << ", "
           << u.data[1] << ", "
           << u.data[2] << ", "
           << u.data[3] << ")";
        return ss.str();
      }
    )
    .def(
      "__getitem__",
      [] (dvec4 const & v, int i) {
        if (i < 0 || 4 <= i) {
          throw std::runtime_error {
            "index for Dvec4 must be in the interval [0, 4)"
          };
        }
        return v.data[i];
      }
    )
    .def(
      "__mul__",
      [] (dvec4 const & x, dmat44 const & A) { return dvec4_dmat44_mul(x, A); }
    )
    .def(
      "dot",
      [] (dvec4 const & u, dvec4 const & v) { return dvec4_dot(u, v); }
    )
    .def(
      "sum",
      [] (dvec4 const & u) { return dvec4_sum(u); }
    )
    .def(
      "__add__",
      [] (dvec4 const & u, dvec4 const & v) { return dvec4_add(u, v); }
    )
    .def(
      "__truediv__",
      [] (dvec4 const & u, dbl a) { return dvec4_dbl_div(u, a); }
    )
    .def_static(
      "m",
      [] (dbl x) { return dvec4_m(x); }
    )
    .def_static(
      "dm",
      [] (dbl x) { return dvec4_dm(x); }
    )
    .def_static(
      "e1",
      [] () { return dvec4_e1(); }
    )
    .def_static(
      "one",
      [] () { return dvec4_one(); }
    )
    .def_static(
      "iota",
      [] () { return dvec4_iota(); }
    )
    ;

  m.def(
    "dot",
    [] (dvec4 const & u, dvec4 const & v) { return dvec4_dot(u, v); }
  );

  m.def(
    "sum",
    [] (dvec4 const & u) { return dvec4_sum(u); }
  );

  py::class_<ivec2>(m, "Ivec2")
    .def(py::init<int, int>())
    .def(py::init(
           [] (std::pair<int, int> const & ind) {
             auto ptr = std::make_unique<ivec2>();
             ptr->i = ind.first;
             ptr->j = ind.second;
             return ptr;
           }
         ))
    .def(py::init(
           [] (dvec2 const & u) {
             return std::make_unique<ivec2>(dvec2_to_ivec2(u));
           }
         ))
    .def(
      "__repr__",
      [] (ivec2 const & u) {
        std::stringstream ss;
        ss << "Ivec2(" << u.i << ", " << u.j << ")";
        return ss.str();
      }
    )
    .def_readwrite("i", &ivec2::i)
    .def_readwrite("j", &ivec2::j)
    ;

  m.def(
    "ccomb",
    [] (dvec2 const & v0, dvec2 const & v1, dbl t) {
      return dvec2_ccomb(v0, v1, t);
    }
  );

  m.def(
    "cproj",
    [] (dvec2 const & u, dvec2 const & v) {
      return dvec2_cproj(u, v);
    }
  );

  m.def(
    "dist",
    [] (dvec2 const & v0, dvec2 const & v1) {
      return dvec2_dist(v0, v1);
    }
  );

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
