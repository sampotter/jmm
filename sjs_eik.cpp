#include <pybind11/pybind11.h>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "bicubic.h"
#include "heap.h"
#include "index.h"
#include "jet.h"
#include "sjs_eik.h"

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

static dbl s_wrapper(void *vp, dvec2 xy);
static dvec2 grad_s_wrapper(void *vp, dvec2 xy);

struct sjs_wrapper {
  sjs * ptr {nullptr};

  std::optional<std::function<dbl(dbl, dbl)>> s;
  std::optional<std::function<std::tuple<dbl, dbl>(dbl, dbl)>> grad_s;

  sjs_wrapper(std::array<int, 2> const & shape,
              std::array<dbl, 2> const & xymin,
              dbl h,
              std::function<dbl(dbl, dbl)> s,
              std::function<std::tuple<dbl, dbl>(dbl, dbl)> grad_s):
    s {s},
    grad_s {grad_s}
  {
    sjs_alloc(&ptr);
    sjs_init(
      ptr,
      ivec2 {shape[0], shape[1]},
      dvec2 {xymin[0], xymin[1]},
      h,
      s_wrapper,
      grad_s_wrapper,
      (void *) this
    );
  }

  ~sjs_wrapper() {
    sjs_deinit(ptr);
    sjs_dealloc(&ptr);
  }
};

static dbl s_wrapper(void *vp, dvec2 xy) {
  sjs_wrapper * swp = (sjs_wrapper *) vp;
  if (!swp->s) {
    throw std::runtime_error {"ERROR: No s function for sjs!"};
  }
  return (*swp->s)(xy.x, xy.y);
}

static dvec2 grad_s_wrapper(void *vp, dvec2 xy) {
  sjs_wrapper * swp = (sjs_wrapper *) vp;
  if (!swp->grad_s) {
    throw std::runtime_error {"ERROR: No grad_s function for sjs!"};
  }
  dvec2 tmp;
  std::tie(tmp.x, tmp.y) = (*swp->grad_s)(xy.x, xy.y);
  return tmp;
}

PYBIND11_MODULE (_sjs_eik, m) {
  m.doc() = R"pbdoc(
_sjs_eik
--------

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
      "restrict",
      [] (bicubic const & B, bicubic_variable var, int edge) {
        return bicubic_restrict(&B, var, edge);
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
        std::optional<int> size;
        if (heap_size(w.ptr) > 0) {
          *size = heap_size(w.ptr);
        }
        return size;
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

  // sjs_eik.h

  py::class_<sjs_wrapper>(m, "StaticJetScheme")
    .def(py::init<
           std::array<int, 2> const &,
           std::array<dbl, 2> const &,
           dbl,
           std::function<dbl(dbl, dbl)> const &,
           std::function<std::tuple<dbl, dbl>(dbl, dbl)> const &
         >())
    .def(
      "step",
      [] (sjs_wrapper const & w) { sjs_step(w.ptr); }
    )
    .def(
      "solve",
      [] (sjs_wrapper const & w) { sjs_solve(w.ptr); }
    )
    .def(
      "add_trial",
      [] (sjs_wrapper const & w, int i, int j, jet_s jet) {
        sjs_add_trial(w.ptr, ivec2 {i, j}, jet);
      }
    )
    .def(
      "add_valid",
      [] (sjs_wrapper const & w, int i, int j, jet_s jet) {
        sjs_add_valid(w.ptr, ivec2 {i, j}, jet);
      }
    )
    .def(
      "make_bd",
      [] (sjs_wrapper const & w, int i, int j) {
        sjs_make_bd(w.ptr, ivec2 {i, j});
      }
    )
    .def_property_readonly(
      "shape",
      [] (sjs_wrapper const &w) { return sjs_get_shape(w.ptr); }
    )
    .def(
      "get_jet",
      [] (sjs_wrapper const & w, int i, int j) {
        return sjs_get_jet(w.ptr, {i, j});
      }
    )
    .def(
      "set_jet",
      [] (sjs_wrapper const & w, int i, int j, jet_s jet) {
        sjs_set_jet(w.ptr, {i, j}, jet);
      }
    )
    .def(
      "get_state",
      [] (sjs_wrapper const & w, int i, int j) {
        return sjs_get_state(w.ptr, {i, j});
      }
    )
    .def_property_readonly(
      "states",
      [] (sjs_wrapper const & w) {
        ivec2 shape = sjs_get_shape(w.ptr);
        return py::array {
          {shape.i, shape.j},
          {sizeof(state)*shape.j, sizeof(state)},
          sjs_get_states_ptr(w.ptr)
        };
      }
    )
    .def(
      "T",
      [] (sjs_wrapper const & w, dbl x, dbl y) {
        return sjs_T(w.ptr, dvec2 {x, y});
      }
    )
    .def(
      "Tx",
      [] (sjs_wrapper const & w, dbl x, dbl y) {
        return sjs_Tx(w.ptr, dvec2 {x, y});
      }
    )
    .def(
      "Ty",
      [] (sjs_wrapper const & w, dbl x, dbl y) {
        return sjs_Ty(w.ptr, dvec2 {x, y});
      }
    )
    .def(
      "Txy",
      [] (sjs_wrapper const & w, dbl x, dbl y) {
        return sjs_Txy(w.ptr, dvec2 {x, y});
      }
    )
    .def(
      "can_build_cell",
      [] (sjs_wrapper const & w, int i, int j) {
        return sjs_can_build_cell(w.ptr, ivec2 {i, j});
      }
    )
    .def(
      "build_cells",
      [] (sjs_wrapper const & w) { sjs_build_cells(w.ptr); }
    )
    .def(
      "get_bicubic",
      [] (sjs_wrapper const & w, int i, int j) {
        return sjs_get_bicubic(w.ptr, ivec2 {i, j});
      }
    )
    .def_property_readonly(
      "bicubics",
      [] (sjs_wrapper const & w) {
        ivec2 shape = sjs_get_shape(w.ptr);
        return py::array {
          {shape.i - 1, shape.j - 1},
          {sizeof(bicubic)*(shape.j - 1), sizeof(bicubic)},
          sjs_get_bicubics_ptr(w.ptr)
        };
      }
    )
    .def_property_readonly(
      "heap",
      [] (sjs_wrapper const & w) {
        return heap_wrapper {sjs_get_heap(w.ptr)};
      }
    )
    ;

  // vec.h

  py::class_<dvec2>(m, "Dvec2")
    .def(py::init<dbl, dbl>())
    .def_readwrite("x", &dvec2::x)
    .def_readwrite("y", &dvec2::y)
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
      "floor",
      [] (dvec2 const & v) {
        return dvec2_floor(v);
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
