#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "def.h"
#include "heap.h"
#include "hermite.h"
#include "jet.h"

struct heap_wrapper {
  heap * ptr {nullptr};
  heap_wrapper(int capacity,
               std::function<dbl(int)> const & value,
               std::function<void(int,int)> const & setpos) {
    heap_alloc(&ptr);
    heap_init(
      ptr,
      capacity,
      ^(int l) { return value(l); },
      ^(int l, int pos) { setpos(l, pos); }
    );
  }
  ~heap_wrapper() {
    heap_deinit(ptr);
    heap_dealloc(&ptr);
  }
};

PYBIND11_MODULE (sjs_eik, m) {
  m.doc() = R"pbdoc(
sjs_eik
-------

TODO!
)pbdoc";

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
    .def(
      "front",
      [] (heap_wrapper const & w) { return heap_front(w.ptr); }
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

  // hermite.h

  py::enum_<bicubic_variable>(m, "BicubicVariable")
    .value("Lambda", bicubic_variable::LAMBDA)
    .value("Mu", bicubic_variable::MU)
    ;

  py::class_<cubic>(m, "Cubic")
    .def(py::init(
           [] (std::array<dbl, 4> const & a) {
             auto ptr = std::make_unique<cubic>();
             std::copy(std::begin(a), std::end(a), ptr->a);
             return ptr;
           }
         ))
    .def(
      "f",
      [] (cubic const & C, dbl lam) { return cubic_f(&C, lam); }
    )
    .def(
      "df",
      [] (cubic const & C, dbl lam) { return cubic_df(&C, lam); }
    )
    ;

  py::class_<bicubic>(m, "Bicubic")
    .def(py::init(
           [] (std::array<std::array<dbl, 4>, 4> const & data) {
             auto ptr = std::make_unique<bicubic>();
             dbl data_arr[4][4];
             for (int i = 0; i < 4; ++i) {
               for (int j = 0; j < 4; ++j) {
                 data_arr[i][j] = data[i][j];
               }
             }
             bicubic_set_A(ptr.get(), data_arr);
             return ptr;
           }
         ))
    .def_property_readonly(
      "A",
      [] (bicubic const & B) {
        std::array<std::array<dbl, 4>, 4> A;
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            A[i][j] = B.A[i][j];
          }
        }
        return A;
      }
    )
    .def(
      "set_A",
      [] (bicubic & B, std::array<std::array<dbl, 4>, 4> const & data) {
        dbl data_arr[4][4];
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            data_arr[i][j] = data[i][j];
          }
        }
        bicubic_set_A(&B, data_arr);
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

  // jet.h

  py::class_<jet>(m, "Jet")
    .def(py::init<dbl, dbl, dbl, dbl>())
    .def_readwrite("f", &jet::f)
    .def_readwrite("fx", &jet::fx)
    .def_readwrite("fy", &jet::fy)
    .def_readwrite("fxy", &jet::fxy)
    ;

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}
