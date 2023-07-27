#pragma once

#include "common.h"
#include "vec.h"

struct field2 {
  dbl(*f)(dbl, dbl, void*);
  void(*grad_f)(dbl, dbl, void*, dbl2);
  void *context;
};

dbl field2_f(field2_s const *field, dbl2 const x);
void field2_grad_f(field2_s const *field, dbl2 const x, dbl2 Df);

struct field3 {
  dbl(*f)(dbl3, void*);
  void(*grad_f)(dbl3, void*, dbl3);
  void *context;
};

dbl field3_f(field3_s const *field, dbl3 const x);
void field3_grad_f(field3_s const *field, dbl3 const x, dbl3 Df);
