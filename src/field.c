#include <jmm/field.h>

#include "macros.h"

dbl field2_f(field2_s const *field, dbl2 const x) {
  return field->f(SPLAT2(x), field->context);
}

void field2_grad_f(field2_s const *field, dbl2 const x, dbl2 df) {
  return field->grad_f(SPLAT2(x), field->context, df);
}
