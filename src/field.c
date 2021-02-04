#include "field.h"

dbl field2_f(field2_s const *field, dvec2 xy) {
  return field->f(xy.x, xy.y, field->context);
}

dvec2 field2_grad_f(field2_s const *field, dvec2 xy) {
  return field->grad_f(xy.x, xy.y, field->context);
}
