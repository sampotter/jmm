#include "edge.h"

#include "macros.h"

edge_s make_edge(size_t l0, size_t l1) {
  return (edge_s) {.l = {MIN(l0, l1), MAX(l0, l1)}};
}
