#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stddef.h>

typedef struct edge {
  size_t l[2];
} edge_s;

edge_s make_edge(size_t l0, size_t l1);

#ifdef __cplusplus
}
#endif
