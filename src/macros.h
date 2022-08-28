#pragma once

#define die() assert(false)

#define MIN(x, y) ({                            \
    __typeof(x) x_ = (x);                       \
    __typeof(y) y_ = (y);                       \
    x_ < y_ ? x_ : y_;                          \
  })

#define MAX(x, y) ({                            \
    __typeof(x) x_ = (x);                       \
    __typeof(y) y_ = (y);                       \
    x_ > y_ ? x_ : y_;                          \
  })

#define SWAP(x, y) do {                         \
    __typeof(x) tmp = x;                        \
    x = y;                                      \
    y = tmp;                                    \
  } while (0)

#define ROTATE3(x, y, z) do {                                           \
    __typeof(x) tmp = x;                                                \
    x = y;                                                              \
    y = z;                                                              \
    z = tmp;                                                            \
  } while (0)

#define SORT2(x, y) do {                        \
    __typeof(x) tmp = MIN(x, y);                \
    y = MAX(x, y);                              \
    x = tmp;                                    \
  } while (0)

#define SORT3(x, y, z) do {                                             \
    SORT2(x, y);                                                        \
    SORT2(x, z);                                                        \
    SORT2(y, z);                                                        \
  } while (0)

#define SORT_UINT2(l) SORT2(l[0], l[1])
#define SORT_UINT3(l) SORT3(l[0], l[1], l[2])

#define SPLAT2(arr) arr[0], arr[1]
