#pragma once

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

#define SWAP(x, y) {                            \
    __typeof(x) tmp = x;                        \
    x = y;                                      \
    y = tmp;                                    \
  }

#define SORT2(x, y) {                           \
    __typeof(x) tmp = MIN(x, y);                \
    y = MAX(x, y);                              \
    x = tmp;                                    \
  }

#define SORT3(x, y, z)  {                                               \
    SORT2(x, y);                                                        \
    SORT2(x, z);                                                        \
    SORT2(y, z);                                                        \
  }
