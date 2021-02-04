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
