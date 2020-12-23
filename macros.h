#pragma once

#define SWAP(x, y) {                            \
    __typeof(x) tmp = x;                        \
    x = y;                                      \
    y = tmp;                                    \
  }
