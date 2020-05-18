#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void npy_write_2d_dbl_array(char const *filename, void *data, int m, int n, int stride);

#ifdef __cplusplus
}
#endif
