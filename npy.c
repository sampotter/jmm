#include "npy.h"

#include "def.h"

#include <stdio.h>

void npy_write_2d_dbl_array(char const *filename, void *data, int m, int n, int stride) {
  FILE *stream = fopen(filename, "w");

  /**
   * Write header
   *
   * See:
   *
   * https://docs.scipy.org/doc/numpy/reference/generated/numpy.lib.format.html
   *   #module-numpy.lib.format
   *
   * for an explanation of how this should work. We just use version
   * 1.0 of the npy file format for simplicity.
   */

  // Write magic
  fprintf(stream, "\x93NUMPY");

  // Write major version
  fputc('\x01', stream);
  fputc('\x00', stream);

  // Write header dictionary
  char buffer[1 << 16];
  unsigned short nbytes = sprintf(
    buffer,
    "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d), }",
    m, n
    );

  int rem = 64 - ((10 + nbytes) % 64);
  int headerlen = nbytes + rem;

  // Write length of header
  fwrite((void *)&headerlen, sizeof(unsigned short), 1, stream);

  // Write header
  fprintf(stream, "%s", buffer);

  // Pad the rest of the header with spaces for alignment
  for (int i = 0; i < rem - 1; ++i) {
    fputc('\x20', stream);
  }
  fputc('\n', stream);

  // Write the data
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      fwrite(&data[stride*(n*i + j)], sizeof(dbl), 1, stream);
    }
  }

  fclose(stream);
}
