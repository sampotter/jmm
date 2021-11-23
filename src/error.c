#include "error.h"

char const *jmm_error_to_string(jmm_error_e error) {
  static char const *error_strings[] = {
    "JMM_ERROR_NONE",
    "JMM_ERROR_BAD_ARGUMENTS",
    "JMM_ERROR_RUNTIME_ERROR"
  };
  return error_strings[error];
}
