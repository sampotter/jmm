#pragma once

typedef enum jmm_error {
  JMM_ERROR_NONE,
  JMM_ERROR_BAD_ARGUMENTS,
  JMM_ERROR_RUNTIME_ERROR
} jmm_error_e;

char const *jmm_error_to_string(jmm_error_e error);
