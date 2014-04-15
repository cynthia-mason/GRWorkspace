/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Fri Apr 11 11:48:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen2.h"
#include "rdivide.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void rdivide(const real_T x[17], const real_T y[17], real_T z[17])
{
  int32_T i;
  for (i = 0; i < 17; i++) {
    z[i] = x[i] / y[i];
  }
}

/* End of code generation (rdivide.c) */
