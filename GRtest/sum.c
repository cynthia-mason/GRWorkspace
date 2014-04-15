/*
 * sum.c
 *
 * Code generation for function 'sum'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen.h"
#include "sum.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
real_T sum(const real_T x[10])
{
  real_T y;
  int32_T k;
  y = x[0];
  for (k = 0; k < 9; k++) {
    y += x[k + 1];
  }

  return y;
}

/* End of code generation (sum.c) */
