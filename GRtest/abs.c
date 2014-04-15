/*
 * abs.c
 *
 * Code generation for function 'abs'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen.h"
#include "abs.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_abs(const real_T x_data[289], const int32_T x_size[2], real_T y_data[289],
           int32_T y_size[2])
{
  int8_T iv0[2];
  int32_T i9;
  int32_T k;
  for (i9 = 0; i9 < 2; i9++) {
    iv0[i9] = (int8_T)x_size[i9];
  }

  y_size[0] = iv0[0];
  y_size[1] = iv0[1];
  i9 = x_size[0] * x_size[1];
  for (k = 0; k <= i9 - 1; k++) {
    y_data[k] = fabs(x_data[k]);
  }
}

/* End of code generation (abs.c) */
