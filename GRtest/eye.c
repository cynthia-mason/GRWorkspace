/*
 * eye.c
 *
 * Code generation for function 'eye'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen.h"
#include "eye.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static real_T rt_roundd_snf(real_T u);

/* Function Definitions */
static real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

void eye(real_T n, real_T I_data[289], int32_T I_size[2])
{
  int32_T q;
  int32_T i;
  real_T d0;
  I_size[0] = (int32_T)n;
  I_size[1] = (int32_T)n;
  q = (int32_T)n * (int32_T)n - 1;
  for (i = 0; i <= q; i++) {
    I_data[i] = 0.0;
  }

  d0 = rt_roundd_snf(n);
  if (d0 < 2.147483648E+9) {
    if (d0 >= -2.147483648E+9) {
      q = (int32_T)d0;
    } else {
      q = MIN_int32_T;
    }
  } else if (d0 >= 2.147483648E+9) {
    q = MAX_int32_T;
  } else {
    q = 0;
  }

  if (q > 0) {
    for (i = 0; i + 1 <= q; i++) {
      I_data[i + I_size[0] * i] = 1.0;
    }
  }
}

/* End of code generation (eye.c) */
