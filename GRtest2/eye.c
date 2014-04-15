/*
 * eye.c
 *
 * Code generation for function 'eye'
 *
 * C source code generated on: Fri Apr 11 11:48:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen2.h"
#include "eye.h"
#include "Function_PCG_Wood_clean_codegen2_emxutil.h"

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

void eye(real_T n, emxArray_real_T *I)
{
  int32_T q;
  int32_T loop_ub;
  real_T d0;
  q = I->size[0] * I->size[1];
  I->size[0] = (int32_T)n;
  I->size[1] = (int32_T)n;
  emxEnsureCapacity((emxArray__common *)I, q, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)n * (int32_T)n - 1;
  for (q = 0; q <= loop_ub; q++) {
    I->data[q] = 0.0;
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
    for (loop_ub = 0; loop_ub + 1 <= q; loop_ub++) {
      I->data[loop_ub + I->size[0] * loop_ub] = 1.0;
    }
  }
}

/* End of code generation (eye.c) */
