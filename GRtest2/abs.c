/*
 * abs.c
 *
 * Code generation for function 'abs'
 *
 * C source code generated on: Fri Apr 11 11:48:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen2.h"
#include "abs.h"
#include "Function_PCG_Wood_clean_codegen2_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  uint16_T uv0[2];
  int32_T i8;
  int32_T k;
  for (i8 = 0; i8 < 2; i8++) {
    uv0[i8] = (uint16_T)x->size[i8];
  }

  i8 = y->size[0] * y->size[1];
  y->size[0] = uv0[0];
  y->size[1] = uv0[1];
  emxEnsureCapacity((emxArray__common *)y, i8, (int32_T)sizeof(real_T));
  i8 = x->size[0] * x->size[1];
  for (k = 0; k <= i8 - 1; k++) {
    y->data[k] = fabs(x->data[k]);
  }
}

/* End of code generation (abs.c) */
