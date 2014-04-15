/*
 * colon.c
 *
 * Code generation for function 'colon'
 *
 * C source code generated on: Fri Apr 11 11:48:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen2.h"
#include "colon.h"
#include "Function_PCG_Wood_clean_codegen2_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void eml_signed_integer_colon(int32_T b, emxArray_int32_T *y)
{
  int32_T n;
  int32_T yk;
  int32_T k;
  if (b < 1) {
    n = 0;
  } else {
    n = b;
  }

  yk = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity((emxArray__common *)y, yk, (int32_T)sizeof(int32_T));
  if (n > 0) {
    y->data[0] = 1;
    yk = 1;
    for (k = 2; k <= n; k++) {
      yk++;
      y->data[k - 1] = yk;
    }
  }
}

/* End of code generation (colon.c) */
