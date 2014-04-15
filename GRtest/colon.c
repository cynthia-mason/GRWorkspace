/*
 * colon.c
 *
 * Code generation for function 'colon'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen.h"
#include "colon.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void eml_signed_integer_colon(int32_T b, int32_T y_data[17], int32_T y_size[2])
{
  int32_T n;
  int32_T yk;
  int32_T k;
  if (b < 1) {
    n = 0;
  } else {
    n = b;
  }

  y_size[0] = 1;
  y_size[1] = n;
  if (n > 0) {
    y_data[0] = 1;
    yk = 1;
    for (k = 2; k <= n; k++) {
      yk++;
      y_data[k - 1] = yk;
    }
  }
}

/* End of code generation (colon.c) */
