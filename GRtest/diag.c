/*
 * diag.c
 *
 * Code generation for function 'diag'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen.h"
#include "diag.h"
#include "Function_PCG_Wood_clean_codegen_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int32_T m;
  int32_T n;
  emxArray_int32_T *r6;
  emxArray_int32_T *r7;
  emxArray_real_T *b_v;
  emxArray_real_T *c_v;
  emxArray_real_T *d_v;
  int32_T stride;
  m = v->size[0];
  n = v->size[1];
  emxInit_int32_T(&r6, 2);
  b_emxInit_int32_T(&r7, 1);
  b_emxInit_real_T(&b_v, 1);
  b_emxInit_real_T(&c_v, 1);
  b_emxInit_real_T(&d_v, 1);
  if ((m == 0) || (n == 0)) {
    n = d->size[0] * d->size[1];
    d->size[0] = 0;
    d->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)d, n, (int32_T)sizeof(real_T));
  } else {
    stride = m + 1;
    if (m <= n) {
      n = m;
    }

    m = stride * (n - 1);
    n = r7->size[0];
    r7->size[0] = m / stride + 1;
    emxEnsureCapacity((emxArray__common *)r7, n, (int32_T)sizeof(int32_T));
    m /= stride;
    for (n = 0; n <= m; n++) {
      r7->data[n] = 1 + stride * n;
    }

    n = r6->size[0] * r6->size[1];
    r6->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)r6, n, (int32_T)sizeof(int32_T));
    m = r7->size[0];
    n = r6->size[0] * r6->size[1];
    r6->size[1] = m;
    emxEnsureCapacity((emxArray__common *)r6, n, (int32_T)sizeof(int32_T));
    m = r7->size[0] - 1;
    for (n = 0; n <= m; n++) {
      r6->data[n] = r7->data[n];
    }

    m = r6->size[1];
    n = b_v->size[0];
    b_v->size[0] = m;
    emxEnsureCapacity((emxArray__common *)b_v, n, (int32_T)sizeof(real_T));
    m--;
    for (n = 0; n <= m; n++) {
      b_v->data[n] = v->data[r6->data[n] - 1];
    }

    m = b_v->size[0];
    n = d->size[0] * d->size[1];
    d->size[0] = m;
    emxEnsureCapacity((emxArray__common *)d, n, (int32_T)sizeof(real_T));
    m = r6->size[1];
    n = c_v->size[0];
    c_v->size[0] = m;
    emxEnsureCapacity((emxArray__common *)c_v, n, (int32_T)sizeof(real_T));
    m--;
    for (n = 0; n <= m; n++) {
      c_v->data[n] = v->data[r6->data[n] - 1];
    }

    n = d->size[0] * d->size[1];
    d->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)d, n, (int32_T)sizeof(real_T));
    m = r6->size[1];
    n = d_v->size[0];
    d_v->size[0] = m;
    emxEnsureCapacity((emxArray__common *)d_v, n, (int32_T)sizeof(real_T));
    m--;
    for (n = 0; n <= m; n++) {
      d_v->data[n] = v->data[r6->data[n] - 1];
    }

    m = d_v->size[0] - 1;
    for (n = 0; n <= m; n++) {
      d->data[n] = v->data[r6->data[n] - 1];
    }
  }

  emxFree_real_T(&d_v);
  emxFree_real_T(&c_v);
  emxFree_real_T(&b_v);
  emxFree_int32_T(&r7);
  emxFree_int32_T(&r6);
}

void diag(const real_T v_data[17], const int32_T v_size[1], real_T d_data[289],
          int32_T d_size[2])
{
  int32_T loop_ub;
  int32_T i8;
  d_size[0] = (int8_T)v_size[0];
  d_size[1] = (int8_T)v_size[0];
  loop_ub = (int8_T)v_size[0] * (int8_T)v_size[0] - 1;
  for (i8 = 0; i8 <= loop_ub; i8++) {
    d_data[i8] = 0.0;
  }

  for (loop_ub = 0; loop_ub + 1 <= v_size[0]; loop_ub++) {
    d_data[loop_ub + d_size[0] * loop_ub] = v_data[loop_ub];
  }
}

/* End of code generation (diag.c) */
