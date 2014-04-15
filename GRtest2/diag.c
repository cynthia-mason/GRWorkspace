/*
 * diag.c
 *
 * Code generation for function 'diag'
 *
 * C source code generated on: Fri Apr 11 11:48:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen2.h"
#include "diag.h"
#include "Function_PCG_Wood_clean_codegen2_emxutil.h"

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
  emxArray_int32_T *r7;
  emxArray_int32_T *r8;
  emxArray_real_T *b_v;
  emxArray_real_T *c_v;
  emxArray_real_T *d_v;
  int32_T stride;
  m = v->size[0];
  n = v->size[1];
  b_emxInit_int32_T(&r7, 2);
  emxInit_int32_T(&r8, 1);
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
    n = r8->size[0];
    r8->size[0] = m / stride + 1;
    emxEnsureCapacity((emxArray__common *)r8, n, (int32_T)sizeof(int32_T));
    m /= stride;
    for (n = 0; n <= m; n++) {
      r8->data[n] = 1 + stride * n;
    }

    n = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)r7, n, (int32_T)sizeof(int32_T));
    m = r8->size[0];
    n = r7->size[0] * r7->size[1];
    r7->size[1] = m;
    emxEnsureCapacity((emxArray__common *)r7, n, (int32_T)sizeof(int32_T));
    m = r8->size[0] - 1;
    for (n = 0; n <= m; n++) {
      r7->data[n] = r8->data[n];
    }

    m = r7->size[1];
    n = b_v->size[0];
    b_v->size[0] = m;
    emxEnsureCapacity((emxArray__common *)b_v, n, (int32_T)sizeof(real_T));
    m--;
    for (n = 0; n <= m; n++) {
      b_v->data[n] = v->data[r7->data[n] - 1];
    }

    m = b_v->size[0];
    n = d->size[0] * d->size[1];
    d->size[0] = m;
    emxEnsureCapacity((emxArray__common *)d, n, (int32_T)sizeof(real_T));
    m = r7->size[1];
    n = c_v->size[0];
    c_v->size[0] = m;
    emxEnsureCapacity((emxArray__common *)c_v, n, (int32_T)sizeof(real_T));
    m--;
    for (n = 0; n <= m; n++) {
      c_v->data[n] = v->data[r7->data[n] - 1];
    }

    n = d->size[0] * d->size[1];
    d->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)d, n, (int32_T)sizeof(real_T));
    m = r7->size[1];
    n = d_v->size[0];
    d_v->size[0] = m;
    emxEnsureCapacity((emxArray__common *)d_v, n, (int32_T)sizeof(real_T));
    m--;
    for (n = 0; n <= m; n++) {
      d_v->data[n] = v->data[r7->data[n] - 1];
    }

    m = d_v->size[0] - 1;
    for (n = 0; n <= m; n++) {
      d->data[n] = v->data[r7->data[n] - 1];
    }
  }

  emxFree_real_T(&d_v);
  emxFree_real_T(&c_v);
  emxFree_real_T(&b_v);
  emxFree_int32_T(&r8);
  emxFree_int32_T(&r7);
}

void diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int32_T nv;
  int32_T j;
  int32_T loop_ub;
  nv = v->size[0];
  j = d->size[0] * d->size[1];
  d->size[0] = (uint16_T)nv;
  emxEnsureCapacity((emxArray__common *)d, j, (int32_T)sizeof(real_T));
  j = d->size[0] * d->size[1];
  d->size[1] = (uint16_T)nv;
  emxEnsureCapacity((emxArray__common *)d, j, (int32_T)sizeof(real_T));
  loop_ub = (uint16_T)nv * (uint16_T)nv - 1;
  for (j = 0; j <= loop_ub; j++) {
    d->data[j] = 0.0;
  }

  for (j = 0; j + 1 <= nv; j++) {
    d->data[j + d->size[0] * j] = v->data[j];
  }
}

/* End of code generation (diag.c) */
