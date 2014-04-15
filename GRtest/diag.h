/*
 * diag.h
 *
 * Code generation for function 'diag'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

#ifndef __DIAG_H__
#define __DIAG_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "Function_PCG_Wood_clean_codegen_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void b_diag(const emxArray_real_T *v, emxArray_real_T *d);
extern void diag(const real_T v_data[17], const int32_T v_size[1], real_T d_data[289], int32_T d_size[2]);
#endif
/* End of code generation (diag.h) */
