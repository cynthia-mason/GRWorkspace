/*
 * Function_PCG_Wood_clean_codegen_emxutil.h
 *
 * Code generation for function 'Function_PCG_Wood_clean_codegen_emxutil'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

#ifndef __FUNCTION_PCG_WOOD_CLEAN_CODEGEN_EMXUTIL_H__
#define __FUNCTION_PCG_WOOD_CLEAN_CODEGEN_EMXUTIL_H__
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
extern void b_emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
extern void b_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
extern void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
#endif
/* End of code generation (Function_PCG_Wood_clean_codegen_emxutil.h) */
