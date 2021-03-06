/*
 * Function_PCG_Wood_clean_codegen2.c
 *
 * Code generation for function 'Function_PCG_Wood_clean_codegen2'
 *
 * C source code generated on: Fri Apr 11 11:48:45 2014
 *
 */

#include "stdio.h"
/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen2.h"
#include "Function_PCG_Wood_clean_codegen2_emxutil.h"
#include "cond.h"
#include "rdivide.h"
#include "mrdivide.h"
#include "eye.h"
#include "diag.h"
#include "abs.h"
#include "sum.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_eml_li_find(const boolean_T x_data[256], const int32_T x_size[1],
  int32_T y_data[256], int32_T y_size[1]);
static void b_eml_xgemm(int32_T m, int32_T n, const emxArray_real_T *A, int32_T
  lda, const real_T B_data[1169], const int32_T B_size[2], emxArray_real_T *C,
  int32_T ldc);
static void c_eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y);
static void c_eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
  lda, const emxArray_real_T *B, int32_T ldb, emxArray_real_T *C, int32_T ldc);
static int32_T compute_nones(const emxArray_boolean_T *x, int32_T n);
static void d_eml_xgemm(int32_T m, int32_T k, const real_T A_data[1169], const
  int32_T A_size[2], int32_T lda, const real_T B[10], int32_T ldb, real_T
  C_data[1169], int32_T C_size[1], int32_T ldc);
static void e_eml_xgemm(int32_T m, int32_T k, const real_T A_data[1169], const
  int32_T A_size[2], int32_T lda, const real_T B[100], int32_T ldb,
  emxArray_real_T *C, int32_T ldc);
static void eml_li_find(const boolean_T x[17], int32_T y_data[17], int32_T
  y_size[1]);
static void eml_xgemm(int32_T n, int32_T k, const real_T A_data[2730], const
                      int32_T A_size[2], const emxArray_real_T *B, int32_T ldb,
                      emxArray_real_T *C);

/* Function Definitions */
static void b_eml_li_find(const boolean_T x_data[256], const int32_T x_size[1],
  int32_T y_data[256], int32_T y_size[1])
{
  int32_T k;
  int32_T i;
  k = 0;
  for (i = 1; i <= x_size[0]; i++) {
    if (x_data[i - 1]) {
      k++;
    }
  }

  y_size[0] = k;
  k = 0;
  for (i = 1; i <= x_size[0]; i++) {
    if (x_data[i - 1]) {
      y_data[k] = i;
      k++;
    }
  }
}

static void b_eml_xgemm(int32_T m, int32_T n, const emxArray_real_T *A, int32_T
  lda, const real_T B_data[1169], const int32_T B_size[2], emxArray_real_T *C,
  int32_T ldc)
{
  emxArray_real_T *b_C;
  int32_T ar;
  int32_T ic;
  int32_T br;
  int32_T ib;
  int32_T ia;
  emxInit_real_T(&b_C, 2);
  if ((m == 0) || (n == 0)) {
  } else {
    ar = 0;
    while ((ldc > 0) && (ar <= 0)) {
      for (ic = 1; ic <= m; ic++) {
        C->data[ic - 1] = 0.0;
      }

      ar = ldc;
    }

    br = 0;
    ar = 0;
    while ((ldc > 0) && (ar <= 0)) {
      ar = -1;
      for (ib = br; ib + 1 <= br + 10; ib++) {
        if (B_data[ib] != 0.0) {
          ia = ar;
          for (ic = 0; ic + 1 <= m; ic++) {
            ia++;
            C->data[ic] += B_data[ib] * A->data[ia];
          }
        }

        ar += lda;
      }

      br += 10;
      ar = ldc;
    }
  }

  emxFree_real_T(&b_C);
}

static void c_eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y)
{
  int32_T n;
  int32_T k;
  int32_T i;
  n = x->size[0] * x->size[1];
  k = compute_nones(x, n);
  i = y->size[0];
  y->size[0] = k;
  emxEnsureCapacity((emxArray__common *)y, i, (int32_T)sizeof(int32_T));
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      y->data[k] = i;
      k++;
    }
  }
}

static void c_eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
  lda, const emxArray_real_T *B, int32_T ldb, emxArray_real_T *C, int32_T ldc)
{
  emxArray_real_T *b_C;
  int32_T c;
  int32_T cr;
  int32_T i15;
  int32_T ic;
  int32_T br;
  int32_T ar;
  int32_T ib;
  int32_T ia;
  int32_T i16;
  emxInit_real_T(&b_C, 2);
  if (m == 0) {
  } else {
    c = ldc * 9;
    cr = 0;
    while ((ldc > 0) && (cr <= c)) {
      i15 = cr + m;
      for (ic = cr; ic + 1 <= i15; ic++) {
        C->data[ic] = 0.0;
      }

      cr += ldc;
    }

    br = 0;
    cr = 0;
    while ((ldc > 0) && (cr <= c)) {
      ar = -1;
      i15 = br + k;
      for (ib = br; ib + 1 <= i15; ib++) {
        if (B->data[ib] != 0.0) {
          ia = ar;
          i16 = cr + m;
          for (ic = cr; ic + 1 <= i16; ic++) {
            ia++;
            C->data[ic] += B->data[ib] * A->data[ia];
          }
        }

        ar += lda;
      }

      br += ldb;
      cr += ldc;
    }
  }

  emxFree_real_T(&b_C);
}

static int32_T compute_nones(const emxArray_boolean_T *x, int32_T n)
{
  int32_T k;
  int32_T i;
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      k++;
    }
  }

  return k;
}

static void d_eml_xgemm(int32_T m, int32_T k, const real_T A_data[1169], const
  int32_T A_size[2], int32_T lda, const real_T B[10], int32_T ldb, real_T
  C_data[1169], int32_T C_size[1], int32_T ldc)
{
  int32_T cr;
  if (m == 0) {
  } else {
    cr = 0;
    while ((ldc > 0) && (cr <= 0)) {
      for (cr = 1; cr <= m; cr++) {
        C_data[cr - 1] = 0.0;
      }

      cr = ldc;
    }
  }
}

static void e_eml_xgemm(int32_T m, int32_T k, const real_T A_data[1169], const
  int32_T A_size[2], int32_T lda, const real_T B[100], int32_T ldb,
  emxArray_real_T *C, int32_T ldc)
{
  int32_T c;
  int32_T cr;
  int32_T i17;
  int32_T ic;
  if (m == 0) {
  } else {
    c = ldc * 9;
    cr = 0;
    while ((ldc > 0) && (cr <= c)) {
      i17 = cr + m;
      for (ic = cr; ic + 1 <= i17; ic++) {
        C->data[ic] = 0.0;
      }

      cr += ldc;
    }
  }
}

static void eml_li_find(const boolean_T x[17], int32_T y_data[17], int32_T
  y_size[1])
{
  int32_T k;
  int32_T i;
  k = 0;
  for (i = 0; i < 17; i++) {
    if (x[i]) {
      k++;
    }
  }

  y_size[0] = k;
  k = 0;
  for (i = 0; i < 17; i++) {
    if (x[i]) {
      y_data[k] = i + 1;
      k++;
    }
  }
}

static void eml_xgemm(int32_T n, int32_T k, const real_T A_data[2730], const
                      int32_T A_size[2], const emxArray_real_T *B, int32_T ldb,
                      emxArray_real_T *C)
{
  emxArray_real_T *b_C;
  int32_T c;
  int32_T cr;
  int32_T ic;
  int32_T br;
  int32_T ar;
  int32_T i12;
  int32_T ib;
  int32_T ia;
  emxInit_real_T(&b_C, 2);
  if (n == 0) {
  } else {
    c = 10 * (n - 1);
    for (cr = 0; cr <= c; cr += 10) {
      for (ic = cr; ic + 1 <= cr + 10; ic++) {
        C->data[ic] = 0.0;
      }
    }

    br = 0;
    for (cr = 0; cr <= c; cr += 10) {
      ar = -1;
      i12 = br + k;
      for (ib = br; ib + 1 <= i12; ib++) {
        if (B->data[ib] != 0.0) {
          ia = ar;
          for (ic = cr; ic + 1 <= cr + 10; ic++) {
            ia++;
            C->data[ic] += B->data[ib] * A_data[ia];
          }
        }

        ar += 10;
      }

      br += ldb;
    }
  }

  emxFree_real_T(&b_C);
}

void Function_PCG_Wood_clean_codegen2(const real_T A[170], const real_T b[17],
  const real_T c[10], real_T b_gamma, real_T x[10], real_T *NumVar, real_T
  *IPMit, real_T *CondNum, real_T *PCGit)
{
  int32_T Binv_size[2];
  real_T Binv_data[1169];
  real_T dx_1[10];
  int32_T i0;
  real_T epsilon;
  int32_T i;
  boolean_T b_A[17];
  int32_T tmp_size[1];
  int32_T tmp_data[17];
  int32_T loop_ub;
  real_T row_v_data[17];
  real_T row_nonneg1_data[17];
  real_T row_nonneg2_data[17];
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T check;
  int32_T b_check;
  emxArray_real_T *D_track_tempVec;
  emxArray_real_T *Dtreetrack;
  emxArray_real_T *u;
  emxArray_real_T *A_11inv;
  emxArray_real_T *A_12;
  emxArray_real_T *A_12_Psi;
  emxArray_boolean_T *A_12_Psi_to_remove;
  emxArray_real_T *Psi_notinv;
  emxArray_real_T *Psi;
  emxArray_real_T *Arrow_inv;
  emxArray_real_T *C;
  emxArray_int32_T *r0;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  emxArray_real_T *a;
  emxArray_real_T *c_y;
  emxArray_real_T *d_y;
  emxArray_real_T *maxval;
  emxArray_real_T *e_y;
  emxArray_real_T *r1;
  emxArray_real_T *r2;
  emxArray_int32_T *r3;
  emxArray_real_T *r4;
  emxArray_real_T *b_Dtreetrack;
  emxArray_real_T *r5;
  emxArray_real_T *c_Dtreetrack;
  emxArray_real_T *r6;
  emxArray_real_T *D;
  emxArray_real_T *D_track;
  emxArray_real_T *b_D;
  emxArray_real_T *c_D;
  emxArray_real_T *d_D;
  boolean_T exitg1;
  real_T w[17];
  int32_T i4;
  int32_T i5;
  real_T e_D[17];
  int32_T D_tree_size[2];
  real_T D_tree_data[289];
  int32_T i6;
  int32_T D_track_size[2];
  real_T D_track_data[256];
  int32_T D_nonneg_size[2];
  real_T D_nonneg_data[289];
  int32_T Atreetrack_size[2];
  int32_T i7;
  real_T Atreetrack_data[2730];
  int32_T ixstop;
  uint32_T outsz[2];
  emxArray_real_T b_D_tree_data;
  int32_T k;
  emxArray_real_T c_D_tree_data;
  real_T A_data[170];
  real_T addaTree[100];
  int32_T ic;
  int32_T br;
  int32_T ar;
  int32_T ib;
  int32_T ia;
  emxArray_real_T b_D_track_data;
  emxArray_real_T c_D_track_data;
  int32_T pos;
  real_T addaTrack[100];
  emxArray_real_T b_D_nonneg_data;
  emxArray_real_T c_D_nonneg_data;
  real_T addaNonneg[100];
  real_T alpha_2;
  boolean_T exitg8;
  real_T diag_ind_data[256];
  int16_T sz[2];
  real_T tracks_i_data[256];
  boolean_T b_tracks_i_data[256];
  int32_T tracks_i_size[1];
  int32_T b_tmp_data[256];
  real_T b_A_data[160];
  int32_T m;
  real_T diagB_data[145];
  boolean_T exitg7;
  boolean_T exitg6;
  boolean_T exitg5;
  boolean_T exitg4;
  real_T f_y;
  real_T r[10];
  int32_T z_size[1];
  real_T z_data[1169];
  real_T p_data[1169];
  boolean_T exitg3;
  real_T a_data[1169];
  real_T g_y[10];
  real_T dx_1_new[10];
  real_T r_new[10];
  int32_T z_new_size[1];
  real_T z_new_data[1169];
  real_T c_A[17];
  boolean_T exitg2;
  real_T h_y;
  real_T b_r;
  boolean_T guard1 = FALSE;
  emxArray_real_T *i_y;

  /* % This selects the best trees to achieve the minimla global routing. */
  /* Input: */
  /*    A = [A_tree;A_edge;A_nonneg] */
  /*    b = [b_tree;b_edge;b_nonneg] */
  /*    c = Objective function */
  /*    gamma = This is a value between 0 and 1 which dictates how large of a */
  /*        step to take when moving throughout the search space in the  */
  /*        Interior Point Method. Float (1,1) */
  /*     */
  /* Output: */
  /*    x = Selection of trees */
  /*    NumVar = Number of variables (int) */
  /*    IPMit = Number of IPM iterations (int) */
  /*    CondNum = Final condition number of MinvADDA */
  /*    PCGit = Number of PAG iterations (int) */
  /*     */
  /* % */
  /* ---------------------------------------------------- */
  /*  Preconditioned conjugate gradiant method */
  /*  THE INV FUNCTION IS NOT USED IN THIS CODE. Diagonal matrices are taken  */
  /*  advantage of. */
  /*  */
  /*  Must load the following: */
  /*  A=[1,1,1,1,0,0,0,0,0,0;0,0,0,0,1,1,1,0,0,0;0,0,0,0,0,0,0,1,1,1;1,0,0,0,1,0,0,1,0,0;0,1,0,0,0,1,0,0,0,0;0,0,0,0,0,0,1,0,0,1;0,0,1,0,0,0,0,0,1,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,-1,0,0,0,0,0,0,0;0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;0,0,0,0,0,0,-1,0,0,0;0,0,0,0,0,0,0,-1,0,0;0,0,0,0,0,0,0,0,-1,0;0,0,0,0,0,0,0,0,0,-1;]; */
  /*  b= [1;1;1;2;2;1;1;0;0;0;0;0;0;0;0;0;0;]; */
  /*  c=[10;8;6;4;5;4;13;16;7;8;]; */
  /*  gamma=0.9; */
  /* ---------------------------------------------------- */
  /* % Initialization  */
  /* A=sparse(rowA,colA); */
  Binv_size[0] = 1;
  Binv_size[1] = 1;
  Binv_data[0] = 0.0;
  for (i0 = 0; i0 < 10; i0++) {
    dx_1[i0] = A[17 * i0];
  }

  epsilon = 1.0 / (2.0 * sum(dx_1));
  *NumVar = 10.0;
  for (i = 0; i < 10; i++) {
    x[i] = epsilon;
    dx_1[i] = epsilon;
  }

  /* Find where the components of A begin and end [A_tree;A_edge;A_nonneg] */
  /* 1) The Tree constraints end */
  for (i0 = 0; i0 < 17; i0++) {
    b_A[i0] = (A[153 + i0] == 1.0);
  }

  eml_li_find(b_A, tmp_data, tmp_size);
  loop_ub = tmp_size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    row_v_data[i0] = 1.0 + (real_T)(tmp_data[i0] - 1);
  }

  /* Find the last row of net constratints (A must be constructed such that there are blocks of cascading ones) */
  /* row_v=find(A(:,colA)==1,1);  */
  /* 2) Where the nonnegativity constraints begin */
  for (i0 = 0; i0 < 17; i0++) {
    b_A[i0] = (A[i0] == -1.0);
  }

  eml_li_find(b_A, tmp_data, tmp_size);
  loop_ub = tmp_size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    row_nonneg1_data[i0] = 1.0 + (real_T)(tmp_data[i0] - 1);
  }

  /* row_nonneg12=find(A(:,1)==-1,1); */
  for (i0 = 0; i0 < 17; i0++) {
    b_A[i0] = (A[153 + i0] == -1.0);
  }

  eml_li_find(b_A, tmp_data, tmp_size);
  loop_ub = tmp_size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    row_nonneg2_data[i0] = 1.0 + (real_T)(tmp_data[i0] - 1);
  }

  /* row_nonneg22=find(A(:,colA)==-1,1); */
  if (row_nonneg1_data[0] > row_nonneg2_data[0]) {
    i0 = 0;
    i1 = 0;
  } else {
    i0 = (int32_T)row_nonneg1_data[0] - 1;
    i1 = (int32_T)row_nonneg2_data[0];
  }

  /* isscalar(track1) */
  /* isscalar(track2) */
  if (row_v_data[0] + 1.0 > row_nonneg1_data[0] - 1.0) {
    i2 = 0;
    i3 = -1;
  } else {
    i2 = (int32_T)row_v_data[0];
    i3 = (int32_T)row_nonneg1_data[0] - 2;
  }

  /* % Calculations */
  check = 1;
  b_check = 0;
  emxInit_real_T(&D_track_tempVec, 2);
  emxInit_real_T(&Dtreetrack, 2);
  emxInit_real_T(&u, 2);
  emxInit_real_T(&A_11inv, 2);
  emxInit_real_T(&A_12, 2);
  emxInit_real_T(&A_12_Psi, 2);
  emxInit_boolean_T(&A_12_Psi_to_remove, 2);
  emxInit_real_T(&Psi_notinv, 2);
  emxInit_real_T(&Psi, 2);
  emxInit_real_T(&Arrow_inv, 2);
  emxInit_real_T(&C, 2);
  emxInit_int32_T(&r0, 1);
  emxInit_real_T(&y, 2);
  emxInit_real_T(&b_y, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&c_y, 2);
  emxInit_real_T(&d_y, 2);
  emxInit_real_T(&maxval, 2);
  emxInit_real_T(&e_y, 2);
  emxInit_real_T(&r1, 2);
  emxInit_real_T(&r2, 2);
  emxInit_int32_T(&r3, 1);
  emxInit_real_T(&r4, 2);
  emxInit_real_T(&b_Dtreetrack, 2);
  emxInit_real_T(&r5, 2);
  emxInit_real_T(&c_Dtreetrack, 2);
  emxInit_real_T(&r6, 2);
  b_emxInit_real_T(&D, 1);
  emxInit_real_T(&D_track, 2);
  b_emxInit_real_T(&b_D, 1);
  b_emxInit_real_T(&c_D, 1);
  b_emxInit_real_T(&d_D, 1);
  exitg1 = FALSE;
  while ((exitg1 == 0U) && (b_check < 20)) {
    check = b_check + 1;

    /*          sprintf('Interior point iteration: %d', check) */
    /*         %% Step 2 */
    for (i4 = 0; i4 < 17; i4++) {
      epsilon = 0.0;
      for (i5 = 0; i5 < 10; i5++) {
        epsilon += A[i4 + 17 * i5] * x[i5];
      }

      w[i4] = b[i4] - epsilon;
    }

    /*         %% Step 3 */
    for (i = 0; i < 17; i++) {
      e_D[i] = 1.0 / w[i];
    }

    /* inverse of slack, therefore the following D's are ALL inverses */
    i4 = d_D->size[0];
    d_D->size[0] = (int32_T)row_v_data[0];
    emxEnsureCapacity((emxArray__common *)d_D, i4, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)row_v_data[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      d_D->data[i4] = e_D[i4];
    }

    diag(d_D, Dtreetrack);
    D_tree_size[0] = Dtreetrack->size[0];
    D_tree_size[1] = Dtreetrack->size[1];
    loop_ub = Dtreetrack->size[0] * Dtreetrack->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      D_tree_data[i4] = Dtreetrack->data[i4];
    }

    if (row_v_data[0] + 1.0 > row_nonneg1_data[0] - 1.0) {
      i4 = 0;
      i5 = -1;
    } else {
      i4 = (int32_T)row_v_data[0];
      i5 = (int32_T)row_nonneg1_data[0] - 2;
    }

    i6 = c_D->size[0];
    c_D->size[0] = (i5 - i4) + 1;
    emxEnsureCapacity((emxArray__common *)c_D, i6, (int32_T)sizeof(real_T));
    loop_ub = i5 - i4;
    for (i5 = 0; i5 <= loop_ub; i5++) {
      c_D->data[i5] = e_D[i4 + i5];
    }

    diag(c_D, Dtreetrack);
    D_track_size[0] = Dtreetrack->size[0];
    D_track_size[1] = Dtreetrack->size[1];
    loop_ub = Dtreetrack->size[0] * Dtreetrack->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      D_track_data[i4] = Dtreetrack->data[i4];
    }

    if (row_nonneg1_data[0] > row_nonneg2_data[0]) {
      i4 = 0;
      i5 = 0;
    } else {
      i4 = (int32_T)row_nonneg1_data[0] - 1;
      i5 = (int32_T)row_nonneg2_data[0];
    }

    i6 = b_D->size[0];
    b_D->size[0] = i5 - i4;
    emxEnsureCapacity((emxArray__common *)b_D, i6, (int32_T)sizeof(real_T));
    loop_ub = (i5 - i4) - 1;
    for (i6 = 0; i6 <= loop_ub; i6++) {
      b_D->data[i6] = e_D[i4 + i6];
    }

    diag(b_D, Dtreetrack);
    D_nonneg_size[0] = Dtreetrack->size[0];
    D_nonneg_size[1] = Dtreetrack->size[1];
    loop_ub = Dtreetrack->size[0] * Dtreetrack->size[1] - 1;
    for (i6 = 0; i6 <= loop_ub; i6++) {
      D_nonneg_data[i6] = Dtreetrack->data[i6];
    }

    /* Construct ADDA */
    Atreetrack_size[0] = 10;
    Atreetrack_size[1] = (int32_T)row_v_data[0];
    loop_ub = (int32_T)row_v_data[0] - 1;
    for (i6 = 0; i6 <= loop_ub; i6++) {
      for (i7 = 0; i7 < 10; i7++) {
        Atreetrack_data[i7 + 10 * i6] = A[i6 + 17 * i7];
      }
    }

    if ((Atreetrack_size[1] == 1) || (D_tree_size[0] == 1)) {
      i6 = u->size[0] * u->size[1];
      u->size[0] = 10;
      u->size[1] = D_tree_size[1];
      emxEnsureCapacity((emxArray__common *)u, i6, (int32_T)sizeof(real_T));
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = D_tree_size[1] - 1;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          u->data[i6 + u->size[0] * i7] = 0.0;
          ixstop = Atreetrack_size[1] - 1;
          for (i = 0; i <= ixstop; i++) {
            u->data[i6 + u->size[0] * i7] += Atreetrack_data[i6 + 10 * i] *
              D_tree_data[i + D_tree_size[0] * i7];
          }
        }
      }
    } else {
      outsz[0] = 10U;
      outsz[1] = (uint32_T)D_tree_size[1];
      i6 = u->size[0] * u->size[1];
      u->size[0] = 10;
      u->size[1] = (int32_T)outsz[1];
      u->size[0] = 10;

      emxEnsureCapacity((emxArray__common *)u, i6, (int32_T)sizeof(real_T));
      loop_ub = u->size[1] - 1;
      for (i6 = 0; i6 <= loop_ub; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          u->data[i7 + u->size[0] * i6] = 0.0;
        }
      }

      b_D_tree_data.data = (real_T *)&D_tree_data;
      b_D_tree_data.size = (int32_T *)&D_tree_size;
      b_D_tree_data.allocatedSize = 289;
      b_D_tree_data.numDimensions = 2;
      b_D_tree_data.canFreeData = FALSE;
      eml_xgemm(D_tree_size[1], Atreetrack_size[1], Atreetrack_data,
                Atreetrack_size, &b_D_tree_data, Atreetrack_size[1], u);
    }

    if ((u->size[1] == 1) || (D_tree_size[0] == 1)) {
      i6 = y->size[0] * y->size[1];
      y->size[0] = 10;
      y->size[1] = D_tree_size[1];
      emxEnsureCapacity((emxArray__common *)y, i6, (int32_T)sizeof(real_T));
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = D_tree_size[1] - 1;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          y->data[i6 + y->size[0] * i7] = 0.0;
          ixstop = u->size[1] - 1;
          for (i = 0; i <= ixstop; i++) {
            y->data[i6 + y->size[0] * i7] += u->data[i6 + u->size[0] * i] *
              D_tree_data[i + D_tree_size[0] * i7];
          }
        }
      }
    } else {
      k = u->size[1];
      outsz[0] = 10U;
      outsz[1] = (uint32_T)D_tree_size[1];
      i6 = y->size[0] * y->size[1];
      y->size[0] = 10;
      y->size[1] = (int32_T)outsz[1];
      y->size[0] = 10;
      emxEnsureCapacity((emxArray__common *)y, i6, (int32_T)sizeof(real_T));
      loop_ub = y->size[1] - 1;
      for (i6 = 0; i6 <= loop_ub; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          y->data[i7 + y->size[0] * i6] = 0.0;
        }
      }

      c_D_tree_data.data = (real_T *)&D_tree_data;
      c_D_tree_data.size = (int32_T *)&D_tree_size;
      c_D_tree_data.allocatedSize = 289;
      c_D_tree_data.numDimensions = 2;
      c_D_tree_data.canFreeData = FALSE;
      eml_xgemm(D_tree_size[1], k, u->data, u->size, &c_D_tree_data, k, y);
    }

    if ((y->size[1] == 1) || ((int32_T)row_v_data[0] == 1)) {
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = (int32_T)row_v_data[0] - 1;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          A_data[i7 + (int32_T)row_v_data[0] * i6] = A[i7 + 17 * i6];
        }
      }

      for (i6 = 0; i6 < 10; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          addaTree[i6 + 10 * i7] = 0.0;
          loop_ub = y->size[1] - 1;
          for (i = 0; i <= loop_ub; i++) {
            addaTree[i6 + 10 * i7] += y->data[i6 + y->size[0] * i] * A_data[i +
              (int32_T)row_v_data[0] * i7];
          }
        }
      }
    } else {
      k = y->size[1];
      memset(&addaTree[0], 0, 100U * sizeof(real_T));
      for (ixstop = 0; ixstop < 92; ixstop += 10) {
        for (ic = ixstop; ic + 1 <= ixstop + 10; ic++) {
          addaTree[ic] = 0.0;
        }
      }

      br = 0;
      for (ixstop = 0; ixstop < 92; ixstop += 10) {
        ar = 0;
        i6 = br + k;
        for (ib = br; ib + 1 <= i6; ib++) {
          if (A[ib % (int32_T)row_v_data[0] + 17 * (ib / (int32_T)row_v_data[0])]
              != 0.0) {
            ia = ar;
            for (ic = ixstop; ic + 1 <= ixstop + 10; ic++) {
              ia++;
              addaTree[ic] += A[ib % (int32_T)row_v_data[0] + 17 * (ib /
                (int32_T)row_v_data[0])] * y->data[ia - 1];
            }
          }

          ar += 10;
        }

        br += k;
      }
    }

    Atreetrack_size[0] = 10;
    Atreetrack_size[1] = (i3 - i2) + 1;
    loop_ub = i3 - i2;
    for (i6 = 0; i6 <= loop_ub; i6++) {
      for (i7 = 0; i7 < 10; i7++) {
        Atreetrack_data[i7 + 10 * i6] = A[(i2 + i6) + 17 * i7];
      }
    }

    if ((Atreetrack_size[1] == 1) || (D_track_size[0] == 1)) {
      i6 = u->size[0] * u->size[1];
      u->size[0] = 10;
      u->size[1] = D_track_size[1];
      emxEnsureCapacity((emxArray__common *)u, i6, (int32_T)sizeof(real_T));
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = D_track_size[1] - 1;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          u->data[i6 + u->size[0] * i7] = 0.0;
          ixstop = Atreetrack_size[1] - 1;
          for (i = 0; i <= ixstop; i++) {
            u->data[i6 + u->size[0] * i7] += Atreetrack_data[i6 + 10 * i] *
              D_track_data[i + D_track_size[0] * i7];
          }
        }
      }
    } else {
      outsz[0] = 10U;
      outsz[1] = (uint32_T)D_track_size[1];
      i6 = u->size[0] * u->size[1];
      u->size[0] = 10;
      u->size[1] = (int32_T)outsz[1];
      u->size[0] = 10;

      emxEnsureCapacity((emxArray__common *)u, i6, (int32_T)sizeof(real_T));
      loop_ub = u->size[1] - 1;
      for (i6 = 0; i6 <= loop_ub; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          u->data[i7 + u->size[0] * i6] = 0.0;
        }
      }

      b_D_track_data.data = (real_T *)&D_track_data;
      b_D_track_data.size = (int32_T *)&D_track_size;
      b_D_track_data.allocatedSize = 256;
      b_D_track_data.numDimensions = 2;
      b_D_track_data.canFreeData = FALSE;
      eml_xgemm(D_track_size[1], Atreetrack_size[1], Atreetrack_data,
                Atreetrack_size, &b_D_track_data, Atreetrack_size[1], u);
    }

    if ((u->size[1] == 1) || (D_track_size[0] == 1)) {
      i6 = y->size[0] * y->size[1];
      y->size[0] = 10;
      y->size[1] = D_track_size[1];
      emxEnsureCapacity((emxArray__common *)y, i6, (int32_T)sizeof(real_T));
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = D_track_size[1] - 1;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          y->data[i6 + y->size[0] * i7] = 0.0;
          ixstop = u->size[1] - 1;
          for (i = 0; i <= ixstop; i++) {
            y->data[i6 + y->size[0] * i7] += u->data[i6 + u->size[0] * i] *
              D_track_data[i + D_track_size[0] * i7];
          }
        }
      }
    } else {
      k = u->size[1];
      outsz[0] = 10U;
      outsz[1] = (uint32_T)D_track_size[1];
      i6 = y->size[0] * y->size[1];
      y->size[0] = 10;
      y->size[1] = (int32_T)outsz[1];
      y->size[0] = 10;
      emxEnsureCapacity((emxArray__common *)y, i6, (int32_T)sizeof(real_T));
      loop_ub = y->size[1] - 1;
      for (i6 = 0; i6 <= loop_ub; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          y->data[i7 + y->size[0] * i6] = 0.0;
        }
      }

      c_D_track_data.data = (real_T *)&D_track_data;
      c_D_track_data.size = (int32_T *)&D_track_size;
      c_D_track_data.allocatedSize = 256;
      c_D_track_data.numDimensions = 2;
      c_D_track_data.canFreeData = FALSE;
      eml_xgemm(D_track_size[1], k, u->data, u->size, &c_D_track_data, k, y);
    }

    if ((y->size[1] == 1) || ((i3 - i2) + 1 == 1)) {
      pos = (i3 - i2) + 1;
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = i3 - i2;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          A_data[i7 + pos * i6] = A[(i2 + i7) + 17 * i6];
        }
      }

      for (i6 = 0; i6 < 10; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          addaTrack[i6 + 10 * i7] = 0.0;
          loop_ub = y->size[1] - 1;
          for (i = 0; i <= loop_ub; i++) {
            addaTrack[i6 + 10 * i7] += y->data[i6 + y->size[0] * i] * A_data[i +
              pos * i7];
          }
        }
      }
    } else {
      k = y->size[1];
      memset(&addaTrack[0], 0, 100U * sizeof(real_T));
      for (ixstop = 0; ixstop < 92; ixstop += 10) {
        for (ic = ixstop; ic + 1 <= ixstop + 10; ic++) {
          addaTrack[ic] = 0.0;
        }
      }

      br = 0;
      for (ixstop = 0; ixstop < 92; ixstop += 10) {
        ar = 0;
        i6 = br + k;
        for (ib = br; ib + 1 <= i6; ib++) {
          pos = (i3 - i2) + 1;
          for (i7 = 0; i7 < 10; i7++) {
            loop_ub = i3 - i2;
            for (i = 0; i <= loop_ub; i++) {
              A_data[i + pos * i7] = A[(i2 + i) + 17 * i7];
            }
          }

          if (A_data[ib] != 0.0) {
            ia = ar;
            for (ic = ixstop; ic + 1 <= ixstop + 10; ic++) {
              ia++;
              pos = (i3 - i2) + 1;
              for (i7 = 0; i7 < 10; i7++) {
                loop_ub = i3 - i2;
                for (i = 0; i <= loop_ub; i++) {
                  A_data[i + pos * i7] = A[(i2 + i) + 17 * i7];
                }
              }

              addaTrack[ic] += A_data[ib] * y->data[ia - 1];
            }
          }

          ar += 10;
        }

        br += k;
      }
    }

    Atreetrack_size[0] = 10;
    Atreetrack_size[1] = i1 - i0;
    loop_ub = (i1 - i0) - 1;
    for (i6 = 0; i6 <= loop_ub; i6++) {
      for (i7 = 0; i7 < 10; i7++) {
        Atreetrack_data[i7 + 10 * i6] = A[(i0 + i6) + 17 * i7];
      }
    }

    if ((Atreetrack_size[1] == 1) || (D_nonneg_size[0] == 1)) {
      i6 = u->size[0] * u->size[1];
      u->size[0] = 10;
      u->size[1] = D_nonneg_size[1];
      emxEnsureCapacity((emxArray__common *)u, i6, (int32_T)sizeof(real_T));
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = D_nonneg_size[1] - 1;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          u->data[i6 + u->size[0] * i7] = 0.0;
          ixstop = Atreetrack_size[1] - 1;
          for (i = 0; i <= ixstop; i++) {
            u->data[i6 + u->size[0] * i7] += Atreetrack_data[i6 + 10 * i] *
              D_nonneg_data[i + D_nonneg_size[0] * i7];
          }
        }
      }
    } else {
      outsz[0] = 10U;
      outsz[1] = (uint32_T)D_nonneg_size[1];
      i6 = u->size[0] * u->size[1];
      u->size[0] = 10;
      u->size[1] = (int32_T)outsz[1];
      u->size[0] = 10;
      u->size[1] = u->size[1];
      emxEnsureCapacity((emxArray__common *)u, i6, (int32_T)sizeof(real_T));
      loop_ub = u->size[1] - 1;
      for (i6 = 0; i6 <= loop_ub; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          u->data[i7 + u->size[0] * i6] = 0.0;
        }
      }

      b_D_nonneg_data.data = (real_T *)&D_nonneg_data;
      b_D_nonneg_data.size = (int32_T *)&D_nonneg_size;
      b_D_nonneg_data.allocatedSize = 289;
      b_D_nonneg_data.numDimensions = 2;
      b_D_nonneg_data.canFreeData = FALSE;
      eml_xgemm(D_nonneg_size[1], Atreetrack_size[1], Atreetrack_data,
                Atreetrack_size, &b_D_nonneg_data, Atreetrack_size[1], u);
    }

    if ((u->size[1] == 1) || (D_nonneg_size[0] == 1)) {
      i6 = y->size[0] * y->size[1];
      y->size[0] = 10;
      y->size[1] = D_nonneg_size[1];
      emxEnsureCapacity((emxArray__common *)y, i6, (int32_T)sizeof(real_T));
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = D_nonneg_size[1] - 1;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          y->data[i6 + y->size[0] * i7] = 0.0;
          ixstop = u->size[1] - 1;
          for (i = 0; i <= ixstop; i++) {
            y->data[i6 + y->size[0] * i7] += u->data[i6 + u->size[0] * i] *
              D_nonneg_data[i + D_nonneg_size[0] * i7];
          }
        }
      }
    } else {
      k = u->size[1];
      outsz[0] = 10U;
      outsz[1] = (uint32_T)D_nonneg_size[1];
      i6 = y->size[0] * y->size[1];
      y->size[0] = 10;
      y->size[1] = (int32_T)outsz[1];
      y->size[0] = 10;
      y->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)y, i6, (int32_T)sizeof(real_T));
      loop_ub = y->size[1] - 1;
      for (i6 = 0; i6 <= loop_ub; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          y->data[i7 + y->size[0] * i6] = 0.0;
        }
      }

      c_D_nonneg_data.data = (real_T *)&D_nonneg_data;
      c_D_nonneg_data.size = (int32_T *)&D_nonneg_size;
      c_D_nonneg_data.allocatedSize = 289;
      c_D_nonneg_data.numDimensions = 2;
      c_D_nonneg_data.canFreeData = FALSE;
      eml_xgemm(D_nonneg_size[1], k, u->data, u->size, &c_D_nonneg_data, k, y);
    }

    if ((y->size[1] == 1) || (i1 - i0 == 1)) {
      pos = i1 - i0;
      for (i6 = 0; i6 < 10; i6++) {
        loop_ub = (i1 - i0) - 1;
        for (i7 = 0; i7 <= loop_ub; i7++) {
          A_data[i7 + pos * i6] = A[(i0 + i7) + 17 * i6];
        }
      }

      for (i6 = 0; i6 < 10; i6++) {
        for (i7 = 0; i7 < 10; i7++) {
          addaNonneg[i6 + 10 * i7] = 0.0;
          loop_ub = y->size[1] - 1;
          for (i = 0; i <= loop_ub; i++) {
            addaNonneg[i6 + 10 * i7] += y->data[i6 + y->size[0] * i] * A_data[i
              + pos * i7];
          }
        }
      }
    } else {
      k = y->size[1];
      memset(&addaNonneg[0], 0, 100U * sizeof(real_T));
      for (ixstop = 0; ixstop < 92; ixstop += 10) {
        for (ic = ixstop; ic + 1 <= ixstop + 10; ic++) {
          addaNonneg[ic] = 0.0;
        }
      }

      br = 0;
      for (ixstop = 0; ixstop < 92; ixstop += 10) {
        ar = 0;
        i6 = br + k;
        for (ib = br; ib + 1 <= i6; ib++) {
          if (A[(i0 + ib % (i1 - i0)) + 17 * (ib / (i1 - i0))] != 0.0) {
            ia = ar;
            for (ic = ixstop; ic + 1 <= ixstop + 10; ic++) {
              ia++;
              addaNonneg[ic] += A[(i0 + ib % (i1 - i0)) + 17 * (ib / (i1 - i0))]
                * y->data[ia - 1];
            }
          }

          ar += 10;
        }

        br += k;
      }
    }

    for (i6 = 0; i6 < 100; i6++) {
      addaTree[i6] = (addaTree[i6] + addaTrack[i6]) + addaNonneg[i6];
    }

    ar = 1;
    pos = i5 - i4;
    alpha_2 = e_D[i4];
    if (pos > 1) {
      if (rtIsNaN(e_D[i4])) {
        ib = 2;
        exitg8 = FALSE;
        while ((exitg8 == 0U) && (ib <= pos)) {
          ar = ib;
          if (!rtIsNaN(e_D[(i4 + ib) - 1])) {
            alpha_2 = e_D[(i4 + ib) - 1];
            exitg8 = TRUE;
          } else {
            ib++;
          }
        }
      }

      if (ar < pos) {
        while (ar + 1 <= pos) {
          if (e_D[i4 + ar] > alpha_2) {
            alpha_2 = e_D[i4 + ar];
          }

          ar++;
        }
      }
    }

    epsilon = alpha_2 / 10.0;
    pos = D_track_size[0] * D_track_size[0];
    loop_ub = pos - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      diag_ind_data[i4] = 0.0;
    }

    pos = 0;
    *PCGit = 1.0;
    while (pos + 1 <= D_track_size[0] * D_track_size[0] + 1) {
      diag_ind_data[pos] = *PCGit;
      pos = (pos + D_track_size[0]) + 1;
      (*PCGit)++;
    }

    /* bool_tracks2=D_track>Breakpoint; */
    pos = D_track_size[0] * D_track_size[1];
    for (i4 = 0; i4 < 2; i4++) {
      sz[i4] = 0;
    }

    sz[0] = (int16_T)(D_track_size[0] * D_track_size[1]);
    sz[1] = 1;
    for (k = 0; k + 1 <= pos; k++) {
      tracks_i_data[k] = D_track_data[k];
    }

    /* tracks_i2=diag_ind(bool_tracks2);%provides row and column(omitted here) indicies of values greater than value (each in vector) */
    tracks_i_size[0] = sz[0];
    loop_ub = sz[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      b_tracks_i_data[i4] = (tracks_i_data[i4] > epsilon);
    }

    b_eml_li_find(b_tracks_i_data, tracks_i_size, b_tmp_data, tmp_size);
    loop_ub = tmp_size[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      tracks_i_data[i4] = diag_ind_data[b_tmp_data[i4] - 1];
    }

    /* [tracks_i,~]=find(D_track>Breakpoint); */
    /* test=all(tracks_i2==tracks_i) */
    /* adds the active tracks */
    /* construct u */
    pos = (i3 - i2) + 1;
    for (i4 = 0; i4 < 10; i4++) {
      loop_ub = i3 - i2;
      for (i5 = 0; i5 <= loop_ub; i5++) {
        b_A_data[i5 + pos * i4] = A[(i2 + i5) + 17 * i4];
      }
    }

    Atreetrack_size[0] = 10;
    Atreetrack_size[1] = (int32_T)row_v_data[0] + tmp_size[0];
    loop_ub = (int32_T)row_v_data[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      for (i5 = 0; i5 < 10; i5++) {
        Atreetrack_data[i5 + 10 * i4] = A[i4 + 17 * i5];
      }
    }

    loop_ub = tmp_size[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      for (i5 = 0; i5 < 10; i5++) {
        Atreetrack_data[i5 + 10 * (i4 + (int32_T)row_v_data[0])] = b_A_data
          [((int32_T)tracks_i_data[i4] + pos * i5) - 1];
      }
    }

    i4 = D_track->size[0] * D_track->size[1];
    D_track->size[0] = tmp_size[0];
    D_track->size[1] = tmp_size[0];
    emxEnsureCapacity((emxArray__common *)D_track, i4, (int32_T)sizeof(real_T));
    loop_ub = tmp_size[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      ixstop = tmp_size[0] - 1;
      for (i5 = 0; i5 <= ixstop; i5++) {
        D_track->data[i5 + D_track->size[0] * i4] = D_track_data[((int32_T)
          tracks_i_data[i5] + D_track_size[0] * ((int32_T)tracks_i_data[i4] - 1))
          - 1];
      }
    }

    b_diag(D_track, r1);
    i4 = D_track_tempVec->size[0] * D_track_tempVec->size[1];
    D_track_tempVec->size[0] = r1->size[0];
    D_track_tempVec->size[1] = r1->size[1];
    emxEnsureCapacity((emxArray__common *)D_track_tempVec, i4, (int32_T)sizeof
                      (real_T));
    loop_ub = r1->size[0] * r1->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      D_track_tempVec->data[i4] = r1->data[i4];
    }

    pos = D_track_tempVec->size[0];
    i4 = D->size[0];
    D->size[0] = (int32_T)row_v_data[0] + pos;
    emxEnsureCapacity((emxArray__common *)D, i4, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)row_v_data[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      D->data[i4] = e_D[i4];
    }

    loop_ub = pos - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      D->data[i4 + (int32_T)row_v_data[0]] = D_track_tempVec->data[i4];
    }

    diag(D, Dtreetrack);
    if ((Atreetrack_size[1] == 1) || (Dtreetrack->size[0] == 1)) {
      i4 = u->size[0] * u->size[1];
      u->size[0] = 10;
      u->size[1] = Dtreetrack->size[1];
      emxEnsureCapacity((emxArray__common *)u, i4, (int32_T)sizeof(real_T));
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = Dtreetrack->size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          u->data[i4 + u->size[0] * i5] = 0.0;
          ixstop = Atreetrack_size[1] - 1;
          for (i6 = 0; i6 <= ixstop; i6++) {
            u->data[i4 + u->size[0] * i5] += Atreetrack_data[i4 + 10 * i6] *
              Dtreetrack->data[i6 + Dtreetrack->size[0] * i5];
          }
        }
      }
    } else {
      outsz[0] = 10U;
      outsz[1] = (uint32_T)Dtreetrack->size[1];
      i4 = u->size[0] * u->size[1];
      u->size[0] = 10;
      u->size[1] = (int32_T)outsz[1];
      u->size[0] = 10;

      emxEnsureCapacity((emxArray__common *)u, i4, (int32_T)sizeof(real_T));
      loop_ub = u->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        for (i5 = 0; i5 < 10; i5++) {
          u->data[i5 + u->size[0] * i4] = 0.0;
        }
      }

      eml_xgemm(Dtreetrack->size[1], Atreetrack_size[1], Atreetrack_data,
                Atreetrack_size, Dtreetrack, Atreetrack_size[1], u);
    }

    /* Prep for step 4: */
    if (!(tmp_size[0] == 0)) {
      if ((D_nonneg_size[1] == 1) || (D_nonneg_size[0] == 1)) {
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = D_nonneg_size[0];
        b_y->size[1] = D_nonneg_size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = D_nonneg_size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = D_nonneg_size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i4 + b_y->size[0] * i5] = 0.0;
            pos = D_nonneg_size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              b_y->data[i4 + b_y->size[0] * i5] += D_nonneg_data[i4 +
                D_nonneg_size[0] * i6] * D_nonneg_data[i6 + D_nonneg_size[0] *
                i5];
            }
          }
        }
      } else {
        k = D_nonneg_size[1];
        outsz[0] = (uint32_T)D_nonneg_size[0];
        outsz[1] = (uint32_T)D_nonneg_size[1];
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = (int32_T)outsz[0];
        b_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        m = D_nonneg_size[0];
        ia = D_nonneg_size[1];
        i4 = b_y->size[0] * b_y->size[1];


        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = b_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i5 + b_y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          pos = m * (ia - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              b_y->data[ic] = 0.0;
            }
          }

          br = 0;
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (D_nonneg_data[ib] != 0.0) {
                ia = ar;
                i5 = ixstop + m;
                for (ic = ixstop; ic + 1 <= i5; ic++) {
                  ia++;
                  b_y->data[ic] += D_nonneg_data[ib] * D_nonneg_data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
          }
        }
      }

      b_diag(b_y, r1);
      pos = r1->size[0];
      m = r1->size[1];
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        diagB_data[i4] = r1->data[i4];
      }

      i4 = r6->size[0] * r6->size[1];
      r6->size[0] = pos;
      r6->size[1] = m;
      emxEnsureCapacity((emxArray__common *)r6, i4, (int32_T)sizeof(real_T));
      loop_ub = pos * m - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        r6->data[i4] = 1.0 / diagB_data[i4];
      }

      b_diag(r6, r1);
      Binv_size[0] = r1->size[0];
      Binv_size[1] = r1->size[1];
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        Binv_data[i4] = r1->data[i4];
      }

      i4 = a->size[0] * a->size[1];
      a->size[0] = u->size[1];
      a->size[1] = 10;
      emxEnsureCapacity((emxArray__common *)a, i4, (int32_T)sizeof(real_T));
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = u->size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          a->data[i5 + a->size[0] * i4] = u->data[i4 + u->size[0] * i5];
        }
      }

      outsz[0] = (uint32_T)a->size[0];
      outsz[1] = (uint32_T)u->size[1];
      i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
      Dtreetrack->size[0] = (int32_T)outsz[0];
      Dtreetrack->size[1] = (int32_T)outsz[1];
      emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                        (real_T));
      m = a->size[0];
      i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
      emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                        (real_T));
      loop_ub = Dtreetrack->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstop = Dtreetrack->size[0] - 1;
        for (i5 = 0; i5 <= ixstop; i5++) {
          Dtreetrack->data[i5 + Dtreetrack->size[0] * i4] = 0.0;
        }
      }

      pos = m * (u->size[1] - 1);
      for (ixstop = 0; ixstop <= pos; ixstop += m) {
        i4 = ixstop + m;
        for (ic = ixstop; ic + 1 <= i4; ic++) {
          Dtreetrack->data[ic] = 0.0;
        }
      }

      br = 0;
      for (ixstop = 0; ixstop <= pos; ixstop += m) {
        ar = 0;
        for (ib = br; ib + 1 <= br + 10; ib++) {
          if (u->data[ib] != 0.0) {
            ia = ar;
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              ia++;
              Dtreetrack->data[ic] += u->data[ib] * a->data[ia - 1];
            }
          }

          ar += m;
        }

        br += 10;
      }

      pos = Dtreetrack->size[0];
      eye((real_T)pos, Dtreetrack);
      i4 = a->size[0] * a->size[1];
      a->size[0] = u->size[1];
      a->size[1] = 10;
      emxEnsureCapacity((emxArray__common *)a, i4, (int32_T)sizeof(real_T));
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = u->size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          a->data[i5 + a->size[0] * i4] = u->data[i4 + u->size[0] * i5];
        }
      }

      if (Binv_size[0] == 1) {
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = a->size[0];
        c_y->size[1] = Binv_size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        loop_ub = a->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Binv_size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            c_y->data[i4 + c_y->size[0] * i5] = 0.0;
            for (i6 = 0; i6 < 10; i6++) {
              c_y->data[i4 + c_y->size[0] * i5] += a->data[i4 + a->size[0] * i6]
                * Binv_data[i6 + i5];
            }
          }
        }
      } else {
        outsz[0] = (uint32_T)a->size[0];
        outsz[1] = (uint32_T)Binv_size[1];
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = (int32_T)outsz[0];
        c_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        m = a->size[0];
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = c_y->size[0];
        c_y->size[1] = c_y->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        loop_ub = c_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = c_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            c_y->data[i5 + c_y->size[0] * i4] = 0.0;
          }
        }

        b_eml_xgemm(m, Binv_size[1], a, m, Binv_data, Binv_size, c_y, m);
      }

      if (c_y->size[1] == 1) {
        i4 = C->size[0] * C->size[1];
        C->size[0] = c_y->size[0];
        C->size[1] = u->size[1];
        emxEnsureCapacity((emxArray__common *)C, i4, (int32_T)sizeof(real_T));
        loop_ub = c_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = u->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            C->data[i4 + C->size[0] * i5] = 0.0;
            pos = c_y->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              C->data[i4 + C->size[0] * i5] += c_y->data[i4 + c_y->size[0] * i6]
                * u->data[i6 + u->size[0] * i5];
            }
          }
        }
      } else {
        outsz[0] = (uint32_T)c_y->size[0];
        outsz[1] = (uint32_T)u->size[1];
        i4 = C->size[0] * C->size[1];
        C->size[0] = (int32_T)outsz[0];
        C->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)C, i4, (int32_T)sizeof(real_T));
        m = c_y->size[0];
        i4 = C->size[0] * C->size[1];
        C->size[0] = C->size[0];
        C->size[1] = C->size[1];
        emxEnsureCapacity((emxArray__common *)C, i4, (int32_T)sizeof(real_T));
        loop_ub = C->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = C->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            C->data[i5 + C->size[0] * i4] = 0.0;
          }
        }

        if (m == 0) {
        } else {
          pos = m * (u->size[1] - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              C->data[ic] = 0.0;
            }
          }
        }
      }

      i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
      emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                        (real_T));
      pos = Dtreetrack->size[0];
      m = Dtreetrack->size[1];
      loop_ub = pos * m - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        Dtreetrack->data[i4] += C->data[i4];
      }

      /* Break appart initial arrow% */
      /* The dimensions of the diagonal portion are euqal to the number of rows that create the trees in A */
      i4 = c_Dtreetrack->size[0] * c_Dtreetrack->size[1];
      c_Dtreetrack->size[0] = (int32_T)row_v_data[0];
      c_Dtreetrack->size[1] = (int32_T)row_v_data[0];
      emxEnsureCapacity((emxArray__common *)c_Dtreetrack, i4, (int32_T)sizeof
                        (real_T));
      loop_ub = (int32_T)row_v_data[0] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstop = (int32_T)row_v_data[0] - 1;
        for (i5 = 0; i5 <= ixstop; i5++) {
          c_Dtreetrack->data[i5 + c_Dtreetrack->size[0] * i4] = Dtreetrack->
            data[i5 + Dtreetrack->size[0] * i4];
        }
      }

      b_diag(c_Dtreetrack, r1);
      i4 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = r1->size[0];
      d_y->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        d_y->data[i4] = r1->data[i4];
      }

      i4 = r5->size[0] * r5->size[1];
      r5->size[0] = d_y->size[0];
      r5->size[1] = d_y->size[1];
      emxEnsureCapacity((emxArray__common *)r5, i4, (int32_T)sizeof(real_T));
      loop_ub = d_y->size[0] * d_y->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        r5->data[i4] = 1.0 / d_y->data[i4];
      }

      b_diag(r5, r1);
      i4 = A_11inv->size[0] * A_11inv->size[1];
      A_11inv->size[0] = r1->size[0];
      A_11inv->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)A_11inv, i4, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        A_11inv->data[i4] = r1->data[i4];
      }

      if ((int32_T)row_v_data[0] + 1 > Dtreetrack->size[0]) {
        i4 = 0;
        i5 = 0;
      } else {
        i4 = (int32_T)row_v_data[0];
        i5 = Dtreetrack->size[0];
      }

      i6 = A_12->size[0] * A_12->size[1];
      A_12->size[0] = (int32_T)row_v_data[0];
      A_12->size[1] = i5 - i4;
      emxEnsureCapacity((emxArray__common *)A_12, i6, (int32_T)sizeof(real_T));
      loop_ub = (i5 - i4) - 1;
      for (i5 = 0; i5 <= loop_ub; i5++) {
        ixstop = (int32_T)row_v_data[0] - 1;
        for (i6 = 0; i6 <= ixstop; i6++) {
          A_12->data[i6 + A_12->size[0] * i5] = Dtreetrack->data[(i4 + i5) +
            Dtreetrack->size[0] * i6];
        }
      }

      if ((int32_T)row_v_data[0] + 1 > Dtreetrack->size[0]) {
        i4 = 0;
        i5 = 0;
      } else {
        i4 = (int32_T)row_v_data[0];
        i5 = Dtreetrack->size[0];
      }

      if ((int32_T)row_v_data[0] + 1 > Dtreetrack->size[1]) {
        i6 = 0;
        i7 = 0;
      } else {
        i6 = (int32_T)row_v_data[0];
        i7 = Dtreetrack->size[1];
      }

      i = A_12_Psi->size[0] * A_12_Psi->size[1];
      A_12_Psi->size[0] = A_12->size[0];
      A_12_Psi->size[1] = A_12->size[1];
      emxEnsureCapacity((emxArray__common *)A_12_Psi, i, (int32_T)sizeof(real_T));
      loop_ub = A_12->size[0] * A_12->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        A_12_Psi->data[i] = A_12->data[i];
      }

      for (i = 0; i < 2; i++) {
        outsz[i] = (uint32_T)A_12->size[i];
      }

      outsz[0] = 1U;
      i = maxval->size[0] * maxval->size[1];
      maxval->size[0] = 1;
      maxval->size[1] = (int32_T)outsz[1];
      emxEnsureCapacity((emxArray__common *)maxval, i, (int32_T)sizeof(real_T));
      ia = A_12->size[0];
      pos = A_12->size[1];
      ib = 0;
      m = -1;
      for (i = 1; i <= pos; i++) {
        ar = ib;
        ixstop = ib + ia;
        alpha_2 = A_12->data[ib];
        if (ia > 1) {
          if (rtIsNaN(A_12->data[ib])) {
            br = ib + 1;
            exitg7 = FALSE;
            while ((exitg7 == 0U) && (br + 1 <= ixstop)) {
              ar = br;
              if (!rtIsNaN(A_12->data[br])) {
                alpha_2 = A_12->data[br];
                exitg7 = TRUE;
              } else {
                br++;
              }
            }
          }

          if (ar + 1 < ixstop) {
            for (br = ar + 1; br + 1 <= ixstop; br++) {
              if (A_12->data[br] > alpha_2) {
                alpha_2 = A_12->data[br];
              }
            }
          }
        }

        m++;
        maxval->data[m] = alpha_2;
        ib += ia;
      }

      ar = 1;
      ia = maxval->size[1];
      alpha_2 = maxval->data[0];
      if (ia > 1) {
        if (rtIsNaN(maxval->data[0])) {
          ib = 2;
          exitg6 = FALSE;
          while ((exitg6 == 0U) && (ib <= ia)) {
            ar = ib;
            if (!rtIsNaN(maxval->data[ib - 1])) {
              alpha_2 = maxval->data[ib - 1];
              exitg6 = TRUE;
            } else {
              ib++;
            }
          }
        }

        if (ar < ia) {
          while (ar + 1 <= ia) {
            if (maxval->data[ar] > alpha_2) {
              alpha_2 = maxval->data[ar];
            }

            ar++;
          }
        }
      }

      for (i = 0; i < 2; i++) {
        outsz[i] = (uint32_T)A_12->size[i];
      }

      outsz[0] = 1U;
      i = maxval->size[0] * maxval->size[1];
      maxval->size[0] = 1;
      maxval->size[1] = (int32_T)outsz[1];
      emxEnsureCapacity((emxArray__common *)maxval, i, (int32_T)sizeof(real_T));
      ia = A_12->size[0];
      pos = A_12->size[1];
      ib = 0;
      m = -1;
      for (i = 1; i <= pos; i++) {
        ar = ib;
        ixstop = ib + ia;
        epsilon = A_12->data[ib];
        if (ia > 1) {
          if (rtIsNaN(A_12->data[ib])) {
            br = ib + 1;
            exitg5 = FALSE;
            while ((exitg5 == 0U) && (br + 1 <= ixstop)) {
              ar = br;
              if (!rtIsNaN(A_12->data[br])) {
                epsilon = A_12->data[br];
                exitg5 = TRUE;
              } else {
                br++;
              }
            }
          }

          if (ar + 1 < ixstop) {
            for (br = ar + 1; br + 1 <= ixstop; br++) {
              if (A_12->data[br] < epsilon) {
                epsilon = A_12->data[br];
              }
            }
          }
        }

        m++;
        maxval->data[m] = epsilon;
        ib += ia;
      }

      ar = 1;
      ia = maxval->size[1];
      epsilon = maxval->data[0];
      if (ia > 1) {
        if (rtIsNaN(maxval->data[0])) {
          ib = 2;
          exitg4 = FALSE;
          while ((exitg4 == 0U) && (ib <= ia)) {
            ar = ib;
            if (!rtIsNaN(maxval->data[ib - 1])) {
              epsilon = maxval->data[ib - 1];
              exitg4 = TRUE;
            } else {
              ib++;
            }
          }
        }

        if (ar < ia) {
          while (ar + 1 <= ia) {
            if (maxval->data[ar] < epsilon) {
              epsilon = maxval->data[ar];
            }

            ar++;
          }
        }
      }

      f_y = fabs(alpha_2 - epsilon) / 100.0;
      b_abs(A_12, r2);
      i = A_12_Psi_to_remove->size[0] * A_12_Psi_to_remove->size[1];
      A_12_Psi_to_remove->size[0] = r2->size[0];
      A_12_Psi_to_remove->size[1] = r2->size[1];
      emxEnsureCapacity((emxArray__common *)A_12_Psi_to_remove, i, (int32_T)
                        sizeof(boolean_T));
      loop_ub = r2->size[0] * r2->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        A_12_Psi_to_remove->data[i] = (r2->data[i] < f_y);
      }

      c_eml_li_find(A_12_Psi_to_remove, r3);
      i = r0->size[0];
      r0->size[0] = r3->size[0];
      emxEnsureCapacity((emxArray__common *)r0, i, (int32_T)sizeof(int32_T));
      loop_ub = r3->size[0] - 1;
      for (i = 0; i <= loop_ub; i++) {
        r0->data[i] = r3->data[i];
      }

      loop_ub = r0->size[0] - 1;
      for (i = 0; i <= loop_ub; i++) {
        A_12_Psi->data[r0->data[i] - 1] = 0.0;
      }

      /* Psi is the difficult portion to inverse */
      i = b_Dtreetrack->size[0] * b_Dtreetrack->size[1];
      b_Dtreetrack->size[0] = (int32_T)row_v_data[0];
      b_Dtreetrack->size[1] = (int32_T)row_v_data[0];
      emxEnsureCapacity((emxArray__common *)b_Dtreetrack, i, (int32_T)sizeof
                        (real_T));
      loop_ub = (int32_T)row_v_data[0] - 1;
      for (i = 0; i <= loop_ub; i++) {
        ixstop = (int32_T)row_v_data[0] - 1;
        for (br = 0; br <= ixstop; br++) {
          b_Dtreetrack->data[br + b_Dtreetrack->size[0] * i] = Dtreetrack->
            data[br + Dtreetrack->size[0] * i];
        }
      }

      b_diag(b_Dtreetrack, r1);
      i = d_y->size[0] * d_y->size[1];
      d_y->size[0] = r1->size[0];
      d_y->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)d_y, i, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        d_y->data[i] = r1->data[i];
      }

      i = Arrow_inv->size[0] * Arrow_inv->size[1];
      Arrow_inv->size[0] = A_12_Psi->size[1];
      Arrow_inv->size[1] = A_12_Psi->size[0];
      emxEnsureCapacity((emxArray__common *)Arrow_inv, i, (int32_T)sizeof(real_T));
      loop_ub = A_12_Psi->size[0] - 1;
      for (i = 0; i <= loop_ub; i++) {
        ixstop = A_12_Psi->size[1] - 1;
        for (br = 0; br <= ixstop; br++) {
          Arrow_inv->data[br + Arrow_inv->size[0] * i] = A_12_Psi->data[i +
            A_12_Psi->size[0] * br];
        }
      }

      i = r4->size[0] * r4->size[1];
      r4->size[0] = d_y->size[0];
      r4->size[1] = d_y->size[1];
      emxEnsureCapacity((emxArray__common *)r4, i, (int32_T)sizeof(real_T));
      loop_ub = d_y->size[0] * d_y->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        r4->data[i] = 1.0 / d_y->data[i];
      }

      b_diag(r4, r1);
      i = Psi_notinv->size[0] * Psi_notinv->size[1];
      Psi_notinv->size[0] = r1->size[0];
      Psi_notinv->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)Psi_notinv, i, (int32_T)sizeof
                        (real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        Psi_notinv->data[i] = r1->data[i];
      }

      if ((Arrow_inv->size[1] == 1) || (Psi_notinv->size[0] == 1)) {
        i = b_y->size[0] * b_y->size[1];
        b_y->size[0] = Arrow_inv->size[0];
        b_y->size[1] = Psi_notinv->size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i, (int32_T)sizeof(real_T));
        loop_ub = Arrow_inv->size[0] - 1;
        for (i = 0; i <= loop_ub; i++) {
          ixstop = Psi_notinv->size[1] - 1;
          for (br = 0; br <= ixstop; br++) {
            b_y->data[i + b_y->size[0] * br] = 0.0;
            pos = Arrow_inv->size[1] - 1;
            for (m = 0; m <= pos; m++) {
              b_y->data[i + b_y->size[0] * br] += Arrow_inv->data[i +
                Arrow_inv->size[0] * m] * Psi_notinv->data[m + Psi_notinv->size
                [0] * br];
            }
          }
        }
      } else {
        k = Arrow_inv->size[1];
        outsz[0] = (uint32_T)Arrow_inv->size[0];
        outsz[1] = (uint32_T)Psi_notinv->size[1];
        i = b_y->size[0] * b_y->size[1];
        b_y->size[0] = (int32_T)outsz[0];
        b_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)b_y, i, (int32_T)sizeof(real_T));
        m = Arrow_inv->size[0];
        ia = Psi_notinv->size[1];
        i = b_y->size[0] * b_y->size[1];


        emxEnsureCapacity((emxArray__common *)b_y, i, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          ixstop = b_y->size[0] - 1;
          for (br = 0; br <= ixstop; br++) {
            b_y->data[br + b_y->size[0] * i] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          ixstop = 0;
          while (ixstop <= 0) {
            for (ic = 1; ic <= m; ic++) {
              b_y->data[ic - 1] = 0.0;
            }

            ixstop = m;
          }

          br = 0;
          ixstop = 0;
          while (ixstop <= 0) {
            ar = 0;
            i = br + k;
            for (ib = br; ib + 1 <= i; ib++) {
              if (Psi_notinv->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= m; ic++) {
                  ia++;
                  b_y->data[ic] += Psi_notinv->data[ib] * Arrow_inv->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
            ixstop = m;
          }
        }
      }

      if ((b_y->size[1] == 1) || (A_12_Psi->size[0] == 1)) {
        i = C->size[0] * C->size[1];
        C->size[0] = b_y->size[0];
        C->size[1] = A_12_Psi->size[1];
        emxEnsureCapacity((emxArray__common *)C, i, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[0] - 1;
        for (i = 0; i <= loop_ub; i++) {
          ixstop = A_12_Psi->size[1] - 1;
          for (br = 0; br <= ixstop; br++) {
            C->data[i + C->size[0] * br] = 0.0;
            pos = b_y->size[1] - 1;
            for (m = 0; m <= pos; m++) {
              C->data[i + C->size[0] * br] += b_y->data[i + b_y->size[0] * m] *
                A_12_Psi->data[m + A_12_Psi->size[0] * br];
            }
          }
        }
      } else {
        outsz[0] = (uint32_T)b_y->size[0];
        outsz[1] = (uint32_T)A_12_Psi->size[1];
        i = C->size[0] * C->size[1];
        C->size[0] = (int32_T)outsz[0];
        C->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)C, i, (int32_T)sizeof(real_T));
        m = b_y->size[0];
        ia = A_12_Psi->size[1];
        i = C->size[0] * C->size[1];
        C->size[0] = C->size[0];
        C->size[1] = C->size[1];
        emxEnsureCapacity((emxArray__common *)C, i, (int32_T)sizeof(real_T));
        loop_ub = C->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          ixstop = C->size[0] - 1;
          for (br = 0; br <= ixstop; br++) {
            C->data[br + C->size[0] * i] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          pos = m * (ia - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i = ixstop + m;
            for (ic = ixstop; ic + 1 <= i; ic++) {
              C->data[ic] = 0.0;
            }
          }
        }
      }

      i = Psi_notinv->size[0] * Psi_notinv->size[1];
      Psi_notinv->size[0] = i5 - i4;
      Psi_notinv->size[1] = i7 - i6;
      emxEnsureCapacity((emxArray__common *)Psi_notinv, i, (int32_T)sizeof
                        (real_T));
      loop_ub = (i7 - i6) - 1;
      for (i7 = 0; i7 <= loop_ub; i7++) {
        ixstop = (i5 - i4) - 1;
        for (i = 0; i <= ixstop; i++) {
          Psi_notinv->data[i + Psi_notinv->size[0] * i7] = Dtreetrack->data[(i4
            + i) + Dtreetrack->size[0] * (i6 + i7)] - C->data[i + C->size[0] *
            i7];
        }
      }

      eye((real_T)Psi_notinv->size[1], Dtreetrack);
      mrdivide(Dtreetrack, Psi_notinv, Psi);

      /* Try to avoid the back slash! */
      /* Inversion of Arrow% */
      if ((A_11inv->size[1] == 1) || (A_12->size[0] == 1)) {
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = A_11inv->size[0];
        b_y->size[1] = A_12->size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = A_11inv->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = A_12->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i4 + b_y->size[0] * i5] = 0.0;
            pos = A_11inv->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              b_y->data[i4 + b_y->size[0] * i5] += A_11inv->data[i4 +
                A_11inv->size[0] * i6] * A_12->data[i6 + A_12->size[0] * i5];
            }
          }
        }
      } else {
        k = A_11inv->size[1];
        outsz[0] = (uint32_T)A_11inv->size[0];
        outsz[1] = (uint32_T)A_12->size[1];
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = (int32_T)outsz[0];
        b_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        m = A_11inv->size[0];
        ia = A_12->size[1];
        i4 = b_y->size[0] * b_y->size[1];


        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = b_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i5 + b_y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          pos = m * (ia - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              b_y->data[ic] = 0.0;
            }
          }

          br = 0;
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (A_12->data[ib] != 0.0) {
                ia = ar;
                i5 = ixstop + m;
                for (ic = ixstop; ic + 1 <= i5; ic++) {
                  ia++;
                  b_y->data[ic] += A_12->data[ib] * A_11inv->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
          }
        }
      }

      if ((b_y->size[1] == 1) || (Psi->size[0] == 1)) {
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = b_y->size[0];
        Dtreetrack->size[1] = Psi->size[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        loop_ub = b_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Psi->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            Dtreetrack->data[i4 + Dtreetrack->size[0] * i5] = 0.0;
            pos = b_y->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              Dtreetrack->data[i4 + Dtreetrack->size[0] * i5] += b_y->data[i4 +
                b_y->size[0] * i6] * Psi->data[i6 + Psi->size[0] * i5];
            }
          }
        }
      } else {
        k = b_y->size[1];
        outsz[0] = (uint32_T)b_y->size[0];
        outsz[1] = (uint32_T)Psi->size[1];
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = (int32_T)outsz[0];
        Dtreetrack->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        m = b_y->size[0];
        ia = Psi->size[1];
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = Dtreetrack->size[0];
        Dtreetrack->size[1] = Dtreetrack->size[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        loop_ub = Dtreetrack->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Dtreetrack->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            Dtreetrack->data[i5 + Dtreetrack->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          pos = m * (ia - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              Dtreetrack->data[ic] = 0.0;
            }
          }

          br = 0;
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (Psi->data[ib] != 0.0) {
                ia = ar;
                i5 = ixstop + m;
                for (ic = ixstop; ic + 1 <= i5; ic++) {
                  ia++;
                  Dtreetrack->data[ic] += Psi->data[ib] * b_y->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
          }
        }
      }

      i4 = Psi_notinv->size[0] * Psi_notinv->size[1];
      Psi_notinv->size[0] = A_12->size[1];
      Psi_notinv->size[1] = A_12->size[0];
      emxEnsureCapacity((emxArray__common *)Psi_notinv, i4, (int32_T)sizeof
                        (real_T));
      loop_ub = A_12->size[0] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstop = A_12->size[1] - 1;
        for (i5 = 0; i5 <= ixstop; i5++) {
          Psi_notinv->data[i5 + Psi_notinv->size[0] * i4] = A_12->data[i4 +
            A_12->size[0] * i5];
        }
      }

      if ((Dtreetrack->size[1] == 1) || (Psi_notinv->size[0] == 1)) {
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = Dtreetrack->size[0];
        b_y->size[1] = Psi_notinv->size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = Dtreetrack->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Psi_notinv->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i4 + b_y->size[0] * i5] = 0.0;
            pos = Dtreetrack->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              b_y->data[i4 + b_y->size[0] * i5] += Dtreetrack->data[i4 +
                Dtreetrack->size[0] * i6] * Psi_notinv->data[i6 +
                Psi_notinv->size[0] * i5];
            }
          }
        }
      } else {
        k = Dtreetrack->size[1];
        outsz[0] = (uint32_T)Dtreetrack->size[0];
        outsz[1] = (uint32_T)Psi_notinv->size[1];
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = (int32_T)outsz[0];
        b_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        m = Dtreetrack->size[0];
        ia = Psi_notinv->size[1];
        i4 = b_y->size[0] * b_y->size[1];

        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = b_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i5 + b_y->size[0] * i4] = 0.0;
          }
        }

        if (m == 0) {
        } else {
          pos = m * (ia - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              b_y->data[ic] = 0.0;
            }
          }

          br = 0;
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (Psi_notinv->data[ib] != 0.0) {
                ia = ar;
                i5 = ixstop + m;
                for (ic = ixstop; ic + 1 <= i5; ic++) {
                  ia++;
                  b_y->data[ic] += Psi_notinv->data[ib] * Dtreetrack->data[ia -
                    1];
                }
              }

              ar += m;
            }

            br += k;
          }
        }
      }

      if ((b_y->size[1] == 1) || (A_11inv->size[0] == 1)) {
        i4 = C->size[0] * C->size[1];
        C->size[0] = b_y->size[0];
        C->size[1] = A_11inv->size[1];
        emxEnsureCapacity((emxArray__common *)C, i4, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = A_11inv->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            C->data[i4 + C->size[0] * i5] = 0.0;
            pos = b_y->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              C->data[i4 + C->size[0] * i5] += b_y->data[i4 + b_y->size[0] * i6]
                * A_11inv->data[i6 + A_11inv->size[0] * i5];
            }
          }
        }
      } else {
        k = b_y->size[1];
        outsz[0] = (uint32_T)b_y->size[0];
        outsz[1] = (uint32_T)A_11inv->size[1];
        i4 = C->size[0] * C->size[1];
        C->size[0] = (int32_T)outsz[0];
        C->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)C, i4, (int32_T)sizeof(real_T));
        m = b_y->size[0];
        ia = A_11inv->size[1];
        i4 = C->size[0] * C->size[1];
        C->size[0] = C->size[0];
        C->size[1] = C->size[1];
        emxEnsureCapacity((emxArray__common *)C, i4, (int32_T)sizeof(real_T));
        loop_ub = C->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = C->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            C->data[i5 + C->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          ixstop = 0;
          while (ixstop <= 0) {
            for (ic = 1; ic <= m; ic++) {
              C->data[ic - 1] = 0.0;
            }

            ixstop = m;
          }

          br = 0;
          ixstop = 0;
          while (ixstop <= 0) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (A_11inv->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= m; ic++) {
                  ia++;
                  C->data[ic] += A_11inv->data[ib] * b_y->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
            ixstop = m;
          }
        }
      }

      i4 = Arrow_inv->size[0] * Arrow_inv->size[1];
      Arrow_inv->size[0] = A_11inv->size[0];
      Arrow_inv->size[1] = A_11inv->size[1];
      emxEnsureCapacity((emxArray__common *)Arrow_inv, i4, (int32_T)sizeof
                        (real_T));
      loop_ub = A_11inv->size[0] * A_11inv->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        Arrow_inv->data[i4] = -A_11inv->data[i4];
      }

      if ((Arrow_inv->size[1] == 1) || (A_12->size[0] == 1)) {
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = Arrow_inv->size[0];
        b_y->size[1] = A_12->size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = Arrow_inv->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = A_12->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i4 + b_y->size[0] * i5] = 0.0;
            pos = Arrow_inv->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              b_y->data[i4 + b_y->size[0] * i5] += Arrow_inv->data[i4 +
                Arrow_inv->size[0] * i6] * A_12->data[i6 + A_12->size[0] * i5];
            }
          }
        }
      } else {
        outsz[0] = (uint32_T)Arrow_inv->size[0];
        outsz[1] = (uint32_T)A_12->size[1];
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = (int32_T)outsz[0];
        b_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        m = Arrow_inv->size[0];
        ia = A_12->size[1];
        i4 = b_y->size[0] * b_y->size[1];

        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = b_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i5 + b_y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          pos = m * (ia - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              b_y->data[ic] = 0.0;
            }
          }
        }
      }

      if ((b_y->size[1] == 1) || (Psi->size[0] == 1)) {
        i4 = Psi_notinv->size[0] * Psi_notinv->size[1];
        Psi_notinv->size[0] = b_y->size[0];
        Psi_notinv->size[1] = Psi->size[1];
        emxEnsureCapacity((emxArray__common *)Psi_notinv, i4, (int32_T)sizeof
                          (real_T));
        loop_ub = b_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Psi->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            Psi_notinv->data[i4 + Psi_notinv->size[0] * i5] = 0.0;
            pos = b_y->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              Psi_notinv->data[i4 + Psi_notinv->size[0] * i5] += b_y->data[i4 +
                b_y->size[0] * i6] * Psi->data[i6 + Psi->size[0] * i5];
            }
          }
        }
      } else {
        k = b_y->size[1];
        outsz[0] = (uint32_T)b_y->size[0];
        outsz[1] = (uint32_T)Psi->size[1];
        i4 = Psi_notinv->size[0] * Psi_notinv->size[1];
        Psi_notinv->size[0] = (int32_T)outsz[0];
        Psi_notinv->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)Psi_notinv, i4, (int32_T)sizeof
                          (real_T));
        m = b_y->size[0];
        ia = Psi->size[1];
        i4 = Psi_notinv->size[0] * Psi_notinv->size[1];
        Psi_notinv->size[0] = Psi_notinv->size[0];
        Psi_notinv->size[1] = Psi_notinv->size[1];
        emxEnsureCapacity((emxArray__common *)Psi_notinv, i4, (int32_T)sizeof
                          (real_T));
        loop_ub = Psi_notinv->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Psi_notinv->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            Psi_notinv->data[i5 + Psi_notinv->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          pos = m * (ia - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              Psi_notinv->data[ic] = 0.0;
            }
          }

          br = 0;
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (Psi->data[ib] != 0.0) {
                ia = ar;
                i5 = ixstop + m;
                for (ic = ixstop; ic + 1 <= i5; ic++) {
                  ia++;
                  Psi_notinv->data[ic] += Psi->data[ib] * b_y->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
          }
        }
      }

      i4 = Arrow_inv->size[0] * Arrow_inv->size[1];
      Arrow_inv->size[0] = A_12->size[1];
      Arrow_inv->size[1] = A_12->size[0];
      emxEnsureCapacity((emxArray__common *)Arrow_inv, i4, (int32_T)sizeof
                        (real_T));
      loop_ub = A_12->size[0] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstop = A_12->size[1] - 1;
        for (i5 = 0; i5 <= ixstop; i5++) {
          Arrow_inv->data[i5 + Arrow_inv->size[0] * i4] = -A_12->data[i4 +
            A_12->size[0] * i5];
        }
      }

      if ((Arrow_inv->size[1] == 1) || (A_11inv->size[0] == 1)) {
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = Arrow_inv->size[0];
        b_y->size[1] = A_11inv->size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = Arrow_inv->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = A_11inv->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i4 + b_y->size[0] * i5] = 0.0;
            pos = Arrow_inv->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              b_y->data[i4 + b_y->size[0] * i5] += Arrow_inv->data[i4 +
                Arrow_inv->size[0] * i6] * A_11inv->data[i6 + A_11inv->size[0] *
                i5];
            }
          }
        }
      } else {
        k = Arrow_inv->size[1];
        outsz[0] = (uint32_T)Arrow_inv->size[0];
        outsz[1] = (uint32_T)A_11inv->size[1];
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = (int32_T)outsz[0];
        b_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        m = Arrow_inv->size[0];
        ia = A_11inv->size[1];
        i4 = b_y->size[0] * b_y->size[1];

        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = b_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i5 + b_y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          ixstop = 0;
          while (ixstop <= 0) {
            for (ic = 1; ic <= m; ic++) {
              b_y->data[ic - 1] = 0.0;
            }

            ixstop = m;
          }

          br = 0;
          ixstop = 0;
          while (ixstop <= 0) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (A_11inv->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= m; ic++) {
                  ia++;
                  b_y->data[ic] += A_11inv->data[ib] * Arrow_inv->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
            ixstop = m;
          }
        }
      }

      if ((Psi->size[1] == 1) || (b_y->size[0] == 1)) {
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = Psi->size[0];
        Dtreetrack->size[1] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        loop_ub = Psi->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = b_y->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            Dtreetrack->data[i4 + Dtreetrack->size[0] * i5] = 0.0;
            pos = Psi->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              Dtreetrack->data[i4 + Dtreetrack->size[0] * i5] += Psi->data[i4 +
                Psi->size[0] * i6] * b_y->data[i6 + b_y->size[0] * i5];
            }
          }
        }
      } else {
        k = Psi->size[1];
        outsz[0] = (uint32_T)Psi->size[0];
        outsz[1] = (uint32_T)b_y->size[1];
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = (int32_T)outsz[0];
        Dtreetrack->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        m = Psi->size[0];
        ia = b_y->size[1];
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = Dtreetrack->size[0];
        Dtreetrack->size[1] = Dtreetrack->size[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        loop_ub = Dtreetrack->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Dtreetrack->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            Dtreetrack->data[i5 + Dtreetrack->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          ixstop = 0;
          while (ixstop <= 0) {
            for (ic = 1; ic <= m; ic++) {
              Dtreetrack->data[ic - 1] = 0.0;
            }

            ixstop = m;
          }

          br = 0;
          ixstop = 0;
          while (ixstop <= 0) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (b_y->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= m; ic++) {
                  ia++;
                  Dtreetrack->data[ic] += b_y->data[ib] * Psi->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
            ixstop = m;
          }
        }
      }

      i4 = Arrow_inv->size[0] * Arrow_inv->size[1];
      Arrow_inv->size[0] = A_11inv->size[0] + Dtreetrack->size[0];
      Arrow_inv->size[1] = A_11inv->size[1] + Psi_notinv->size[1];
      emxEnsureCapacity((emxArray__common *)Arrow_inv, i4, (int32_T)sizeof
                        (real_T));
      loop_ub = A_11inv->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstop = A_11inv->size[0] - 1;
        for (i5 = 0; i5 <= ixstop; i5++) {
          Arrow_inv->data[i5 + Arrow_inv->size[0] * i4] = A_11inv->data[i5 +
            A_11inv->size[0] * i4] + C->data[i5 + C->size[0] * i4];
        }
      }

      loop_ub = Psi_notinv->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstop = Psi_notinv->size[0] - 1;
        for (i5 = 0; i5 <= ixstop; i5++) {
          Arrow_inv->data[i5 + Arrow_inv->size[0] * (i4 + A_11inv->size[1])] =
            Psi_notinv->data[i5 + Psi_notinv->size[0] * i4];
        }
      }

      loop_ub = Dtreetrack->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstop = Dtreetrack->size[0] - 1;
        for (i5 = 0; i5 <= ixstop; i5++) {
          Arrow_inv->data[(i5 + A_11inv->size[0]) + Arrow_inv->size[0] * i4] =
            Dtreetrack->data[i5 + Dtreetrack->size[0] * i4];
        }
      }

      loop_ub = Psi->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstop = Psi->size[0] - 1;
        for (i5 = 0; i5 <= ixstop; i5++) {
          Arrow_inv->data[(i5 + A_11inv->size[0]) + Arrow_inv->size[0] * (i4 +
            Dtreetrack->size[1])] = Psi->data[i5 + Psi->size[0] * i4];
        }
      }

      if (Binv_size[1] == 1) {
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = Binv_size[0];
        b_y->size[1] = u->size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = Binv_size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = u->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i4 + b_y->size[0] * i5] = 0.0;
            i6 = 0;
            while (i6 <= 0) {
              b_y->data[i4 + b_y->size[0] * i5] += Binv_data[i4] * u->data
                [u->size[0] * i5];
              i6 = 1;
            }
          }
        }
      } else {
        outsz[0] = (uint32_T)Binv_size[0];
        outsz[1] = (uint32_T)u->size[1];
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = (int32_T)outsz[0];
        b_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        i4 = b_y->size[0] * b_y->size[1];

        emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = b_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            b_y->data[i5 + b_y->size[0] * i4] = 0.0;
          }
        }

        if (Binv_size[0] == 0) {
        } else {
          pos = Binv_size[0] * (u->size[1] - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += Binv_size[0]) {
            i4 = ixstop + Binv_size[0];
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              b_y->data[ic] = 0.0;
            }
          }
        }
      }

      if ((b_y->size[1] == 1) || (Arrow_inv->size[0] == 1)) {
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = b_y->size[0];
        Dtreetrack->size[1] = Arrow_inv->size[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        loop_ub = b_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Arrow_inv->size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            Dtreetrack->data[i4 + Dtreetrack->size[0] * i5] = 0.0;
            pos = b_y->size[1] - 1;
            for (i6 = 0; i6 <= pos; i6++) {
              Dtreetrack->data[i4 + Dtreetrack->size[0] * i5] += b_y->data[i4 +
                b_y->size[0] * i6] * Arrow_inv->data[i6 + Arrow_inv->size[0] *
                i5];
            }
          }
        }
      } else {
        k = b_y->size[1];
        outsz[0] = (uint32_T)b_y->size[0];
        outsz[1] = (uint32_T)Arrow_inv->size[1];
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = (int32_T)outsz[0];
        Dtreetrack->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        m = b_y->size[0];
        ia = Arrow_inv->size[1];
        i4 = Dtreetrack->size[0] * Dtreetrack->size[1];
        Dtreetrack->size[0] = Dtreetrack->size[0];
        Dtreetrack->size[1] = Dtreetrack->size[1];
        emxEnsureCapacity((emxArray__common *)Dtreetrack, i4, (int32_T)sizeof
                          (real_T));
        loop_ub = Dtreetrack->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Dtreetrack->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            Dtreetrack->data[i5 + Dtreetrack->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (ia == 0)) {
        } else {
          pos = m * (ia - 1);
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            i4 = ixstop + m;
            for (ic = ixstop; ic + 1 <= i4; ic++) {
              Dtreetrack->data[ic] = 0.0;
            }
          }

          br = 0;
          for (ixstop = 0; ixstop <= pos; ixstop += m) {
            ar = 0;
            i4 = br + k;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (Arrow_inv->data[ib] != 0.0) {
                ia = ar;
                i5 = ixstop + m;
                for (ic = ixstop; ic + 1 <= i5; ic++) {
                  ia++;
                  Dtreetrack->data[ic] += Arrow_inv->data[ib] * b_y->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
          }
        }
      }

      i4 = a->size[0] * a->size[1];
      a->size[0] = u->size[1];
      a->size[1] = 10;
      emxEnsureCapacity((emxArray__common *)a, i4, (int32_T)sizeof(real_T));
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = u->size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          a->data[i5 + a->size[0] * i4] = u->data[i4 + u->size[0] * i5];
        }
      }

      if ((Dtreetrack->size[1] == 1) || (a->size[0] == 1)) {
        i4 = e_y->size[0] * e_y->size[1];
        e_y->size[0] = Dtreetrack->size[0];
        e_y->size[1] = 10;
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        loop_ub = Dtreetrack->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          for (i5 = 0; i5 < 10; i5++) {
            e_y->data[i4 + e_y->size[0] * i5] = 0.0;
            ixstop = Dtreetrack->size[1] - 1;
            for (i6 = 0; i6 <= ixstop; i6++) {
              e_y->data[i4 + e_y->size[0] * i5] += Dtreetrack->data[i4 +
                Dtreetrack->size[0] * i6] * a->data[i6 + a->size[0] * i5];
            }
          }
        }
      } else {
        k = Dtreetrack->size[1];
        outsz[0] = (uint32_T)Dtreetrack->size[0];
        outsz[1] = 10U;
        i4 = e_y->size[0] * e_y->size[1];
        e_y->size[0] = (int32_T)outsz[0];
        e_y->size[1] = 10;
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        m = Dtreetrack->size[0];
        i4 = e_y->size[0] * e_y->size[1];
        e_y->size[0] = e_y->size[0];
        e_y->size[1] = 10;
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        for (i4 = 0; i4 < 10; i4++) {
          loop_ub = e_y->size[0] - 1;
          for (i5 = 0; i5 <= loop_ub; i5++) {
            e_y->data[i5 + e_y->size[0] * i4] = 0.0;
          }
        }

        c_eml_xgemm(m, k, Dtreetrack, m, a, k, e_y, m);
      }

      if (Binv_size[0] == 1) {
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = e_y->size[0];
        c_y->size[1] = Binv_size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        loop_ub = e_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = Binv_size[1] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            c_y->data[i4 + c_y->size[0] * i5] = 0.0;
            for (i6 = 0; i6 < 10; i6++) {
              c_y->data[i4 + c_y->size[0] * i5] += e_y->data[i4 + e_y->size[0] *
                i6] * Binv_data[i6 + i5];
            }
          }
        }
      } else {
        outsz[0] = (uint32_T)e_y->size[0];
        outsz[1] = (uint32_T)Binv_size[1];
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = (int32_T)outsz[0];
        c_y->size[1] = (int32_T)outsz[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        m = e_y->size[0];
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = c_y->size[0];
        c_y->size[1] = c_y->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        loop_ub = c_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstop = c_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstop; i5++) {
            c_y->data[i5 + c_y->size[0] * i4] = 0.0;
          }
        }

        b_eml_xgemm(m, Binv_size[1], e_y, m, Binv_data, Binv_size, c_y, m);
      }

      loop_ub = Binv_size[0] * Binv_size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        epsilon = Binv_data[i4] - c_y->data[i4];
        Binv_data[i4] = epsilon;
      }

      for (i4 = 0; i4 < 10; i4++) {
        epsilon = 0.0;
        for (i5 = 0; i5 < 10; i5++) {
          epsilon += addaTree[i4 + 10 * i5] * dx_1[i5];
        }

        r[i4] = c[i4] - epsilon;
      }

      if (Binv_size[1] == 1) {
        z_size[0] = Binv_size[0];
        loop_ub = Binv_size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          z_data[i4] = 0.0;
          for (i5 = 0; i5 < 10; i5++) {
            z_data[i4] += Binv_data[i4 + Binv_size[0] * i5] * r[i5];
          }
        }
      } else {
        outsz[0] = (uint32_T)Binv_size[0];
        outsz[1] = 1U;
        z_size[0] = (int32_T)outsz[0];
        pos = z_size[0];
        loop_ub = pos - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          z_data[i4] = 0.0;
        }

        d_eml_xgemm(Binv_size[0], 0, Binv_data, Binv_size, Binv_size[0], r, 0,
                    z_data, z_size, Binv_size[0]);
      }

      m = z_size[0];
      loop_ub = z_size[0] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        p_data[i4] = z_data[i4];
      }

      /* condition_num_MinvADDA=cond(Minv*ADDA); */
      /* sprintf('Condition number for Minv*ADDA: %g', condition_num_MinvADDA) */
      /* %Step 4 */
      /* PCG */
      *PCGit = 1.0;
      i = 0;
      exitg3 = FALSE;
      while ((exitg3 == 0U) && (i < 10)) {
        *PCGit = 1.0 + (real_T)i;
        if (z_size[0] == 1) {
          f_y = 0.0;
          for (i4 = 0; i4 < 10; i4++) {
            f_y += r[i4] * z_data[i4];
          }
        } else {
          f_y = 0.0;
          for (k = 0; k < 10; k++) {
            f_y += r[k] * z_data[k];
          }
        }

        loop_ub = m - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          a_data[i4] = p_data[i4];
        }

        if (m == 1) {
          for (i4 = 0; i4 < 10; i4++) {
            g_y[i4] = 0.0;
            for (i5 = 0; i5 < 10; i5++) {
              g_y[i4] += a_data[i5] * addaTree[i5 + 10 * i4];
            }
          }
        } else {
          memset(&g_y[0], 0, 10U * sizeof(real_T));
          for (ixstop = 0; ixstop < 10; ixstop++) {
            for (ic = ixstop; ic + 1 <= ixstop + 1; ic++) {
              g_y[ic] = 0.0;
            }
          }

          br = 0;
          for (ixstop = 0; ixstop < 10; ixstop++) {
            ar = 0;
            i4 = br + m;
            for (ib = br; ib + 1 <= i4; ib++) {
              if (addaTree[ib] != 0.0) {
                ia = ar;
                for (ic = ixstop; ic + 1 <= ixstop + 1; ic++) {
                  ia++;
                  g_y[ic] += addaTree[ib] * a_data[ia - 1];
                }
              }

              ar++;
            }

            br += m;
          }
        }

        if (m == 1) {
          epsilon = 0.0;
          for (i4 = 0; i4 < 10; i4++) {
            epsilon += g_y[i4] * p_data[i4];
          }
        } else {
          epsilon = 0.0;
          for (k = 0; k < 10; k++) {
            epsilon += g_y[k] * p_data[k];
          }
        }

        epsilon = f_y / epsilon;

        /* dx_1_new2=dx_1+(alpha_1*p); */
        for (pos = 0; pos < 10; pos++) {
          dx_1_new[pos] = dx_1[pos] + epsilon * p_data[pos];
        }

        /*                          dx_1_new2=dx_1+alpha_1*p; */
        /*                          test=all(dx_1_new2==dx_1_new) */
        for (i4 = 0; i4 < 100; i4++) {
          addaTrack[i4] = epsilon * addaTree[i4];
        }

        if (m == 1) {
          for (i4 = 0; i4 < 10; i4++) {
            r_new[i4] = 0.0;
            for (i5 = 0; i5 < 10; i5++) {
              r_new[i4] += addaTrack[i4 + 10 * i5] * p_data[i5];
            }
          }
        } else {
          memset(&r_new[0], 0, 10U * sizeof(real_T));
          if (10 == m) {
            for (i4 = 0; i4 < 10; i4++) {
              r_new[i4] = 0.0;
              for (i5 = 0; i5 < 10; i5++) {
                r_new[i4] += addaTrack[i4 + 10 * i5] * p_data[i5];
              }
            }
          } else {
            memset(&r_new[0], 0, 10U * sizeof(real_T));
            ar = 0;
            for (ib = 0; ib < 10; ib++) {
              if (p_data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic < 10; ic++) {
                  ia++;
                  r_new[ic] += p_data[ib] * addaTrack[ia - 1];
                }
              }

              ar += 10;
            }
          }
        }

        f_y = 0.0;
        for (i4 = 0; i4 < 10; i4++) {
          epsilon = r[i4] - r_new[i4];
          f_y += epsilon * epsilon;
          r_new[i4] = epsilon;
        }

        if (sqrt(f_y) < 1.0E-5) {
          /* norm of rnew */
          /*                              sprintf('Iteration count for PCG: %d, with value', i, sqrt(r_new'*r_new)) */
          exitg3 = TRUE;
        } else {
          if (Binv_size[1] == 1) {
            z_new_size[0] = Binv_size[0];
            loop_ub = Binv_size[0] - 1;
            for (i4 = 0; i4 <= loop_ub; i4++) {
              z_new_data[i4] = 0.0;
              for (i5 = 0; i5 < 10; i5++) {
                z_new_data[i4] += Binv_data[i4 + Binv_size[0] * i5] * r_new[i5];
              }
            }
          } else {
            outsz[0] = (uint32_T)Binv_size[0];
            outsz[1] = 1U;
            z_new_size[0] = (int32_T)outsz[0];
            pos = z_new_size[0];
            loop_ub = pos - 1;
            for (i4 = 0; i4 <= loop_ub; i4++) {
              z_new_data[i4] = 0.0;
            }

            d_eml_xgemm(Binv_size[0], 0, Binv_data, Binv_size, Binv_size[0],
                        r_new, 0, z_new_data, z_new_size, Binv_size[0]);
          }

          loop_ub = z_new_size[0] - 1;
          for (i4 = 0; i4 <= loop_ub; i4++) {
            a_data[i4] = z_new_data[i4];
          }

          if (z_new_size[0] == 1) {
            f_y = 0.0;
            for (i4 = 0; i4 < 10; i4++) {
              f_y += a_data[i4] * r_new[i4];
            }
          } else {
            f_y = 0.0;
            if (z_new_size[0] < 1) {
            } else {
              for (k = 0; k + 1 <= z_new_size[0]; k++) {
                f_y += a_data[k] * r_new[k];
              }
            }
          }

          loop_ub = z_size[0] - 1;
          for (i4 = 0; i4 <= loop_ub; i4++) {
            a_data[i4] = z_data[i4];
          }

          if (z_size[0] == 1) {
            epsilon = 0.0;
            for (i4 = 0; i4 < 10; i4++) {
              epsilon += a_data[i4] * r[i4];
            }
          } else {
            epsilon = 0.0;
            if (z_size[0] < 1) {
            } else {
              for (k = 0; k + 1 <= z_size[0]; k++) {
                epsilon += a_data[k] * r[k];
              }
            }
          }

          epsilon = f_y / epsilon;
          m = z_new_size[0];
          loop_ub = z_new_size[0] - 1;
          for (i4 = 0; i4 <= loop_ub; i4++) {
            p_data[i4] = z_new_data[i4] + epsilon * p_data[i4];
          }

          for (pos = 0; pos < 10; pos++) {
            r[pos] = r_new[pos];
            dx_1[pos] = dx_1_new[pos];
          }

          z_size[0] = z_new_size[0];
          loop_ub = z_new_size[0] - 1;
          for (i4 = 0; i4 <= loop_ub; i4++) {
            z_data[i4] = z_new_data[i4];
          }

          i++;
        }
      }
    }

    /*         %% Step 5 Scaling */
    /*         %% Step 6 Calculating the interior point from the projection */
    for (i = 0; i < 17; i++) {
      e_D[i] = -w[i];
    }

    for (i4 = 0; i4 < 10; i4++) {
      for (i5 = 0; i5 < 17; i5++) {
        A_data[i5 + 17 * i4] = -A[i5 + 17 * i4];
      }
    }

    for (i4 = 0; i4 < 17; i4++) {
      c_A[i4] = 0.0;
      for (i5 = 0; i5 < 10; i5++) {
        c_A[i4] += A_data[i4 + 17 * i5] * dx_1[i5];
      }
    }

    rdivide(e_D, c_A, w);
    for (i = 0; i < 17; i++) {
      b_A[i] = (w[i] < 0.0);
    }

    eml_li_find(b_A, tmp_data, tmp_size);
    loop_ub = tmp_size[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      w[tmp_data[i4] - 1] = 1.0E+6;
    }

    ar = 1;
    alpha_2 = w[0];
    if (rtIsNaN(w[0])) {
      ib = 2;
      exitg2 = FALSE;
      while ((exitg2 == 0U) && (ib < 18)) {
        ar = ib;
        if (!rtIsNaN(w[ib - 1])) {
          alpha_2 = w[ib - 1];
          exitg2 = TRUE;
        } else {
          ib++;
        }
      }
    }

    if (ar < 17) {
      while (ar + 1 < 18) {
        if (w[ar] < alpha_2) {
          alpha_2 = w[ar];
        }

        ar++;
      }
    }

    alpha_2 *= b_gamma;

    /*         %% Step 1 */
    f_y = 0.0;
    epsilon = 0.0;
    h_y = 0.0;
    for (i = 0; i < 10; i++) {
      b_r = x[i] + alpha_2 * dx_1[i];
      f_y += c[i] * b_r;
      epsilon += c[i] * x[i];
      h_y += c[i] * x[i];
      r[i] = b_r;
    }

    guard1 = FALSE;
    if ((f_y - epsilon) / h_y > 0.01) {
      memcpy(&x[0], &r[0], 10U * sizeof(real_T));
      guard1 = TRUE;
    } else {
      f_y = 0.0;
      epsilon = 0.0;
      h_y = 0.0;
      for (k = 0; k < 10; k++) {
        f_y += c[k] * r[k];
        epsilon += c[k] * x[k];
        h_y += c[k] * x[k];
      }

      if ((f_y - epsilon) / h_y <= 0.01) {
        exitg1 = TRUE;
      } else {
        guard1 = TRUE;
      }
    }

    if (guard1 == TRUE) {
      b_check++;
    }
  }

  emxFree_real_T(&d_D);
  emxFree_real_T(&c_D);
  emxFree_real_T(&b_D);
  emxFree_real_T(&D_track);
  emxFree_real_T(&D);
  emxFree_real_T(&r6);
  emxFree_real_T(&c_Dtreetrack);
  emxFree_real_T(&r5);
  emxFree_real_T(&b_Dtreetrack);
  emxFree_real_T(&r4);
  emxFree_int32_T(&r3);
  emxFree_real_T(&r2);
  emxFree_real_T(&r1);
  emxFree_real_T(&e_y);
  emxFree_real_T(&maxval);
  emxFree_real_T(&d_y);
  emxFree_real_T(&c_y);
  emxFree_real_T(&a);
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
  emxFree_int32_T(&r0);
  emxFree_real_T(&C);
  emxFree_real_T(&Arrow_inv);
  emxFree_real_T(&Psi);
  emxFree_real_T(&Psi_notinv);
  emxFree_boolean_T(&A_12_Psi_to_remove);
  emxFree_real_T(&A_12_Psi);
  emxFree_real_T(&A_12);
  emxFree_real_T(&A_11inv);
  emxFree_real_T(&u);
  emxFree_real_T(&Dtreetrack);
  emxFree_real_T(&D_track_tempVec);
  *IPMit = (real_T)check;
  emxInit_real_T(&i_y, 2);
  if (Binv_size[1] == 1) {
    i0 = i_y->size[0] * i_y->size[1];
    i_y->size[0] = Binv_size[0];
    i_y->size[1] = 10;
    emxEnsureCapacity((emxArray__common *)i_y, i0, (int32_T)sizeof(real_T));
    loop_ub = Binv_size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      for (i1 = 0; i1 < 10; i1++) {
        i_y->data[i0 + i_y->size[0] * i1] = 0.0;
        i2 = 0;
        while (i2 <= 0) {
          i_y->data[i0 + i_y->size[0] * i1] += Binv_data[i0] * addaTree[10 * i1];
          i2 = 1;
        }
      }
    }
  } else {
    outsz[0] = (uint32_T)Binv_size[0];
    outsz[1] = 10U;
    i0 = i_y->size[0] * i_y->size[1];
    i_y->size[0] = (int32_T)outsz[0];
    i_y->size[1] = 10;
    emxEnsureCapacity((emxArray__common *)i_y, i0, (int32_T)sizeof(real_T));
    i0 = i_y->size[0] * i_y->size[1];
    i_y->size[0] = i_y->size[0];
    i_y->size[1] = 10;
    emxEnsureCapacity((emxArray__common *)i_y, i0, (int32_T)sizeof(real_T));
    for (i0 = 0; i0 < 10; i0++) {
      loop_ub = i_y->size[0] - 1;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        i_y->data[i1 + i_y->size[0] * i0] = 0.0;
      }
    }

    e_eml_xgemm(Binv_size[0], 0, Binv_data, Binv_size, Binv_size[0], addaTree, 0,
                i_y, Binv_size[0]);
  }

  *CondNum = cond(i_y);

  /* CondNum=condition_num_MinvADDA; */
  emxFree_real_T(&i_y);
}

/* End of code generation (Function_PCG_Wood_clean_codegen2.c) */
