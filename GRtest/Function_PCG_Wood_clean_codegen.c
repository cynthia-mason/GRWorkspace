/*
 * Function_PCG_Wood_clean_codegen.c
 *
 * Code generation for function 'Function_PCG_Wood_clean_codegen'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

#include "math.h"
/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen.h"
#include "Function_PCG_Wood_clean_codegen_emxutil.h"
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
static void b_eml_xgemm(int32_T m, int32_T n, int32_T k, const emxArray_real_T
  *A, int32_T lda, const real_T B_data[170], const int32_T B_size[2], int32_T
  ldb, emxArray_real_T *C, int32_T ldc);
static void c_eml_li_find(const boolean_T x_data[289], const int32_T x_size[2],
  int32_T y_data[289], int32_T y_size[1]);
static void c_eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
  lda, const real_T B_data[170], const int32_T B_size[2], int32_T ldb,
  emxArray_real_T *C, int32_T ldc);
static void d_eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
  lda, const real_T B[10], int32_T ldb, emxArray_real_T *C, int32_T ldc);
static void e_eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
  lda, const real_T B[100], int32_T ldb, emxArray_real_T *C, int32_T ldc);
static void eml_li_find(const boolean_T x[17], int32_T y_data[17], int32_T
  y_size[1]);
static void eml_xgemm(int32_T m, int32_T n, const emxArray_real_T *A, int32_T
                      lda, const emxArray_real_T *B, emxArray_real_T *C, int32_T
                      ldc);

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

static void b_eml_xgemm(int32_T m, int32_T n, int32_T k, const emxArray_real_T
  *A, int32_T lda, const real_T B_data[170], const int32_T B_size[2], int32_T
  ldb, emxArray_real_T *C, int32_T ldc)
{
  emxArray_real_T *b_C;
  int32_T c;
  int32_T cr;
  int32_T i12;
  int32_T ic;
  emxInit_real_T(&b_C, 2);
  if ((m == 0) || (n == 0)) {
  } else {
    c = ldc * (n - 1);
    cr = 0;
    while ((ldc > 0) && (cr <= c)) {
      i12 = cr + m;
      for (ic = cr; ic + 1 <= i12; ic++) {
        C->data[ic] = 0.0;
      }

      cr += ldc;
    }
  }

  emxFree_real_T(&b_C);
}

static void c_eml_li_find(const boolean_T x_data[289], const int32_T x_size[2],
  int32_T y_data[289], int32_T y_size[1])
{
  int32_T n;
  int32_T k;
  int32_T i;
  n = x_size[0] * x_size[1];
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x_data[i - 1]) {
      k++;
    }
  }

  y_size[0] = k;
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x_data[i - 1]) {
      y_data[k] = i;
      k++;
    }
  }
}

static void c_eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
  lda, const real_T B_data[170], const int32_T B_size[2], int32_T ldb,
  emxArray_real_T *C, int32_T ldc)
{
  emxArray_real_T *b_C;
  int32_T c;
  int32_T cr;
  int32_T i16;
  int32_T ic;
  int32_T br;
  int32_T ar;
  int32_T ib;
  int32_T ia;
  int32_T i17;
  emxInit_real_T(&b_C, 2);
  if (m == 0) {
  } else {
    c = ldc * 9;
    cr = 0;
    while ((ldc > 0) && (cr <= c)) {
      i16 = cr + m;
      for (ic = cr; ic + 1 <= i16; ic++) {
        C->data[ic] = 0.0;
      }

      cr += ldc;
    }

    br = 0;
    cr = 0;
    while ((ldc > 0) && (cr <= c)) {
      ar = -1;
      i16 = br + k;
      for (ib = br; ib + 1 <= i16; ib++) {
        if (B_data[ib] != 0.0) {
          ia = ar;
          i17 = cr + m;
          for (ic = cr; ic + 1 <= i17; ic++) {
            ia++;
            C->data[ic] += B_data[ib] * A->data[ia];
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

static void d_eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
  lda, const real_T B[10], int32_T ldb, emxArray_real_T *C, int32_T ldc)
{
  emxArray_real_T *b_C;
  int32_T cr;
  b_emxInit_real_T(&b_C, 1);
  if (m == 0) {
  } else {
    cr = 0;
    while ((ldc > 0) && (cr <= 0)) {
      for (cr = 1; cr <= m; cr++) {
        C->data[cr - 1] = 0.0;
      }

      cr = ldc;
    }
  }

  emxFree_real_T(&b_C);
}

static void e_eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
  lda, const real_T B[100], int32_T ldb, emxArray_real_T *C, int32_T ldc)
{
  int32_T c;
  int32_T cr;
  int32_T i18;
  int32_T ic;
  if (m == 0) {
  } else {
    c = ldc * 9;
    cr = 0;
    while ((ldc > 0) && (cr <= c)) {
      i18 = cr + m;
      for (ic = cr; ic + 1 <= i18; ic++) {
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

static void eml_xgemm(int32_T m, int32_T n, const emxArray_real_T *A, int32_T
                      lda, const emxArray_real_T *B, emxArray_real_T *C, int32_T
                      ldc)
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
        if (B->data[ib] != 0.0) {
          ia = ar;
          for (ic = 0; ic + 1 <= m; ic++) {
            ia++;
            C->data[ic] += B->data[ib] * A->data[ia];
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

void Function_PCG_Wood_clean_codegen(const real_T A[170], const real_T b[17],
  const real_T c[10], real_T b_gamma, real_T x[10], real_T *NumVar, real_T
  *IPMit, real_T *CondNum, real_T *PCGit)
{
  emxArray_real_T *Binv;
  int32_T i0;
  real_T dx_1[10];
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
  emxArray_real_T *A_11inv;
  emxArray_real_T *Arrow_inv;
  emxArray_real_T *z;
  emxArray_real_T *p;
  emxArray_real_T *z_new;
  emxArray_real_T *r0;
  emxArray_real_T *b_b;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  emxArray_real_T *c_y;
  emxArray_real_T *a;
  emxArray_real_T *d_y;
  emxArray_real_T *c_b;
  emxArray_real_T *e_y;
  emxArray_real_T *f_y;
  emxArray_real_T *b_a;
  emxArray_real_T *r1;
  emxArray_real_T *r2;
  emxArray_real_T *D_tree;
  emxArray_real_T *r3;
  emxArray_real_T *b_D_tree;
  emxArray_real_T *r4;
  emxArray_real_T *r5;
  emxArray_real_T *D_track;
  boolean_T exitg1;
  real_T w[17];
  int32_T i4;
  int32_T i5;
  real_T D[17];
  real_T D_data[17];
  int32_T D_size[1];
  int32_T D_tree_size[2];
  real_T D_tree_data[289];
  int32_T b_D_size[1];
  int32_T b_tmp_size[2];
  real_T b_tmp_data[289];
  int32_T D_track_size[2];
  real_T D_track_data[256];
  int32_T c_D_size[1];
  int32_T D_nonneg_size[2];
  real_T D_nonneg_data[289];
  real_T a_data[2730];
  int32_T u_size[2];
  real_T u_data[170];
  int32_T ixstart;
  int32_T i6;
  int32_T outsz[2];
  int32_T n;
  int32_T cr;
  int32_T ic;
  int32_T ar;
  int32_T ib;
  int32_T ia;
  int32_T iy;
  real_T y_data[170];
  real_T A_data[170];
  real_T ADDA[100];
  int32_T ixstop;
  real_T g_y[100];
  real_T h_y[100];
  emxArray_real_T b_D_track_data;
  int32_T varargin_1_size[2];
  real_T varargin_1_data[2176];
  real_T alpha_2;
  boolean_T exitg8;
  real_T Breakpoint_data[1];
  real_T diag_ind_data[256];
  int16_T sz[2];
  real_T tracks_i_data[256];
  real_T i_y;
  boolean_T b_tracks_i_data[256];
  int32_T tracks_i_size[1];
  int32_T c_tmp_data[256];
  emxArray_real_T b_D_tree_data;
  real_T d_tmp_data[145];
  real_T b_A_data[160];
  int32_T a_size[2];
  real_T b_a_data[170];
  int32_T m;
  int32_T A_12_size[2];
  real_T A_12_data[289];
  int32_T i7;
  boolean_T exitg7;
  boolean_T exitg6;
  boolean_T exitg5;
  boolean_T exitg4;
  int32_T A_12_Psi_to_remove_size[2];
  boolean_T A_12_Psi_to_remove_data[289];
  int32_T e_tmp_data[289];
  int32_T f_tmp_data[289];
  real_T b_D_nonneg_data[289];
  int32_T b_D_nonneg_size[2];
  real_T r[10];
  boolean_T exitg3;
  real_T j_y[10];
  real_T dx_1_new[10];
  real_T r_new[10];
  boolean_T exitg2;
  real_T k_y;
  real_T b_r;
  boolean_T guard1 = FALSE;
  emxInit_real_T(&Binv, 2);

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
  i0 = Binv->size[0] * Binv->size[1];
  Binv->size[0] = 1;
  Binv->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)Binv, i0, (int32_T)sizeof(real_T));
  Binv->data[0] = 0.0;
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
  emxInit_real_T(&A_11inv, 2);
  emxInit_real_T(&Arrow_inv, 2);
  b_emxInit_real_T(&z, 1);
  b_emxInit_real_T(&p, 1);
  b_emxInit_real_T(&z_new, 1);
  emxInit_real_T(&r0, 2);
  emxInit_real_T(&b_b, 2);
  emxInit_real_T(&y, 2);
  emxInit_real_T(&b_y, 2);
  emxInit_real_T(&c_y, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&d_y, 2);
  emxInit_real_T(&c_b, 2);
  emxInit_real_T(&e_y, 2);
  emxInit_real_T(&f_y, 2);
  emxInit_real_T(&b_a, 2);
  emxInit_real_T(&r1, 2);
  emxInit_real_T(&r2, 2);
  emxInit_real_T(&D_tree, 2);
  emxInit_real_T(&r3, 2);
  emxInit_real_T(&b_D_tree, 2);
  emxInit_real_T(&r4, 2);
  emxInit_real_T(&r5, 2);
  emxInit_real_T(&D_track, 2);
  exitg1 = FALSE;
  while ((exitg1 == 0U) && (b_check < 20)) {
    check = b_check + 1;

    /* sprintf('Interior point iteration: %d', check) */
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
      D[i] = 1.0 / w[i];
    }

    /* inverse of slack, therefore the following D's are ALL inverses */
    D_size[0] = (int32_T)row_v_data[0];
    loop_ub = (int32_T)row_v_data[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      D_data[i4] = D[i4];
    }

    diag(D_data, D_size, D_tree_data, D_tree_size);
    if (row_v_data[0] + 1.0 > row_nonneg1_data[0] - 1.0) {
      i4 = 0;
      i5 = -1;
    } else {
      i4 = (int32_T)row_v_data[0];
      i5 = (int32_T)row_nonneg1_data[0] - 2;
    }

    b_D_size[0] = (i5 - i4) + 1;
    loop_ub = i5 - i4;
    for (i5 = 0; i5 <= loop_ub; i5++) {
      D_data[i5] = D[i4 + i5];
    }

    diag(D_data, b_D_size, b_tmp_data, b_tmp_size);
    D_track_size[0] = b_tmp_size[0];
    D_track_size[1] = b_tmp_size[1];
    loop_ub = b_tmp_size[0] * b_tmp_size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      D_track_data[i4] = b_tmp_data[i4];
    }

    if (row_nonneg1_data[0] > row_nonneg2_data[0]) {
      i4 = 0;
      i5 = 0;
    } else {
      i4 = (int32_T)row_nonneg1_data[0] - 1;
      i5 = (int32_T)row_nonneg2_data[0];
    }

    c_D_size[0] = i5 - i4;
    loop_ub = (i5 - i4) - 1;
    for (i5 = 0; i5 <= loop_ub; i5++) {
      D_data[i5] = D[i4 + i5];
    }

    diag(D_data, c_D_size, D_nonneg_data, D_nonneg_size);
    loop_ub = (int32_T)row_v_data[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      for (i5 = 0; i5 < 10; i5++) {
        a_data[i5 + 10 * i4] = A[i4 + 17 * i5];
      }
    }

    if (((int32_T)row_v_data[0] == 1) || (D_tree_size[0] == 1)) {
      u_size[0] = 10;
      u_size[1] = D_tree_size[1];
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = D_tree_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          u_data[i4 + 10 * i5] = 0.0;
          ixstart = (int32_T)row_v_data[0] - 1;
          for (i6 = 0; i6 <= ixstart; i6++) {
            u_data[i4 + 10 * i5] += a_data[i4 + 10 * i6] * D_tree_data[i6 +
              D_tree_size[0] * i5];
          }
        }
      }
    } else {
      outsz[0] = 10;
      outsz[1] = D_tree_size[1];
      u_size[0] = 10;
      u_size[1] = outsz[1];
      n = D_tree_size[1];
      u_size[0] = 10;
      loop_ub = u_size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        memset(&u_data[10 * i4], 0, 10U * sizeof(real_T));
      }

      if (n == 0) {
      } else {
        ixstart = 10 * (n - 1);
        for (cr = 1; cr - 1 <= ixstart; cr += 10) {
          for (ic = cr; ic <= cr + 9; ic++) {
            u_data[ic - 1] = 0.0;
          }
        }

        n = 0;
        for (cr = 0; cr <= ixstart; cr += 10) {
          ar = 0;
          i4 = n + (int32_T)row_v_data[0];
          for (ib = n; ib + 1 <= i4; ib++) {
            if (D_tree_data[ib] != 0.0) {
              ia = ar;
              for (ic = cr; ic + 1 <= cr + 10; ic++) {
                ia++;
                u_data[ic] += D_tree_data[ib] * a_data[ia - 1];
              }
            }

            ar += 10;
          }

          n += (int32_T)row_v_data[0];
        }
      }
    }

    if ((u_size[1] == 1) || (D_tree_size[0] == 1)) {
      iy = D_tree_size[1];
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = D_tree_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          y_data[i4 + 10 * i5] = 0.0;
          ixstart = u_size[1] - 1;
          for (i6 = 0; i6 <= ixstart; i6++) {
            y_data[i4 + 10 * i5] += u_data[i4 + 10 * i6] * D_tree_data[i6 +
              D_tree_size[0] * i5];
          }
        }
      }
    } else {
      outsz[0] = 10;
      outsz[1] = D_tree_size[1];
      iy = outsz[1];
      n = D_tree_size[1];
      loop_ub = outsz[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        memset(&y_data[10 * i4], 0, 10U * sizeof(real_T));
      }

      if (n == 0) {
      } else {
        ixstart = 10 * (n - 1);
        for (cr = 1; cr - 1 <= ixstart; cr += 10) {
          for (ic = cr; ic <= cr + 9; ic++) {
            y_data[ic - 1] = 0.0;
          }
        }

        n = 0;
        for (cr = 0; cr <= ixstart; cr += 10) {
          ar = 0;
          i4 = n + u_size[1];
          for (ib = n; ib + 1 <= i4; ib++) {
            if (D_tree_data[ib] != 0.0) {
              ia = ar;
              for (ic = cr; ic + 1 <= cr + 10; ic++) {
                ia++;
                y_data[ic] += D_tree_data[ib] * u_data[ia - 1];
              }
            }

            ar += 10;
          }

          n += u_size[1];
        }
      }
    }

    if ((iy == 1) || ((int32_T)row_v_data[0] == 1)) {
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = (int32_T)row_v_data[0] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          A_data[i5 + (int32_T)row_v_data[0] * i4] = A[i5 + 17 * i4];
        }
      }

      for (i4 = 0; i4 < 10; i4++) {
        for (i5 = 0; i5 < 10; i5++) {
          ADDA[i4 + 10 * i5] = 0.0;
          loop_ub = iy - 1;
          for (i6 = 0; i6 <= loop_ub; i6++) {
            ADDA[i4 + 10 * i5] += y_data[i4 + 10 * i6] * A_data[i6 + (int32_T)
              row_v_data[0] * i5];
          }
        }
      }
    } else {
      memset(&ADDA[0], 0, 100U * sizeof(real_T));
      for (cr = 0; cr < 92; cr += 10) {
        for (ic = cr; ic + 1 <= cr + 10; ic++) {
          ADDA[ic] = 0.0;
        }
      }

      n = 0;
      for (cr = 0; cr < 92; cr += 10) {
        ar = 0;
        i4 = n + iy;
        for (ib = n; ib + 1 <= i4; ib++) {
          if (A[ib % (int32_T)row_v_data[0] + 17 * (ib / (int32_T)row_v_data[0])]
              != 0.0) {
            ia = ar;
            for (ic = cr; ic + 1 <= cr + 10; ic++) {
              ia++;
              ADDA[ic] += A[ib % (int32_T)row_v_data[0] + 17 * (ib / (int32_T)
                row_v_data[0])] * y_data[ia - 1];
            }
          }

          ar += 10;
        }

        n += iy;
      }
    }

    ixstop = (i3 - i2) + 1;
    loop_ub = i3 - i2;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      for (i5 = 0; i5 < 10; i5++) {
        a_data[i5 + 10 * i4] = A[(i2 + i4) + 17 * i5];
      }
    }

    if ((ixstop == 1) || (D_track_size[0] == 1)) {
      u_size[0] = 10;
      u_size[1] = D_track_size[1];
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = D_track_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          u_data[i4 + 10 * i5] = 0.0;
          ixstart = ixstop - 1;
          for (i6 = 0; i6 <= ixstart; i6++) {
            u_data[i4 + 10 * i5] += a_data[i4 + 10 * i6] * D_track_data[i6 +
              D_track_size[0] * i5];
          }
        }
      }
    } else {
      outsz[0] = 10;
      outsz[1] = D_track_size[1];
      u_size[0] = 10;
      u_size[1] = outsz[1];
      n = D_track_size[1];
      u_size[0] = 10;
      loop_ub = u_size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        memset(&u_data[10 * i4], 0, 10U * sizeof(real_T));
      }

      if (n == 0) {
      } else {
        ixstart = 10 * (n - 1);
        for (cr = 1; cr - 1 <= ixstart; cr += 10) {
          for (ic = cr; ic <= cr + 9; ic++) {
            u_data[ic - 1] = 0.0;
          }
        }

        n = 0;
        for (cr = 0; cr <= ixstart; cr += 10) {
          ar = 0;
          i4 = n + ixstop;
          for (ib = n; ib + 1 <= i4; ib++) {
            if (D_track_data[ib] != 0.0) {
              ia = ar;
              for (ic = cr; ic + 1 <= cr + 10; ic++) {
                ia++;
                u_data[ic] += D_track_data[ib] * a_data[ia - 1];
              }
            }

            ar += 10;
          }

          n += ixstop;
        }
      }
    }

    if ((u_size[1] == 1) || (D_track_size[0] == 1)) {
      iy = D_track_size[1];
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = D_track_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          y_data[i4 + 10 * i5] = 0.0;
          ixstart = u_size[1] - 1;
          for (i6 = 0; i6 <= ixstart; i6++) {
            y_data[i4 + 10 * i5] += u_data[i4 + 10 * i6] * D_track_data[i6 +
              D_track_size[0] * i5];
          }
        }
      }
    } else {
      outsz[0] = 10;
      outsz[1] = D_track_size[1];
      iy = outsz[1];
      n = D_track_size[1];
      loop_ub = outsz[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        memset(&y_data[10 * i4], 0, 10U * sizeof(real_T));
      }

      if (n == 0) {
      } else {
        ixstart = 10 * (n - 1);
        for (cr = 1; cr - 1 <= ixstart; cr += 10) {
          for (ic = cr; ic <= cr + 9; ic++) {
            y_data[ic - 1] = 0.0;
          }
        }

        n = 0;
        for (cr = 0; cr <= ixstart; cr += 10) {
          ar = 0;
          i4 = n + u_size[1];
          for (ib = n; ib + 1 <= i4; ib++) {
            if (D_track_data[ib] != 0.0) {
              ia = ar;
              for (ic = cr; ic + 1 <= cr + 10; ic++) {
                ia++;
                y_data[ic] += D_track_data[ib] * u_data[ia - 1];
              }
            }

            ar += 10;
          }

          n += u_size[1];
        }
      }
    }

    if ((iy == 1) || ((i3 - i2) + 1 == 1)) {
      ixstart = (i3 - i2) + 1;
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = i3 - i2;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          A_data[i5 + ixstart * i4] = A[(i2 + i5) + 17 * i4];
        }
      }

      for (i4 = 0; i4 < 10; i4++) {
        for (i5 = 0; i5 < 10; i5++) {
          g_y[i4 + 10 * i5] = 0.0;
          loop_ub = iy - 1;
          for (i6 = 0; i6 <= loop_ub; i6++) {
            g_y[i4 + 10 * i5] += y_data[i4 + 10 * i6] * A_data[i6 + ixstart * i5];
          }
        }
      }
    } else {
      memset(&g_y[0], 0, 100U * sizeof(real_T));
      for (cr = 0; cr < 92; cr += 10) {
        for (ic = cr; ic + 1 <= cr + 10; ic++) {
          g_y[ic] = 0.0;
        }
      }

      n = 0;
      for (cr = 0; cr < 92; cr += 10) {
        ar = 0;
        i4 = n + iy;
        for (ib = n; ib + 1 <= i4; ib++) {
          ixstart = (i3 - i2) + 1;
          for (i5 = 0; i5 < 10; i5++) {
            loop_ub = i3 - i2;
            for (i6 = 0; i6 <= loop_ub; i6++) {
              A_data[i6 + ixstart * i5] = A[(i2 + i6) + 17 * i5];
            }
          }

          if (A_data[ib] != 0.0) {
            ia = ar;
            for (ic = cr; ic + 1 <= cr + 10; ic++) {
              ia++;
              ixstart = (i3 - i2) + 1;
              for (i5 = 0; i5 < 10; i5++) {
                loop_ub = i3 - i2;
                for (i6 = 0; i6 <= loop_ub; i6++) {
                  A_data[i6 + ixstart * i5] = A[(i2 + i6) + 17 * i5];
                }
              }

              g_y[ic] += A_data[ib] * y_data[ia - 1];
            }
          }

          ar += 10;
        }

        n += iy;
      }
    }

    ixstop = i1 - i0;
    loop_ub = (i1 - i0) - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      for (i5 = 0; i5 < 10; i5++) {
        a_data[i5 + 10 * i4] = A[(i0 + i4) + 17 * i5];
      }
    }

    if ((ixstop == 1) || (D_nonneg_size[0] == 1)) {
      u_size[0] = 10;
      u_size[1] = D_nonneg_size[1];
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = D_nonneg_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          u_data[i4 + 10 * i5] = 0.0;
          ixstart = ixstop - 1;
          for (i6 = 0; i6 <= ixstart; i6++) {
            u_data[i4 + 10 * i5] += a_data[i4 + 10 * i6] * D_nonneg_data[i6 +
              D_nonneg_size[0] * i5];
          }
        }
      }
    } else {
      outsz[0] = 10;
      outsz[1] = D_nonneg_size[1];
      u_size[0] = 10;
      u_size[1] = outsz[1];
      u_size[0] = 10;
      loop_ub = u_size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        memset(&u_data[10 * i4], 0, 10U * sizeof(real_T));
      }

      if (D_nonneg_size[1] == 0) {
      } else {
        ixstart = 10 * (D_nonneg_size[1] - 1);
        for (cr = 1; cr - 1 <= ixstart; cr += 10) {
          for (ic = cr; ic <= cr + 9; ic++) {
            u_data[ic - 1] = 0.0;
          }
        }

        n = 0;
        for (cr = 0; cr <= ixstart; cr += 10) {
          ar = 0;
          i4 = n + ixstop;
          for (ib = n; ib + 1 <= i4; ib++) {
            if (D_nonneg_data[ib] != 0.0) {
              ia = ar;
              for (ic = cr; ic + 1 <= cr + 10; ic++) {
                ia++;
                u_data[ic] += D_nonneg_data[ib] * a_data[ia - 1];
              }
            }

            ar += 10;
          }

          n += ixstop;
        }
      }
    }

    if ((u_size[1] == 1) || (D_nonneg_size[0] == 1)) {
      iy = D_nonneg_size[1];
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = D_nonneg_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          y_data[i4 + 10 * i5] = 0.0;
          ixstart = u_size[1] - 1;
          for (i6 = 0; i6 <= ixstart; i6++) {
            y_data[i4 + 10 * i5] += u_data[i4 + 10 * i6] * D_nonneg_data[i6 +
              D_nonneg_size[0] * i5];
          }
        }
      }
    } else {
      outsz[0] = 10;
      outsz[1] = D_nonneg_size[1];
      iy = outsz[1];
      loop_ub = outsz[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        memset(&y_data[10 * i4], 0, 10U * sizeof(real_T));
      }

      if (D_nonneg_size[1] == 0) {
      } else {
        ixstart = 10 * (D_nonneg_size[1] - 1);
        for (cr = 1; cr - 1 <= ixstart; cr += 10) {
          for (ic = cr; ic <= cr + 9; ic++) {
            y_data[ic - 1] = 0.0;
          }
        }

        n = 0;
        for (cr = 0; cr <= ixstart; cr += 10) {
          ar = 0;
          i4 = n + u_size[1];
          for (ib = n; ib + 1 <= i4; ib++) {
            if (D_nonneg_data[ib] != 0.0) {
              ia = ar;
              for (ic = cr; ic + 1 <= cr + 10; ic++) {
                ia++;
                y_data[ic] += D_nonneg_data[ib] * u_data[ia - 1];
              }
            }

            ar += 10;
          }

          n += u_size[1];
        }
      }
    }

    if ((iy == 1) || (i1 - i0 == 1)) {
      ixstart = i1 - i0;
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = (i1 - i0) - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          A_data[i5 + ixstart * i4] = A[(i0 + i5) + 17 * i4];
        }
      }

      for (i4 = 0; i4 < 10; i4++) {
        for (i5 = 0; i5 < 10; i5++) {
          h_y[i4 + 10 * i5] = 0.0;
          loop_ub = iy - 1;
          for (i6 = 0; i6 <= loop_ub; i6++) {
            h_y[i4 + 10 * i5] += y_data[i4 + 10 * i6] * A_data[i6 + ixstart * i5];
          }
        }
      }
    } else {
      memset(&h_y[0], 0, 100U * sizeof(real_T));
      for (cr = 0; cr < 92; cr += 10) {
        for (ic = cr; ic + 1 <= cr + 10; ic++) {
          h_y[ic] = 0.0;
        }
      }

      n = 0;
      for (cr = 0; cr < 92; cr += 10) {
        ar = 0;
        i4 = n + iy;
        for (ib = n; ib + 1 <= i4; ib++) {
          if (A[(i0 + ib % (i1 - i0)) + 17 * (ib / (i1 - i0))] != 0.0) {
            ia = ar;
            for (ic = cr; ic + 1 <= cr + 10; ic++) {
              ia++;
              h_y[ic] += A[(i0 + ib % (i1 - i0)) + 17 * (ib / (i1 - i0))] *
                y_data[ia - 1];
            }
          }

          ar += 10;
        }

        n += iy;
      }
    }

    for (i4 = 0; i4 < 100; i4++) {
      ADDA[i4] = (ADDA[i4] + g_y[i4]) + h_y[i4];
    }

    b_D_track_data.data = (real_T *)&D_track_data;
    b_D_track_data.size = (int32_T *)&D_track_size;
    b_D_track_data.allocatedSize = 256;
    b_D_track_data.numDimensions = 2;
    b_D_track_data.canFreeData = FALSE;
    b_diag(&b_D_track_data, r1);
    varargin_1_size[0] = r1->size[0];
    varargin_1_size[1] = r1->size[1];
    loop_ub = r1->size[0] * r1->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      varargin_1_data[i4] = r1->data[i4];
    }

    for (i4 = 0; i4 < 2; i4++) {
      outsz[i4] = varargin_1_size[i4];
    }

    outsz[0] = 1;
    n = 0;
    i = 1;
    while (i <= varargin_1_size[1]) {
      ixstart = n;
      ixstop = n + varargin_1_size[0];
      alpha_2 = varargin_1_data[n];
      if (varargin_1_size[0] > 1) {
        if (rtIsNaN(varargin_1_data[n])) {
          cr = n + 1;
          exitg8 = FALSE;
          while ((exitg8 == 0U) && (cr + 1 <= ixstop)) {
            ixstart = cr;
            if (!rtIsNaN(varargin_1_data[cr])) {
              alpha_2 = varargin_1_data[cr];
              exitg8 = TRUE;
            } else {
              cr++;
            }
          }
        }

        if (ixstart + 1 < ixstop) {
          for (cr = ixstart + 1; cr + 1 <= ixstop; cr++) {
            if (varargin_1_data[cr] > alpha_2) {
              alpha_2 = varargin_1_data[cr];
            }
          }
        }
      }

      D[0] = alpha_2;
      n += varargin_1_size[0];
      i = 2;
    }

    loop_ub = outsz[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      Breakpoint_data[i4] = D[i4] / 10.0;
    }

    n = D_track_size[0] * D_track_size[0];
    loop_ub = n - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      diag_ind_data[i4] = 0.0;
    }

    n = 0;
    *PCGit = 1.0;
    while (n + 1 <= D_track_size[0] * D_track_size[0] + 1) {
      diag_ind_data[n] = *PCGit;
      n = (n + D_track_size[0]) + 1;
      (*PCGit)++;
    }

    /* bool_tracks2=D_track>Breakpoint; */
    n = D_track_size[0] * D_track_size[1];
    for (i4 = 0; i4 < 2; i4++) {
      sz[i4] = 0;
    }

    sz[0] = (int16_T)(D_track_size[0] * D_track_size[1]);
    sz[1] = 1;
    for (iy = 0; iy + 1 <= n; iy++) {
      tracks_i_data[iy] = D_track_data[iy];
    }

    n = D_track_size[0] * D_track_size[1];
    if ((outsz[1] == 1) || (n == 1)) {
      i_y = 0.0;
      loop_ub = outsz[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        i_y += Breakpoint_data[i4];
      }
    } else {
      i_y = 0.0;
    }

    /* tracks_i2=diag_ind(bool_tracks2);%provides row and column(omitted here) indicies of values greater than value (each in vector) */
    tracks_i_size[0] = sz[0];
    loop_ub = sz[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      b_tracks_i_data[i4] = (tracks_i_data[i4] > i_y);
    }

    b_eml_li_find(b_tracks_i_data, tracks_i_size, c_tmp_data, tmp_size);
    loop_ub = tmp_size[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      tracks_i_data[i4] = diag_ind_data[c_tmp_data[i4] - 1];
    }

    /* [tracks_i,~]=find(D_track>Breakpoint); */
    /* test=all(tracks_i2==tracks_i) */
    /* adds the active tracks */
    b_D_tree_data.data = (real_T *)&D_tree_data;
    b_D_tree_data.size = (int32_T *)&D_tree_size;
    b_D_tree_data.allocatedSize = 289;
    b_D_tree_data.numDimensions = 2;
    b_D_tree_data.canFreeData = FALSE;
    b_diag(&b_D_tree_data, r1);
    n = r1->size[0];
    iy = r1->size[1];
    loop_ub = r1->size[0] * r1->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      d_tmp_data[i4] = r1->data[i4];
    }

    i4 = D_track->size[0] * D_track->size[1];
    D_track->size[0] = tmp_size[0];
    D_track->size[1] = tmp_size[0];
    emxEnsureCapacity((emxArray__common *)D_track, i4, (int32_T)sizeof(real_T));
    loop_ub = tmp_size[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      ixstart = tmp_size[0] - 1;
      for (i5 = 0; i5 <= ixstart; i5++) {
        D_track->data[i5 + D_track->size[0] * i4] = D_track_data[((int32_T)
          tracks_i_data[i5] + D_track_size[0] * ((int32_T)tracks_i_data[i4] - 1))
          - 1];
      }
    }

    b_diag(D_track, r1);
    i4 = r0->size[0] * r0->size[1];
    r0->size[0] = r1->size[0];
    r0->size[1] = r1->size[1];
    emxEnsureCapacity((emxArray__common *)r0, i4, (int32_T)sizeof(real_T));
    loop_ub = r1->size[0] * r1->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      r0->data[i4] = r1->data[i4];
    }

    ixstart = (i3 - i2) + 1;
    for (i4 = 0; i4 < 10; i4++) {
      loop_ub = i3 - i2;
      for (i5 = 0; i5 <= loop_ub; i5++) {
        b_A_data[i5 + ixstart * i4] = A[(i2 + i5) + 17 * i4];
      }
    }

    ixstop = (int32_T)row_v_data[0] + tmp_size[0];
    loop_ub = (int32_T)row_v_data[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      for (i5 = 0; i5 < 10; i5++) {
        a_data[i5 + 10 * i4] = A[i4 + 17 * i5];
      }
    }

    loop_ub = tmp_size[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      for (i5 = 0; i5 < 10; i5++) {
        a_data[i5 + 10 * (i4 + (int32_T)row_v_data[0])] = b_A_data[((int32_T)
          tracks_i_data[i4] + ixstart * i5) - 1];
      }
    }

    i4 = r5->size[0] * r5->size[1];
    r5->size[0] = n + r0->size[0];
    r5->size[1] = iy;
    emxEnsureCapacity((emxArray__common *)r5, i4, (int32_T)sizeof(real_T));
    loop_ub = iy - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      ixstart = n - 1;
      for (i5 = 0; i5 <= ixstart; i5++) {
        r5->data[i5 + r5->size[0] * i4] = d_tmp_data[i5 + n * i4];
      }
    }

    loop_ub = r0->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      ixstart = r0->size[0] - 1;
      for (i5 = 0; i5 <= ixstart; i5++) {
        r5->data[(i5 + n) + r5->size[0] * i4] = r0->data[i5 + r0->size[0] * i4];
      }
    }

    b_diag(r5, r1);
    i4 = b_b->size[0] * b_b->size[1];
    b_b->size[0] = r1->size[0];
    b_b->size[1] = r1->size[1];
    emxEnsureCapacity((emxArray__common *)b_b, i4, (int32_T)sizeof(real_T));
    loop_ub = r1->size[0] * r1->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      b_b->data[i4] = r1->data[i4];
    }

    if ((ixstop == 1) || (b_b->size[0] == 1)) {
      u_size[0] = 10;
      u_size[1] = b_b->size[1];
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = b_b->size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          u_data[i4 + 10 * i5] = 0.0;
          ixstart = ixstop - 1;
          for (i6 = 0; i6 <= ixstart; i6++) {
            u_data[i4 + 10 * i5] += a_data[i4 + 10 * i6] * b_b->data[i6 +
              b_b->size[0] * i5];
          }
        }
      }
    } else {
      outsz[0] = 10;
      outsz[1] = b_b->size[1];
      u_size[0] = 10;
      u_size[1] = outsz[1];
      n = b_b->size[1];
      u_size[0] = 10;
      loop_ub = u_size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        memset(&u_data[10 * i4], 0, 10U * sizeof(real_T));
      }

      if (n == 0) {
      } else {
        for (ic = 1; ic < 11; ic++) {
          u_data[ic - 1] = 0.0;
        }

        ar = 0;
        i4 = ixstop;
        for (ib = 0; ib + 1 <= i4; ib++) {
          if (b_b->data[ib] != 0.0) {
            ia = ar;
            for (ic = 0; ic + 1 < 11; ic++) {
              ia++;
              u_data[ic] += b_b->data[ib] * a_data[ia - 1];
            }
          }

          ar += 10;
        }
      }
    }

    /* Prep for step 4: */
    if (!(tmp_size[0] == 0)) {
      if ((D_nonneg_size[1] == 1) || (D_nonneg_size[0] == 1)) {
        i4 = y->size[0] * y->size[1];
        y->size[0] = D_nonneg_size[0];
        y->size[1] = D_nonneg_size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = D_nonneg_size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = D_nonneg_size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i4 + y->size[0] * i5] = 0.0;
            n = D_nonneg_size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              y->data[i4 + y->size[0] * i5] += D_nonneg_data[i4 + D_nonneg_size
                [0] * i6] * D_nonneg_data[i6 + D_nonneg_size[0] * i5];
            }
          }
        }
      } else {
        outsz[0] = D_nonneg_size[0];
        outsz[1] = D_nonneg_size[1];
        i4 = y->size[0] * y->size[1];
        y->size[0] = outsz[0];
        y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        i4 = y->size[0] * y->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i5 + y->size[0] * i4] = 0.0;
          }
        }

        if ((D_nonneg_size[0] == 0) || (D_nonneg_size[1] == 0)) {
        } else {
          ixstart = D_nonneg_size[0] * (D_nonneg_size[1] - 1);
          for (cr = 0; cr <= ixstart; cr += D_nonneg_size[0]) {
            i4 = cr + D_nonneg_size[0];
            for (ic = cr; ic + 1 <= i4; ic++) {
              y->data[ic] = 0.0;
            }
          }

          n = 0;
          for (cr = 0; cr <= ixstart; cr += D_nonneg_size[0]) {
            ar = 0;
            i4 = n + D_nonneg_size[1];
            for (ib = n; ib + 1 <= i4; ib++) {
              if (D_nonneg_data[ib] != 0.0) {
                ia = ar;
                i5 = cr + D_nonneg_size[0];
                for (ic = cr; ic + 1 <= i5; ic++) {
                  ia++;
                  y->data[ic] += D_nonneg_data[ib] * D_nonneg_data[ia - 1];
                }
              }

              ar += D_nonneg_size[0];
            }

            n += D_nonneg_size[1];
          }
        }
      }

      b_diag(y, r1);
      i4 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = r1->size[0];
      b_y->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        b_y->data[i4] = r1->data[i4];
      }

      i4 = r4->size[0] * r4->size[1];
      r4->size[0] = b_y->size[0];
      r4->size[1] = b_y->size[1];
      emxEnsureCapacity((emxArray__common *)r4, i4, (int32_T)sizeof(real_T));
      loop_ub = b_y->size[0] * b_y->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        r4->data[i4] = 1.0 / b_y->data[i4];
      }

      b_diag(r4, r1);
      i4 = Binv->size[0] * Binv->size[1];
      Binv->size[0] = r1->size[0];
      Binv->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)Binv, i4, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        Binv->data[i4] = r1->data[i4];
      }

      a_size[0] = u_size[1];
      a_size[1] = 10;
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = u_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          b_a_data[i5 + a_size[0] * i4] = u_data[i4 + 10 * i5];
        }
      }

      outsz[0] = a_size[0];
      outsz[1] = u_size[1];
      D_nonneg_size[0] = outsz[0];
      D_nonneg_size[1] = outsz[1];
      loop_ub = D_nonneg_size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstart = D_nonneg_size[0] - 1;
        for (i5 = 0; i5 <= ixstart; i5++) {
          D_nonneg_data[i5 + D_nonneg_size[0] * i4] = 0.0;
        }
      }

      if ((a_size[0] == 0) || (u_size[1] == 0)) {
      } else {
        D_nonneg_data[0] = 0.0;
        ar = 0;
        for (ib = 0; ib + 1 < 11; ib++) {
          if (u_data[ib] != 0.0) {
            D_nonneg_data[0] += u_data[ib] * b_a_data[ar];
          }

          ar++;
        }
      }

      eye((real_T)D_nonneg_size[0], D_tree_data, D_tree_size);
      i4 = a->size[0] * a->size[1];
      a->size[0] = u_size[1];
      a->size[1] = 10;
      emxEnsureCapacity((emxArray__common *)a, i4, (int32_T)sizeof(real_T));
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = u_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          a->data[i5 + a->size[0] * i4] = u_data[i4 + 10 * i5];
        }
      }

      if (Binv->size[0] == 1) {
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = a->size[0];
        c_y->size[1] = Binv->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        loop_ub = a->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = Binv->size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            c_y->data[i4 + c_y->size[0] * i5] = 0.0;
            for (i6 = 0; i6 < 10; i6++) {
              c_y->data[i4 + c_y->size[0] * i5] += a->data[i4 + a->size[0] * i6]
                * Binv->data[i6 + Binv->size[0] * i5];
            }
          }
        }
      } else {
        outsz[0] = a->size[0];
        outsz[1] = Binv->size[1];
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = outsz[0];
        c_y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        m = a->size[0];
        i4 = c_y->size[0] * c_y->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        loop_ub = c_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = c_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            c_y->data[i5 + c_y->size[0] * i4] = 0.0;
          }
        }

        eml_xgemm(m, Binv->size[1], a, m, Binv, c_y, m);
      }

      if (c_y->size[1] == 1) {
        i4 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = c_y->size[0];
        d_y->size[1] = u_size[1];
        emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
        loop_ub = c_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = u_size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            d_y->data[i4 + d_y->size[0] * i5] = 0.0;
            n = c_y->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              d_y->data[i4 + d_y->size[0] * i5] += c_y->data[i4 + c_y->size[0] *
                i6] * u_data[i6 + 10 * i5];
            }
          }
        }
      } else {
        outsz[0] = c_y->size[0];
        outsz[1] = u_size[1];
        i4 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = outsz[0];
        d_y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
        m = c_y->size[0];
        i4 = d_y->size[0] * d_y->size[1];
        emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
        loop_ub = d_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = d_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            d_y->data[i5 + d_y->size[0] * i4] = 0.0;
          }
        }

        b_eml_xgemm(m, u_size[1], 0, c_y, m, u_data, u_size, 0, d_y, m);
      }

      n = D_tree_size[0];
      iy = D_tree_size[1];
      loop_ub = n * iy - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        D_tree_data[i4] += d_y->data[i4];
      }

      /* Break appart initial arrow% */
      /* The dimensions of the diagonal portion are euqal to the number of rows that create the trees in A */
      i4 = b_D_tree->size[0] * b_D_tree->size[1];
      b_D_tree->size[0] = (int32_T)row_v_data[0];
      b_D_tree->size[1] = (int32_T)row_v_data[0];
      emxEnsureCapacity((emxArray__common *)b_D_tree, i4, (int32_T)sizeof(real_T));
      loop_ub = (int32_T)row_v_data[0] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstart = (int32_T)row_v_data[0] - 1;
        for (i5 = 0; i5 <= ixstart; i5++) {
          b_D_tree->data[i5 + b_D_tree->size[0] * i4] = D_tree_data[i5 +
            D_tree_size[0] * i4];
        }
      }

      b_diag(b_D_tree, r1);
      i4 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = r1->size[0];
      b_y->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, i4, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        b_y->data[i4] = r1->data[i4];
      }

      i4 = r3->size[0] * r3->size[1];
      r3->size[0] = b_y->size[0];
      r3->size[1] = b_y->size[1];
      emxEnsureCapacity((emxArray__common *)r3, i4, (int32_T)sizeof(real_T));
      loop_ub = b_y->size[0] * b_y->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        r3->data[i4] = 1.0 / b_y->data[i4];
      }

      b_diag(r3, r1);
      i4 = A_11inv->size[0] * A_11inv->size[1];
      A_11inv->size[0] = r1->size[0];
      A_11inv->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)A_11inv, i4, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        A_11inv->data[i4] = r1->data[i4];
      }

      if ((int32_T)row_v_data[0] + 1 > D_tree_size[0]) {
        i4 = 0;
        i5 = 0;
      } else {
        i4 = (int32_T)row_v_data[0];
        i5 = D_tree_size[0];
      }

      A_12_size[0] = (int32_T)row_v_data[0];
      A_12_size[1] = i5 - i4;
      loop_ub = (i5 - i4) - 1;
      for (i5 = 0; i5 <= loop_ub; i5++) {
        ixstart = (int32_T)row_v_data[0] - 1;
        for (i6 = 0; i6 <= ixstart; i6++) {
          A_12_data[i6 + A_12_size[0] * i5] = D_tree_data[(i4 + i5) +
            D_tree_size[0] * i6];
        }
      }

      if ((int32_T)row_v_data[0] + 1 > D_tree_size[0]) {
        i4 = 0;
        i5 = 0;
      } else {
        i4 = (int32_T)row_v_data[0];
        i5 = D_tree_size[0];
      }

      if ((int32_T)row_v_data[0] + 1 > D_tree_size[1]) {
        i6 = 0;
        i7 = 0;
      } else {
        i6 = (int32_T)row_v_data[0];
        i7 = D_tree_size[1];
      }

      D_nonneg_size[0] = A_12_size[0];
      D_nonneg_size[1] = A_12_size[1];
      loop_ub = A_12_size[0] * A_12_size[1] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        D_nonneg_data[ixstop] = A_12_data[ixstop];
      }

      for (ixstop = 0; ixstop < 2; ixstop++) {
        outsz[ixstop] = A_12_size[ixstop];
      }

      outsz[0] = 1;
      n = 1;
      iy = -1;
      for (i = 1; i <= A_12_size[1]; i++) {
        ixstart = n;
        ixstop = (n + A_12_size[0]) - 1;
        alpha_2 = A_12_data[n - 1];
        if (A_12_size[0] > 1) {
          if (rtIsNaN(A_12_data[n - 1])) {
            cr = n;
            exitg7 = FALSE;
            while ((exitg7 == 0U) && (cr + 1 <= ixstop)) {
              ixstart = cr + 1;
              if (!rtIsNaN(A_12_data[cr])) {
                alpha_2 = A_12_data[cr];
                exitg7 = TRUE;
              } else {
                cr++;
              }
            }
          }

          if (ixstart < ixstop) {
            while (ixstart + 1 <= ixstop) {
              if (A_12_data[ixstart] > alpha_2) {
                alpha_2 = A_12_data[ixstart];
              }

              ixstart++;
            }
          }
        }

        iy++;
        D[iy] = alpha_2;
        n += A_12_size[0];
      }

      ixstart = 1;
      alpha_2 = D[0];
      if (outsz[1] > 1) {
        if (rtIsNaN(D[0])) {
          n = 2;
          exitg6 = FALSE;
          while ((exitg6 == 0U) && (n <= outsz[1])) {
            ixstart = n;
            if (!rtIsNaN(D[n - 1])) {
              alpha_2 = D[n - 1];
              exitg6 = TRUE;
            } else {
              n++;
            }
          }
        }

        if (ixstart < outsz[1]) {
          while (ixstart + 1 <= outsz[1]) {
            if (D[ixstart] > alpha_2) {
              alpha_2 = D[ixstart];
            }

            ixstart++;
          }
        }
      }

      for (ixstop = 0; ixstop < 2; ixstop++) {
        outsz[ixstop] = A_12_size[ixstop];
      }

      outsz[0] = 1;
      n = 0;
      iy = -1;
      for (i = 1; i <= A_12_size[1]; i++) {
        ixstart = n;
        ixstop = n + A_12_size[0];
        epsilon = A_12_data[n];
        if (A_12_size[0] > 1) {
          if (rtIsNaN(A_12_data[n])) {
            cr = n + 1;
            exitg5 = FALSE;
            while ((exitg5 == 0U) && (cr + 1 <= ixstop)) {
              ixstart = cr;
              if (!rtIsNaN(A_12_data[cr])) {
                epsilon = A_12_data[cr];
                exitg5 = TRUE;
              } else {
                cr++;
              }
            }
          }

          if (ixstart + 1 < ixstop) {
            for (cr = ixstart + 1; cr + 1 <= ixstop; cr++) {
              if (A_12_data[cr] < epsilon) {
                epsilon = A_12_data[cr];
              }
            }
          }
        }

        iy++;
        D[iy] = epsilon;
        n += A_12_size[0];
      }

      ixstart = 1;
      epsilon = D[0];
      if (outsz[1] > 1) {
        if (rtIsNaN(D[0])) {
          n = 2;
          exitg4 = FALSE;
          while ((exitg4 == 0U) && (n <= outsz[1])) {
            ixstart = n;
            if (!rtIsNaN(D[n - 1])) {
              epsilon = D[n - 1];
              exitg4 = TRUE;
            } else {
              n++;
            }
          }
        }

        if (ixstart < outsz[1]) {
          while (ixstart + 1 <= outsz[1]) {
            if (D[ixstart] < epsilon) {
              epsilon = D[ixstart];
            }

            ixstart++;
          }
        }
      }

      i_y = fabs(alpha_2 - epsilon) / 100.0;
      b_abs(A_12_data, A_12_size, b_tmp_data, b_tmp_size);
      A_12_Psi_to_remove_size[0] = b_tmp_size[0];
      A_12_Psi_to_remove_size[1] = b_tmp_size[1];
      loop_ub = b_tmp_size[0] * b_tmp_size[1] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        A_12_Psi_to_remove_data[ixstop] = (b_tmp_data[ixstop] < i_y);
      }

      c_eml_li_find(A_12_Psi_to_remove_data, A_12_Psi_to_remove_size, e_tmp_data,
                    tmp_size);
      loop_ub = tmp_size[0] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        f_tmp_data[ixstop] = e_tmp_data[ixstop];
      }

      loop_ub = tmp_size[0] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        D_nonneg_data[f_tmp_data[ixstop] - 1] = 0.0;
      }

      /* Psi is the difficult portion to inverse */
      ixstop = D_tree->size[0] * D_tree->size[1];
      D_tree->size[0] = (int32_T)row_v_data[0];
      D_tree->size[1] = (int32_T)row_v_data[0];
      emxEnsureCapacity((emxArray__common *)D_tree, ixstop, (int32_T)sizeof
                        (real_T));
      loop_ub = (int32_T)row_v_data[0] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        ixstart = (int32_T)row_v_data[0] - 1;
        for (i = 0; i <= ixstart; i++) {
          D_tree->data[i + D_tree->size[0] * ixstop] = D_tree_data[i +
            D_tree_size[0] * ixstop];
        }
      }

      b_diag(D_tree, r1);
      ixstop = b_y->size[0] * b_y->size[1];
      b_y->size[0] = r1->size[0];
      b_y->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, ixstop, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        b_y->data[ixstop] = r1->data[ixstop];
      }

      ixstop = d_y->size[0] * d_y->size[1];
      d_y->size[0] = D_nonneg_size[1];
      d_y->size[1] = D_nonneg_size[0];
      emxEnsureCapacity((emxArray__common *)d_y, ixstop, (int32_T)sizeof(real_T));
      loop_ub = D_nonneg_size[0] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        ixstart = D_nonneg_size[1] - 1;
        for (i = 0; i <= ixstart; i++) {
          d_y->data[i + d_y->size[0] * ixstop] = D_nonneg_data[ixstop +
            D_nonneg_size[0] * i];
        }
      }

      ixstop = r2->size[0] * r2->size[1];
      r2->size[0] = b_y->size[0];
      r2->size[1] = b_y->size[1];
      emxEnsureCapacity((emxArray__common *)r2, ixstop, (int32_T)sizeof(real_T));
      loop_ub = b_y->size[0] * b_y->size[1] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        r2->data[ixstop] = 1.0 / b_y->data[ixstop];
      }

      b_diag(r2, r1);
      ixstop = c_b->size[0] * c_b->size[1];
      c_b->size[0] = r1->size[0];
      c_b->size[1] = r1->size[1];
      emxEnsureCapacity((emxArray__common *)c_b, ixstop, (int32_T)sizeof(real_T));
      loop_ub = r1->size[0] * r1->size[1] - 1;
      for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
        c_b->data[ixstop] = r1->data[ixstop];
      }

      if ((d_y->size[1] == 1) || (c_b->size[0] == 1)) {
        ixstop = y->size[0] * y->size[1];
        y->size[0] = d_y->size[0];
        y->size[1] = c_b->size[1];
        emxEnsureCapacity((emxArray__common *)y, ixstop, (int32_T)sizeof(real_T));
        loop_ub = d_y->size[0] - 1;
        for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
          ixstart = c_b->size[1] - 1;
          for (i = 0; i <= ixstart; i++) {
            y->data[ixstop + y->size[0] * i] = 0.0;
            n = d_y->size[1] - 1;
            for (iy = 0; iy <= n; iy++) {
              y->data[ixstop + y->size[0] * i] += d_y->data[ixstop + d_y->size[0]
                * iy] * c_b->data[iy + c_b->size[0] * i];
            }
          }
        }
      } else {
        iy = d_y->size[1];
        outsz[0] = d_y->size[0];
        outsz[1] = c_b->size[1];
        ixstop = y->size[0] * y->size[1];
        y->size[0] = outsz[0];
        y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)y, ixstop, (int32_T)sizeof(real_T));
        m = d_y->size[0];
        n = c_b->size[1];
        ixstop = y->size[0] * y->size[1];
        emxEnsureCapacity((emxArray__common *)y, ixstop, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
          ixstart = y->size[0] - 1;
          for (i = 0; i <= ixstart; i++) {
            y->data[i + y->size[0] * ixstop] = 0.0;
          }
        }

        if ((m == 0) || (n == 0)) {
        } else {
          cr = 0;
          while (cr <= 0) {
            for (ic = 1; ic <= m; ic++) {
              y->data[ic - 1] = 0.0;
            }

            cr = m;
          }

          n = 0;
          cr = 0;
          while (cr <= 0) {
            ar = 0;
            ixstop = n + iy;
            for (ib = n; ib + 1 <= ixstop; ib++) {
              if (c_b->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= m; ic++) {
                  ia++;
                  y->data[ic] += c_b->data[ib] * d_y->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
            cr = m;
          }
        }
      }

      if ((y->size[1] == 1) || (D_nonneg_size[0] == 1)) {
        ixstop = e_y->size[0] * e_y->size[1];
        e_y->size[0] = y->size[0];
        e_y->size[1] = D_nonneg_size[1];
        emxEnsureCapacity((emxArray__common *)e_y, ixstop, (int32_T)sizeof
                          (real_T));
        loop_ub = y->size[0] - 1;
        for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
          ixstart = D_nonneg_size[1] - 1;
          for (i = 0; i <= ixstart; i++) {
            e_y->data[ixstop + e_y->size[0] * i] = 0.0;
            n = y->size[1] - 1;
            for (iy = 0; iy <= n; iy++) {
              e_y->data[ixstop + e_y->size[0] * i] += y->data[ixstop + y->size[0]
                * iy] * D_nonneg_data[iy + D_nonneg_size[0] * i];
            }
          }
        }
      } else {
        iy = y->size[1];
        outsz[0] = y->size[0];
        outsz[1] = D_nonneg_size[1];
        ixstop = e_y->size[0] * e_y->size[1];
        e_y->size[0] = outsz[0];
        e_y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)e_y, ixstop, (int32_T)sizeof
                          (real_T));
        m = y->size[0];
        ixstop = e_y->size[0] * e_y->size[1];
        emxEnsureCapacity((emxArray__common *)e_y, ixstop, (int32_T)sizeof
                          (real_T));
        loop_ub = e_y->size[1] - 1;
        for (ixstop = 0; ixstop <= loop_ub; ixstop++) {
          ixstart = e_y->size[0] - 1;
          for (i = 0; i <= ixstart; i++) {
            e_y->data[i + e_y->size[0] * ixstop] = 0.0;
          }
        }

        if ((m == 0) || (D_nonneg_size[1] == 0)) {
        } else {
          ixstart = m * (D_nonneg_size[1] - 1);
          for (cr = 0; cr <= ixstart; cr += m) {
            ixstop = cr + m;
            for (ic = cr; ic + 1 <= ixstop; ic++) {
              e_y->data[ic] = 0.0;
            }
          }

          n = 0;
          for (cr = 0; cr <= ixstart; cr += m) {
            ar = 0;
            ixstop = n + iy;
            for (ib = n; ib + 1 <= ixstop; ib++) {
              if (D_nonneg_data[ib] != 0.0) {
                ia = ar;
                i = cr + m;
                for (ic = cr; ic + 1 <= i; ic++) {
                  ia++;
                  e_y->data[ic] += D_nonneg_data[ib] * y->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
          }
        }
      }

      D_nonneg_size[0] = i5 - i4;
      D_nonneg_size[1] = i7 - i6;
      loop_ub = (i7 - i6) - 1;
      for (i7 = 0; i7 <= loop_ub; i7++) {
        ixstart = (i5 - i4) - 1;
        for (ixstop = 0; ixstop <= ixstart; ixstop++) {
          D_nonneg_data[ixstop + D_nonneg_size[0] * i7] = D_tree_data[(i4 +
            ixstop) + D_tree_size[0] * (i6 + i7)] - e_y->data[ixstop + e_y->
            size[0] * i7];
        }
      }

      eye((real_T)D_nonneg_size[1], b_tmp_data, b_tmp_size);
      b_D_nonneg_size[0] = D_nonneg_size[0];
      b_D_nonneg_size[1] = D_nonneg_size[1];
      loop_ub = D_nonneg_size[0] * D_nonneg_size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        b_D_nonneg_data[i4] = D_nonneg_data[i4];
      }

      mrdivide(b_tmp_data, b_tmp_size, b_D_nonneg_data, b_D_nonneg_size,
               D_nonneg_data, D_nonneg_size);

      /* Try to avoid the back slash! */
      /* Inversion of Arrow% */
      if ((A_11inv->size[1] == 1) || (A_12_size[0] == 1)) {
        i4 = y->size[0] * y->size[1];
        y->size[0] = A_11inv->size[0];
        y->size[1] = A_12_size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = A_11inv->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = A_12_size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i4 + y->size[0] * i5] = 0.0;
            n = A_11inv->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              y->data[i4 + y->size[0] * i5] += A_11inv->data[i4 + A_11inv->size
                [0] * i6] * A_12_data[i6 + A_12_size[0] * i5];
            }
          }
        }
      } else {
        iy = A_11inv->size[1];
        outsz[0] = A_11inv->size[0];
        outsz[1] = A_12_size[1];
        i4 = y->size[0] * y->size[1];
        y->size[0] = outsz[0];
        y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        m = A_11inv->size[0];
        i4 = y->size[0] * y->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i5 + y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (A_12_size[1] == 0)) {
        } else {
          ixstart = m * (A_12_size[1] - 1);
          for (cr = 0; cr <= ixstart; cr += m) {
            i4 = cr + m;
            for (ic = cr; ic + 1 <= i4; ic++) {
              y->data[ic] = 0.0;
            }
          }

          n = 0;
          for (cr = 0; cr <= ixstart; cr += m) {
            ar = 0;
            i4 = n + iy;
            for (ib = n; ib + 1 <= i4; ib++) {
              if (A_12_data[ib] != 0.0) {
                ia = ar;
                i5 = cr + m;
                for (ic = cr; ic + 1 <= i5; ic++) {
                  ia++;
                  y->data[ic] += A_12_data[ib] * A_11inv->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
          }
        }
      }

      if ((y->size[1] == 1) || (D_nonneg_size[0] == 1)) {
        i4 = e_y->size[0] * e_y->size[1];
        e_y->size[0] = y->size[0];
        e_y->size[1] = D_nonneg_size[1];
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = D_nonneg_size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            e_y->data[i4 + e_y->size[0] * i5] = 0.0;
            n = y->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              e_y->data[i4 + e_y->size[0] * i5] += y->data[i4 + y->size[0] * i6]
                * D_nonneg_data[i6 + D_nonneg_size[0] * i5];
            }
          }
        }
      } else {
        iy = y->size[1];
        outsz[0] = y->size[0];
        outsz[1] = D_nonneg_size[1];
        i4 = e_y->size[0] * e_y->size[1];
        e_y->size[0] = outsz[0];
        e_y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        m = y->size[0];
        i4 = e_y->size[0] * e_y->size[1];
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        loop_ub = e_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = e_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            e_y->data[i5 + e_y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (D_nonneg_size[1] == 0)) {
        } else {
          ixstart = m * (D_nonneg_size[1] - 1);
          for (cr = 0; cr <= ixstart; cr += m) {
            i4 = cr + m;
            for (ic = cr; ic + 1 <= i4; ic++) {
              e_y->data[ic] = 0.0;
            }
          }

          n = 0;
          for (cr = 0; cr <= ixstart; cr += m) {
            ar = 0;
            i4 = n + iy;
            for (ib = n; ib + 1 <= i4; ib++) {
              if (D_nonneg_data[ib] != 0.0) {
                ia = ar;
                i5 = cr + m;
                for (ic = cr; ic + 1 <= i5; ic++) {
                  ia++;
                  e_y->data[ic] += D_nonneg_data[ib] * y->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
          }
        }
      }

      i4 = c_b->size[0] * c_b->size[1];
      c_b->size[0] = A_12_size[1];
      c_b->size[1] = A_12_size[0];
      emxEnsureCapacity((emxArray__common *)c_b, i4, (int32_T)sizeof(real_T));
      loop_ub = A_12_size[0] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstart = A_12_size[1] - 1;
        for (i5 = 0; i5 <= ixstart; i5++) {
          c_b->data[i5 + c_b->size[0] * i4] = A_12_data[i4 + A_12_size[0] * i5];
        }
      }

      if ((e_y->size[1] == 1) || (c_b->size[0] == 1)) {
        i4 = y->size[0] * y->size[1];
        y->size[0] = e_y->size[0];
        y->size[1] = c_b->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = e_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = c_b->size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i4 + y->size[0] * i5] = 0.0;
            n = e_y->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              y->data[i4 + y->size[0] * i5] += e_y->data[i4 + e_y->size[0] * i6]
                * c_b->data[i6 + c_b->size[0] * i5];
            }
          }
        }
      } else {
        iy = e_y->size[1];
        outsz[0] = e_y->size[0];
        outsz[1] = c_b->size[1];
        i4 = y->size[0] * y->size[1];
        y->size[0] = outsz[0];
        y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        m = e_y->size[0];
        n = c_b->size[1];
        i4 = y->size[0] * y->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i5 + y->size[0] * i4] = 0.0;
          }
        }

        if (m == 0) {
        } else {
          ixstart = m * (n - 1);
          for (cr = 0; cr <= ixstart; cr += m) {
            i4 = cr + m;
            for (ic = cr; ic + 1 <= i4; ic++) {
              y->data[ic] = 0.0;
            }
          }

          n = 0;
          for (cr = 0; cr <= ixstart; cr += m) {
            ar = 0;
            i4 = n + iy;
            for (ib = n; ib + 1 <= i4; ib++) {
              if (c_b->data[ib] != 0.0) {
                ia = ar;
                i5 = cr + m;
                for (ic = cr; ic + 1 <= i5; ic++) {
                  ia++;
                  y->data[ic] += c_b->data[ib] * e_y->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
          }
        }
      }

      if ((y->size[1] == 1) || (A_11inv->size[0] == 1)) {
        i4 = e_y->size[0] * e_y->size[1];
        e_y->size[0] = y->size[0];
        e_y->size[1] = A_11inv->size[1];
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = A_11inv->size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            e_y->data[i4 + e_y->size[0] * i5] = 0.0;
            n = y->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              e_y->data[i4 + e_y->size[0] * i5] += y->data[i4 + y->size[0] * i6]
                * A_11inv->data[i6 + A_11inv->size[0] * i5];
            }
          }
        }
      } else {
        iy = y->size[1];
        outsz[0] = y->size[0];
        outsz[1] = A_11inv->size[1];
        i4 = e_y->size[0] * e_y->size[1];
        e_y->size[0] = outsz[0];
        e_y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        m = y->size[0];
        n = A_11inv->size[1];
        i4 = e_y->size[0] * e_y->size[1];
        emxEnsureCapacity((emxArray__common *)e_y, i4, (int32_T)sizeof(real_T));
        loop_ub = e_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = e_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            e_y->data[i5 + e_y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (n == 0)) {
        } else {
          cr = 0;
          while (cr <= 0) {
            for (ic = 1; ic <= m; ic++) {
              e_y->data[ic - 1] = 0.0;
            }

            cr = m;
          }

          n = 0;
          cr = 0;
          while (cr <= 0) {
            ar = 0;
            i4 = n + iy;
            for (ib = n; ib + 1 <= i4; ib++) {
              if (A_11inv->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= m; ic++) {
                  ia++;
                  e_y->data[ic] += A_11inv->data[ib] * y->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
            cr = m;
          }
        }
      }

      i4 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = A_11inv->size[0];
      d_y->size[1] = A_11inv->size[1];
      emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
      loop_ub = A_11inv->size[0] * A_11inv->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        d_y->data[i4] = -A_11inv->data[i4];
      }

      if ((d_y->size[1] == 1) || (A_12_size[0] == 1)) {
        i4 = y->size[0] * y->size[1];
        y->size[0] = d_y->size[0];
        y->size[1] = A_12_size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = d_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = A_12_size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i4 + y->size[0] * i5] = 0.0;
            n = d_y->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              y->data[i4 + y->size[0] * i5] += d_y->data[i4 + d_y->size[0] * i6]
                * A_12_data[i6 + A_12_size[0] * i5];
            }
          }
        }
      } else {
        outsz[0] = d_y->size[0];
        outsz[1] = A_12_size[1];
        i4 = y->size[0] * y->size[1];
        y->size[0] = outsz[0];
        y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        m = d_y->size[0];
        i4 = y->size[0] * y->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i5 + y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (A_12_size[1] == 0)) {
        } else {
          ixstart = m * (A_12_size[1] - 1);
          for (cr = 0; cr <= ixstart; cr += m) {
            i4 = cr + m;
            for (ic = cr; ic + 1 <= i4; ic++) {
              y->data[ic] = 0.0;
            }
          }
        }
      }

      if ((y->size[1] == 1) || (D_nonneg_size[0] == 1)) {
        i4 = c_b->size[0] * c_b->size[1];
        c_b->size[0] = y->size[0];
        c_b->size[1] = D_nonneg_size[1];
        emxEnsureCapacity((emxArray__common *)c_b, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = D_nonneg_size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            c_b->data[i4 + c_b->size[0] * i5] = 0.0;
            n = y->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              c_b->data[i4 + c_b->size[0] * i5] += y->data[i4 + y->size[0] * i6]
                * D_nonneg_data[i6 + D_nonneg_size[0] * i5];
            }
          }
        }
      } else {
        iy = y->size[1];
        outsz[0] = y->size[0];
        outsz[1] = D_nonneg_size[1];
        i4 = c_b->size[0] * c_b->size[1];
        c_b->size[0] = outsz[0];
        c_b->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)c_b, i4, (int32_T)sizeof(real_T));
        m = y->size[0];
        i4 = c_b->size[0] * c_b->size[1];
        emxEnsureCapacity((emxArray__common *)c_b, i4, (int32_T)sizeof(real_T));
        loop_ub = c_b->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = c_b->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            c_b->data[i5 + c_b->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (D_nonneg_size[1] == 0)) {
        } else {
          ixstart = m * (D_nonneg_size[1] - 1);
          for (cr = 0; cr <= ixstart; cr += m) {
            i4 = cr + m;
            for (ic = cr; ic + 1 <= i4; ic++) {
              c_b->data[ic] = 0.0;
            }
          }

          n = 0;
          for (cr = 0; cr <= ixstart; cr += m) {
            ar = 0;
            i4 = n + iy;
            for (ib = n; ib + 1 <= i4; ib++) {
              if (D_nonneg_data[ib] != 0.0) {
                ia = ar;
                i5 = cr + m;
                for (ic = cr; ic + 1 <= i5; ic++) {
                  ia++;
                  c_b->data[ic] += D_nonneg_data[ib] * y->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
          }
        }
      }

      i4 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = A_12_size[1];
      d_y->size[1] = A_12_size[0];
      emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
      loop_ub = A_12_size[0] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstart = A_12_size[1] - 1;
        for (i5 = 0; i5 <= ixstart; i5++) {
          d_y->data[i5 + d_y->size[0] * i4] = -A_12_data[i4 + A_12_size[0] * i5];
        }
      }

      if ((d_y->size[1] == 1) || (A_11inv->size[0] == 1)) {
        i4 = y->size[0] * y->size[1];
        y->size[0] = d_y->size[0];
        y->size[1] = A_11inv->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = d_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = A_11inv->size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i4 + y->size[0] * i5] = 0.0;
            n = d_y->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              y->data[i4 + y->size[0] * i5] += d_y->data[i4 + d_y->size[0] * i6]
                * A_11inv->data[i6 + A_11inv->size[0] * i5];
            }
          }
        }
      } else {
        iy = d_y->size[1];
        outsz[0] = d_y->size[0];
        outsz[1] = A_11inv->size[1];
        i4 = y->size[0] * y->size[1];
        y->size[0] = outsz[0];
        y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        m = d_y->size[0];
        n = A_11inv->size[1];
        i4 = y->size[0] * y->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i5 + y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (n == 0)) {
        } else {
          cr = 0;
          while (cr <= 0) {
            for (ic = 1; ic <= m; ic++) {
              y->data[ic - 1] = 0.0;
            }

            cr = m;
          }

          n = 0;
          cr = 0;
          while (cr <= 0) {
            ar = 0;
            i4 = n + iy;
            for (ib = n; ib + 1 <= i4; ib++) {
              if (A_11inv->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= m; ic++) {
                  ia++;
                  y->data[ic] += A_11inv->data[ib] * d_y->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
            cr = m;
          }
        }
      }

      if ((D_nonneg_size[1] == 1) || (y->size[0] == 1)) {
        i4 = f_y->size[0] * f_y->size[1];
        f_y->size[0] = D_nonneg_size[0];
        f_y->size[1] = y->size[1];
        emxEnsureCapacity((emxArray__common *)f_y, i4, (int32_T)sizeof(real_T));
        loop_ub = D_nonneg_size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = y->size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            f_y->data[i4 + f_y->size[0] * i5] = 0.0;
            n = D_nonneg_size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              f_y->data[i4 + f_y->size[0] * i5] += D_nonneg_data[i4 +
                D_nonneg_size[0] * i6] * y->data[i6 + y->size[0] * i5];
            }
          }
        }
      } else {
        outsz[0] = D_nonneg_size[0];
        outsz[1] = y->size[1];
        i4 = f_y->size[0] * f_y->size[1];
        f_y->size[0] = outsz[0];
        f_y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)f_y, i4, (int32_T)sizeof(real_T));
        n = y->size[1];
        i4 = f_y->size[0] * f_y->size[1];
        emxEnsureCapacity((emxArray__common *)f_y, i4, (int32_T)sizeof(real_T));
        loop_ub = f_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = f_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            f_y->data[i5 + f_y->size[0] * i4] = 0.0;
          }
        }

        if ((D_nonneg_size[0] == 0) || (n == 0)) {
        } else {
          cr = 0;
          while (cr <= 0) {
            for (ic = 1; ic <= D_nonneg_size[0]; ic++) {
              f_y->data[ic - 1] = 0.0;
            }

            cr = D_nonneg_size[0];
          }

          n = 0;
          cr = 0;
          while (cr <= 0) {
            ar = 0;
            i4 = n + D_nonneg_size[1];
            for (ib = n; ib + 1 <= i4; ib++) {
              if (y->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= D_nonneg_size[0]; ic++) {
                  ia++;
                  f_y->data[ic] += y->data[ib] * D_nonneg_data[ia - 1];
                }
              }

              ar += D_nonneg_size[0];
            }

            n += D_nonneg_size[1];
            cr = D_nonneg_size[0];
          }
        }
      }

      i4 = Arrow_inv->size[0] * Arrow_inv->size[1];
      Arrow_inv->size[0] = A_11inv->size[0] + f_y->size[0];
      Arrow_inv->size[1] = A_11inv->size[1] + c_b->size[1];
      emxEnsureCapacity((emxArray__common *)Arrow_inv, i4, (int32_T)sizeof
                        (real_T));
      loop_ub = A_11inv->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstart = A_11inv->size[0] - 1;
        for (i5 = 0; i5 <= ixstart; i5++) {
          Arrow_inv->data[i5 + Arrow_inv->size[0] * i4] = A_11inv->data[i5 +
            A_11inv->size[0] * i4] + e_y->data[i5 + e_y->size[0] * i4];
        }
      }

      loop_ub = c_b->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstart = c_b->size[0] - 1;
        for (i5 = 0; i5 <= ixstart; i5++) {
          Arrow_inv->data[i5 + Arrow_inv->size[0] * (i4 + A_11inv->size[1])] =
            c_b->data[i5 + c_b->size[0] * i4];
        }
      }

      loop_ub = f_y->size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstart = f_y->size[0] - 1;
        for (i5 = 0; i5 <= ixstart; i5++) {
          Arrow_inv->data[(i5 + A_11inv->size[0]) + Arrow_inv->size[0] * i4] =
            f_y->data[i5 + f_y->size[0] * i4];
        }
      }

      loop_ub = D_nonneg_size[1] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        ixstart = D_nonneg_size[0] - 1;
        for (i5 = 0; i5 <= ixstart; i5++) {
          Arrow_inv->data[(i5 + A_11inv->size[0]) + Arrow_inv->size[0] * (i4 +
            f_y->size[1])] = D_nonneg_data[i5 + D_nonneg_size[0] * i4];
        }
      }

      if (Binv->size[1] == 1) {
        i4 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = Binv->size[0];
        d_y->size[1] = u_size[1];
        emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
        loop_ub = Binv->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = u_size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            d_y->data[i4 + d_y->size[0] * i5] = 0.0;
            n = Binv->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              d_y->data[i4 + d_y->size[0] * i5] += Binv->data[i4 + Binv->size[0]
                * i6] * u_data[i6 + 10 * i5];
            }
          }
        }
      } else {
        outsz[0] = Binv->size[0];
        outsz[1] = u_size[1];
        i4 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = outsz[0];
        d_y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
        m = Binv->size[0];
        i4 = d_y->size[0] * d_y->size[1];
        emxEnsureCapacity((emxArray__common *)d_y, i4, (int32_T)sizeof(real_T));
        loop_ub = d_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = d_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            d_y->data[i5 + d_y->size[0] * i4] = 0.0;
          }
        }

        b_eml_xgemm(m, u_size[1], 0, Binv, m, u_data, u_size, 0, d_y, m);
      }

      if ((d_y->size[1] == 1) || (Arrow_inv->size[0] == 1)) {
        i4 = y->size[0] * y->size[1];
        y->size[0] = d_y->size[0];
        y->size[1] = Arrow_inv->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = d_y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = Arrow_inv->size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i4 + y->size[0] * i5] = 0.0;
            n = d_y->size[1] - 1;
            for (i6 = 0; i6 <= n; i6++) {
              y->data[i4 + y->size[0] * i5] += d_y->data[i4 + d_y->size[0] * i6]
                * Arrow_inv->data[i6 + Arrow_inv->size[0] * i5];
            }
          }
        }
      } else {
        iy = d_y->size[1];
        outsz[0] = d_y->size[0];
        outsz[1] = Arrow_inv->size[1];
        i4 = y->size[0] * y->size[1];
        y->size[0] = outsz[0];
        y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        m = d_y->size[0];
        n = Arrow_inv->size[1];
        i4 = y->size[0] * y->size[1];
        emxEnsureCapacity((emxArray__common *)y, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            y->data[i5 + y->size[0] * i4] = 0.0;
          }
        }

        if ((m == 0) || (n == 0)) {
        } else {
          ixstart = m * (n - 1);
          for (cr = 0; cr <= ixstart; cr += m) {
            i4 = cr + m;
            for (ic = cr; ic + 1 <= i4; ic++) {
              y->data[ic] = 0.0;
            }
          }

          n = 0;
          for (cr = 0; cr <= ixstart; cr += m) {
            ar = 0;
            i4 = n + iy;
            for (ib = n; ib + 1 <= i4; ib++) {
              if (Arrow_inv->data[ib] != 0.0) {
                ia = ar;
                i5 = cr + m;
                for (ic = cr; ic + 1 <= i5; ic++) {
                  ia++;
                  y->data[ic] += Arrow_inv->data[ib] * d_y->data[ia - 1];
                }
              }

              ar += m;
            }

            n += iy;
          }
        }
      }

      a_size[0] = u_size[1];
      a_size[1] = 10;
      for (i4 = 0; i4 < 10; i4++) {
        loop_ub = u_size[1] - 1;
        for (i5 = 0; i5 <= loop_ub; i5++) {
          b_a_data[i5 + a_size[0] * i4] = u_data[i4 + 10 * i5];
        }
      }

      if ((y->size[1] == 1) || (a_size[0] == 1)) {
        i4 = a->size[0] * a->size[1];
        a->size[0] = y->size[0];
        a->size[1] = 10;
        emxEnsureCapacity((emxArray__common *)a, i4, (int32_T)sizeof(real_T));
        loop_ub = y->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          for (i5 = 0; i5 < 10; i5++) {
            a->data[i4 + a->size[0] * i5] = 0.0;
            ixstart = y->size[1] - 1;
            for (i6 = 0; i6 <= ixstart; i6++) {
              a->data[i4 + a->size[0] * i5] += y->data[i4 + y->size[0] * i6] *
                b_a_data[i6 + a_size[0] * i5];
            }
          }
        }
      } else {
        iy = y->size[1];
        outsz[0] = y->size[0];
        outsz[1] = 10;
        i4 = a->size[0] * a->size[1];
        a->size[0] = outsz[0];
        a->size[1] = 10;
        emxEnsureCapacity((emxArray__common *)a, i4, (int32_T)sizeof(real_T));
        m = y->size[0];
        i4 = a->size[0] * a->size[1];

        a->size[1] = 10;
        emxEnsureCapacity((emxArray__common *)a, i4, (int32_T)sizeof(real_T));
        for (i4 = 0; i4 < 10; i4++) {
          loop_ub = a->size[0] - 1;
          for (i5 = 0; i5 <= loop_ub; i5++) {
            a->data[i5 + a->size[0] * i4] = 0.0;
          }
        }

        c_eml_xgemm(m, iy, y, m, b_a_data, a_size, iy, a, m);
      }

      if (Binv->size[0] == 1) {
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = a->size[0];
        c_y->size[1] = Binv->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        loop_ub = a->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = Binv->size[1] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            c_y->data[i4 + c_y->size[0] * i5] = 0.0;
            for (i6 = 0; i6 < 10; i6++) {
              c_y->data[i4 + c_y->size[0] * i5] += a->data[i4 + a->size[0] * i6]
                * Binv->data[i6 + Binv->size[0] * i5];
            }
          }
        }
      } else {
        outsz[0] = a->size[0];
        outsz[1] = Binv->size[1];
        i4 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = outsz[0];
        c_y->size[1] = outsz[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        m = a->size[0];
        i4 = c_y->size[0] * c_y->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i4, (int32_T)sizeof(real_T));
        loop_ub = c_y->size[1] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          ixstart = c_y->size[0] - 1;
          for (i5 = 0; i5 <= ixstart; i5++) {
            c_y->data[i5 + c_y->size[0] * i4] = 0.0;
          }
        }

        eml_xgemm(m, Binv->size[1], a, m, Binv, c_y, m);
      }

      i4 = Binv->size[0] * Binv->size[1];

      emxEnsureCapacity((emxArray__common *)Binv, i4, (int32_T)sizeof(real_T));
      n = Binv->size[0];
      iy = Binv->size[1];
      loop_ub = n * iy - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        Binv->data[i4] -= c_y->data[i4];
      }

      for (i4 = 0; i4 < 10; i4++) {
        epsilon = 0.0;
        for (i5 = 0; i5 < 10; i5++) {
          epsilon += ADDA[i4 + 10 * i5] * dx_1[i5];
        }

        r[i4] = c[i4] - epsilon;
      }

      if (Binv->size[1] == 1) {
        i4 = z->size[0];
        z->size[0] = Binv->size[0];
        emxEnsureCapacity((emxArray__common *)z, i4, (int32_T)sizeof(real_T));
        loop_ub = Binv->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          z->data[i4] = 0.0;
          for (i5 = 0; i5 < 10; i5++) {
            z->data[i4] += Binv->data[i4 + Binv->size[0] * i5] * r[i5];
          }
        }
      } else {
        outsz[0] = Binv->size[0];
        outsz[1] = 1;
        i4 = z->size[0];
        z->size[0] = outsz[0];
        emxEnsureCapacity((emxArray__common *)z, i4, (int32_T)sizeof(real_T));
        m = Binv->size[0];
        n = z->size[0];
        i4 = z->size[0];
        z->size[0] = n;
        emxEnsureCapacity((emxArray__common *)z, i4, (int32_T)sizeof(real_T));
        loop_ub = n - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          z->data[i4] = 0.0;
        }

        d_eml_xgemm(m, 0, Binv, m, r, 0, z, m);
      }

      i4 = p->size[0];
      p->size[0] = z->size[0];
      emxEnsureCapacity((emxArray__common *)p, i4, (int32_T)sizeof(real_T));
      loop_ub = z->size[0] - 1;
      for (i4 = 0; i4 <= loop_ub; i4++) {
        p->data[i4] = z->data[i4];
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
        if (z->size[0] == 1) {
          i_y = 0.0;
          for (i4 = 0; i4 < 10; i4++) {
            i_y += r[i4] * z->data[i4];
          }
        } else {
          i_y = 0.0;
          for (iy = 0; iy < 10; iy++) {
            i_y += r[iy] * z->data[iy];
          }
        }

        i4 = b_a->size[0] * b_a->size[1];
        b_a->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)b_a, i4, (int32_T)sizeof(real_T));
        n = p->size[0];
        i4 = b_a->size[0] * b_a->size[1];
        b_a->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_a, i4, (int32_T)sizeof(real_T));
        loop_ub = p->size[0] - 1;
        for (i4 = 0; i4 <= loop_ub; i4++) {
          b_a->data[i4] = p->data[i4];
        }

        if (b_a->size[1] == 1) {
          for (i4 = 0; i4 < 10; i4++) {
            j_y[i4] = 0.0;
            for (i5 = 0; i5 < 10; i5++) {
              epsilon = j_y[i4] + b_a->data[i5] * ADDA[i5 + 10 * i4];
              j_y[i4] = epsilon;
            }
          }
        } else {
          iy = b_a->size[1];
          memset(&j_y[0], 0, 10U * sizeof(real_T));
          for (cr = 0; cr < 10; cr++) {
            for (ic = cr; ic + 1 <= cr + 1; ic++) {
              j_y[ic] = 0.0;
            }
          }

          n = 0;
          for (cr = 0; cr < 10; cr++) {
            ar = 0;
            i4 = n + iy;
            for (ib = n; ib + 1 <= i4; ib++) {
              if (ADDA[ib] != 0.0) {
                ia = ar;
                for (ic = cr; ic + 1 <= cr + 1; ic++) {
                  ia++;
                  j_y[ic] += ADDA[ib] * b_a->data[ia - 1];
                }
              }

              ar++;
            }

            n += iy;
          }
        }

        if (p->size[0] == 1) {
          epsilon = 0.0;
          for (i4 = 0; i4 < 10; i4++) {
            epsilon += j_y[i4] * p->data[i4];
          }
        } else {
          epsilon = 0.0;
          for (iy = 0; iy < 10; iy++) {
            epsilon += j_y[iy] * p->data[iy];
          }
        }

        epsilon = i_y / epsilon;

        /* dx_1_new2=dx_1+(alpha_1*p); */
        for (n = 0; n < 10; n++) {
          dx_1_new[n] = dx_1[n] + epsilon * p->data[n];
        }

        /*                          dx_1_new2=dx_1+alpha_1*p; */
        /*                          test=all(dx_1_new2==dx_1_new) */
        for (i4 = 0; i4 < 100; i4++) {
          g_y[i4] = epsilon * ADDA[i4];
        }

        if (p->size[0] == 1) {
          for (i4 = 0; i4 < 10; i4++) {
            r_new[i4] = 0.0;
            for (i5 = 0; i5 < 10; i5++) {
              epsilon = r_new[i4] + g_y[i4 + 10 * i5] * p->data[i5];
              r_new[i4] = epsilon;
            }
          }
        } else {
          memset(&r_new[0], 0, 10U * sizeof(real_T));
          if (10 == p->size[0]) {
            for (i4 = 0; i4 < 10; i4++) {
              r_new[i4] = 0.0;
              for (i5 = 0; i5 < 10; i5++) {
                epsilon = r_new[i4] + g_y[i4 + 10 * i5] * p->data[i5];
                r_new[i4] = epsilon;
              }
            }
          } else {
            memset(&r_new[0], 0, 10U * sizeof(real_T));
            ar = 0;
            for (ib = 0; ib < 10; ib++) {
              if (p->data[ib] != 0.0) {
                ia = ar;
                for (ic = 0; ic < 10; ic++) {
                  ia++;
                  epsilon = r_new[ic] + p->data[ib] * g_y[ia - 1];
                  r_new[ic] = epsilon;
                }
              }

              ar += 10;
            }
          }
        }

        i_y = 0.0;
        for (i4 = 0; i4 < 10; i4++) {
          epsilon = r[i4] - r_new[i4];
          i_y += epsilon * epsilon;
          r_new[i4] = epsilon;
        }

        if (sqrt(i_y) < 1.0E-5) {
          /* norm of rnew */
          /* sprintf('Iteration count for PCG: %d, with value', i, sqrt(r_new'*r_new)) */
          exitg3 = TRUE;
        } else {
          if (Binv->size[1] == 1) {
            i4 = z_new->size[0];
            z_new->size[0] = Binv->size[0];
            emxEnsureCapacity((emxArray__common *)z_new, i4, (int32_T)sizeof
                              (real_T));
            loop_ub = Binv->size[0] - 1;
            for (i4 = 0; i4 <= loop_ub; i4++) {
              z_new->data[i4] = 0.0;
              for (i5 = 0; i5 < 10; i5++) {
                z_new->data[i4] += Binv->data[i4 + Binv->size[0] * i5] *
                  r_new[i5];
              }
            }
          } else {
            outsz[0] = Binv->size[0];
            outsz[1] = 1;
            i4 = z_new->size[0];
            z_new->size[0] = outsz[0];
            emxEnsureCapacity((emxArray__common *)z_new, i4, (int32_T)sizeof
                              (real_T));
            m = Binv->size[0];
            n = z_new->size[0];
            i4 = z_new->size[0];
            z_new->size[0] = n;
            emxEnsureCapacity((emxArray__common *)z_new, i4, (int32_T)sizeof
                              (real_T));
            loop_ub = n - 1;
            for (i4 = 0; i4 <= loop_ub; i4++) {
              z_new->data[i4] = 0.0;
            }

            d_eml_xgemm(m, 0, Binv, m, r_new, 0, z_new, m);
          }

          i4 = b_a->size[0] * b_a->size[1];
          b_a->size[0] = 1;
          emxEnsureCapacity((emxArray__common *)b_a, i4, (int32_T)sizeof(real_T));
          n = z_new->size[0];
          i4 = b_a->size[0] * b_a->size[1];
          b_a->size[1] = n;
          emxEnsureCapacity((emxArray__common *)b_a, i4, (int32_T)sizeof(real_T));
          loop_ub = z_new->size[0] - 1;
          for (i4 = 0; i4 <= loop_ub; i4++) {
            b_a->data[i4] = z_new->data[i4];
          }

          if (b_a->size[1] == 1) {
            i_y = 0.0;
            for (i4 = 0; i4 < 10; i4++) {
              i_y += b_a->data[i4] * r_new[i4];
            }
          } else {
            iy = b_a->size[1];
            i_y = 0.0;
            if (iy < 1) {
            } else {
              for (n = 0; n + 1 <= iy; n++) {
                i_y += b_a->data[n] * r_new[n];
              }
            }
          }

          i4 = b_a->size[0] * b_a->size[1];
          b_a->size[0] = 1;
          emxEnsureCapacity((emxArray__common *)b_a, i4, (int32_T)sizeof(real_T));
          n = z->size[0];
          i4 = b_a->size[0] * b_a->size[1];
          b_a->size[1] = n;
          emxEnsureCapacity((emxArray__common *)b_a, i4, (int32_T)sizeof(real_T));
          loop_ub = z->size[0] - 1;
          for (i4 = 0; i4 <= loop_ub; i4++) {
            b_a->data[i4] = z->data[i4];
          }

          if (b_a->size[1] == 1) {
            epsilon = 0.0;
            for (i4 = 0; i4 < 10; i4++) {
              epsilon += b_a->data[i4] * r[i4];
            }
          } else {
            iy = b_a->size[1];
            epsilon = 0.0;
            if (iy < 1) {
            } else {
              for (n = 0; n + 1 <= iy; n++) {
                epsilon += b_a->data[n] * r[n];
              }
            }
          }

          epsilon = i_y / epsilon;
          i4 = p->size[0];
          p->size[0] = z_new->size[0];
          emxEnsureCapacity((emxArray__common *)p, i4, (int32_T)sizeof(real_T));
          loop_ub = z_new->size[0] - 1;
          for (i4 = 0; i4 <= loop_ub; i4++) {
            p->data[i4] = z_new->data[i4] + epsilon * p->data[i4];
          }

          for (n = 0; n < 10; n++) {
            r[n] = r_new[n];
            dx_1[n] = dx_1_new[n];
          }

          i4 = z->size[0];
          z->size[0] = z_new->size[0];
          emxEnsureCapacity((emxArray__common *)z, i4, (int32_T)sizeof(real_T));
          loop_ub = z_new->size[0] - 1;
          for (i4 = 0; i4 <= loop_ub; i4++) {
            z->data[i4] = z_new->data[i4];
          }

          i++;
        }
      }
    }

    /*         %% Step 5 Scaling */
    /*         %% Step 6 Calculating the interior point from the projection */
    for (i = 0; i < 17; i++) {
      D[i] = -w[i];
    }

    for (i4 = 0; i4 < 10; i4++) {
      for (i5 = 0; i5 < 17; i5++) {
        A_data[i5 + 17 * i4] = -A[i5 + 17 * i4];
      }
    }

    for (i4 = 0; i4 < 17; i4++) {
      D_data[i4] = 0.0;
      for (i5 = 0; i5 < 10; i5++) {
        D_data[i4] += A_data[i4 + 17 * i5] * dx_1[i5];
      }
    }

    rdivide(D, D_data, w);
    for (i = 0; i < 17; i++) {
      b_A[i] = (w[i] < 0.0);
    }

    eml_li_find(b_A, tmp_data, tmp_size);
    loop_ub = tmp_size[0] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
      w[tmp_data[i4] - 1] = 1.0E+6;
    }

    ixstart = 1;
    alpha_2 = w[0];
    if (rtIsNaN(w[0])) {
      n = 2;
      exitg2 = FALSE;
      while ((exitg2 == 0U) && (n < 18)) {
        ixstart = n;
        if (!rtIsNaN(w[n - 1])) {
          alpha_2 = w[n - 1];
          exitg2 = TRUE;
        } else {
          n++;
        }
      }
    }

    if (ixstart < 17) {
      while (ixstart + 1 < 18) {
        if (w[ixstart] < alpha_2) {
          alpha_2 = w[ixstart];
        }

        ixstart++;
      }
    }

    alpha_2 *= b_gamma;

    /*         %% Step 1 */
    i_y = 0.0;
    epsilon = 0.0;
    k_y = 0.0;
    for (i = 0; i < 10; i++) {
      b_r = x[i] + alpha_2 * dx_1[i];
      i_y += c[i] * b_r;
      epsilon += c[i] * x[i];
      k_y += c[i] * x[i];
      r[i] = b_r;
    }

    guard1 = FALSE;
    if ((i_y - epsilon) / k_y > 0.01) {
      memcpy(&x[0], &r[0], 10U * sizeof(real_T));
      guard1 = TRUE;
    } else {
      i_y = 0.0;
      epsilon = 0.0;
      k_y = 0.0;
      for (iy = 0; iy < 10; iy++) {
        i_y += c[iy] * r[iy];
        epsilon += c[iy] * x[iy];
        k_y += c[iy] * x[iy];
      }

      if ((i_y - epsilon) / k_y <= 0.01) {
        exitg1 = TRUE;
      } else {
        guard1 = TRUE;
      }
    }

    if (guard1 == TRUE) {
      b_check++;
    }
  }

  emxFree_real_T(&D_track);
  emxFree_real_T(&r5);
  emxFree_real_T(&r4);
  emxFree_real_T(&b_D_tree);
  emxFree_real_T(&r3);
  emxFree_real_T(&D_tree);
  emxFree_real_T(&r2);
  emxFree_real_T(&r1);
  emxFree_real_T(&b_a);
  emxFree_real_T(&f_y);
  emxFree_real_T(&e_y);
  emxFree_real_T(&c_b);
  emxFree_real_T(&d_y);
  emxFree_real_T(&c_y);
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
  emxFree_real_T(&b_b);
  emxFree_real_T(&r0);
  emxFree_real_T(&z_new);
  emxFree_real_T(&p);
  emxFree_real_T(&z);
  emxFree_real_T(&Arrow_inv);
  emxFree_real_T(&A_11inv);
  *IPMit = (real_T)check;
  if (Binv->size[1] == 1) {
    i0 = a->size[0] * a->size[1];
    a->size[0] = Binv->size[0];
    a->size[1] = 10;
    emxEnsureCapacity((emxArray__common *)a, i0, (int32_T)sizeof(real_T));
    loop_ub = Binv->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      for (i1 = 0; i1 < 10; i1++) {
        a->data[i0 + a->size[0] * i1] = 0.0;
        ixstart = Binv->size[1] - 1;
        for (i2 = 0; i2 <= ixstart; i2++) {
          a->data[i0 + a->size[0] * i1] += Binv->data[i0 + Binv->size[0] * i2] *
            ADDA[i2 + 10 * i1];
        }
      }
    }
  } else {
    outsz[0] = Binv->size[0];
    outsz[1] = 10;
    i0 = a->size[0] * a->size[1];
    a->size[0] = outsz[0];
    a->size[1] = 10;
    emxEnsureCapacity((emxArray__common *)a, i0, (int32_T)sizeof(real_T));
    m = Binv->size[0];
    i0 = a->size[0] * a->size[1];

    a->size[1] = 10;
    emxEnsureCapacity((emxArray__common *)a, i0, (int32_T)sizeof(real_T));
    for (i0 = 0; i0 < 10; i0++) {
      loop_ub = a->size[0] - 1;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        a->data[i1 + a->size[0] * i0] = 0.0;
      }
    }

    e_eml_xgemm(m, 0, Binv, m, ADDA, 0, a, m);
  }

  emxFree_real_T(&Binv);
  *CondNum = cond(a);

  /* CondNum=condition_num_MinvADDA; */
  emxFree_real_T(&a);
}

/* End of code generation (Function_PCG_Wood_clean_codegen.c) */
