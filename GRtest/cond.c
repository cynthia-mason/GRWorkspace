/*
 * cond.c
 *
 * Code generation for function 'cond'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen.h"
#include "cond.h"
#include "Function_PCG_Wood_clean_codegen_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_eml_xaxpy(int32_T n, real_T a, const emxArray_real_T *x, int32_T
  ix0, emxArray_real_T *y, int32_T iy0);
static void b_eml_xscal(int32_T n, real_T a, real_T x[10], int32_T ix0);
static void c_eml_xaxpy(int32_T n, real_T a, const emxArray_real_T *x, int32_T
  ix0, emxArray_real_T *y, int32_T iy0);
static real_T c_eml_xnrm2(int32_T n, const emxArray_real_T *x, int32_T ix0);
static real_T d_eml_xnrm2(int32_T n, const real_T x[10], int32_T ix0);
static real_T eml_div(real_T x, real_T y);
static void eml_xaxpy(int32_T n, real_T a, int32_T ix0, emxArray_real_T *y,
                      int32_T iy0);
static real_T eml_xdotc(int32_T n, const emxArray_real_T *x, int32_T ix0, const
  emxArray_real_T *y, int32_T iy0);
static void eml_xgesvd(const emxArray_real_T *A, real_T S_data[10], int32_T
  S_size[1]);
static void eml_xrotg(real_T *a, real_T *b, real_T *c, real_T *s);
static void eml_xscal(int32_T n, real_T a, emxArray_real_T *x, int32_T ix0);

/* Function Definitions */
static void b_eml_xaxpy(int32_T n, real_T a, const emxArray_real_T *x, int32_T
  ix0, emxArray_real_T *y, int32_T iy0)
{
  emxArray_real_T *b_y;
  int32_T ix;
  int32_T iy;
  int32_T k;
  b_emxInit_real_T(&b_y, 1);
  if ((n < 1) || (a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k <= n - 1; k++) {
      y->data[iy] += a * x->data[ix];
      ix++;
      iy++;
    }
  }

  emxFree_real_T(&b_y);
}

static void b_eml_xscal(int32_T n, real_T a, real_T x[10], int32_T ix0)
{
  int32_T i20;
  int32_T k;
  i20 = (ix0 + n) - 1;
  for (k = ix0; k <= i20; k++) {
    x[k - 1] *= a;
  }
}

static void c_eml_xaxpy(int32_T n, real_T a, const emxArray_real_T *x, int32_T
  ix0, emxArray_real_T *y, int32_T iy0)
{
  emxArray_real_T *b_y;
  int32_T ix;
  int32_T iy;
  int32_T k;
  emxInit_real_T(&b_y, 2);
  if ((n < 1) || (a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k <= n - 1; k++) {
      y->data[iy] += a * x->data[ix];
      ix++;
      iy++;
    }
  }

  emxFree_real_T(&b_y);
}

static real_T c_eml_xnrm2(int32_T n, const emxArray_real_T *x, int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  int32_T k;
  real_T absxk;
  real_T t;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = fabs(x->data[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x->data[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

static real_T d_eml_xnrm2(int32_T n, const real_T x[10], int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  int32_T k;
  real_T absxk;
  real_T t;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = fabs(x[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

static real_T eml_div(real_T x, real_T y)
{
  return x / y;
}

static void eml_xaxpy(int32_T n, real_T a, int32_T ix0, emxArray_real_T *y,
                      int32_T iy0)
{
  emxArray_real_T *b_y;
  int32_T ix;
  int32_T iy;
  int32_T k;
  emxInit_real_T(&b_y, 2);
  if ((n < 1) || (a == 0.0)) {
  } else {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k <= n - 1; k++) {
      y->data[iy] += a * y->data[ix];
      ix++;
      iy++;
    }
  }

  emxFree_real_T(&b_y);
}

static real_T eml_xdotc(int32_T n, const emxArray_real_T *x, int32_T ix0, const
  emxArray_real_T *y, int32_T iy0)
{
  real_T d;
  int32_T ix;
  int32_T iy;
  int32_T k;
  d = 0.0;
  if (n < 1) {
  } else {
    ix = ix0;
    iy = iy0;
    for (k = 1; k <= n; k++) {
      d += x->data[ix - 1] * y->data[iy - 1];
      ix++;
      iy++;
    }
  }

  return d;
}

static void eml_xgesvd(const emxArray_real_T *A, real_T S_data[10], int32_T
  S_size[1])
{
  emxArray_real_T *b_A;
  int32_T mm;
  int32_T kase;
  int32_T n;
  int32_T minnp;
  real_T s_data[10];
  real_T e[10];
  emxArray_real_T *work;
  int32_T nrt;
  int32_T nct;
  int32_T q;
  int32_T qp1;
  int32_T iter;
  real_T ztest0;
  int32_T qs;
  int32_T mm1;
  int32_T m;
  real_T ztest;
  real_T tiny;
  real_T snorm;
  boolean_T exitg3;
  boolean_T exitg2;
  real_T sn;
  real_T sm;
  real_T varargin_1[5];
  boolean_T exitg1;
  real_T sqds;
  real_T b;
  emxInit_real_T(&b_A, 2);
  mm = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = 10;
  emxEnsureCapacity((emxArray__common *)b_A, mm, (int32_T)sizeof(real_T));
  kase = A->size[0] * A->size[1] - 1;
  for (mm = 0; mm <= kase; mm++) {
    b_A->data[mm] = A->data[mm];
  }

  n = A->size[0];
  if (n <= 10) {
    minnp = n;
  } else {
    minnp = 10;
  }

  kase = n + 1;
  if (kase <= 10) {
  } else {
    kase = 10;
  }

  kase--;
  for (mm = 0; mm <= kase; mm++) {
    s_data[mm] = 0.0;
  }

  memset(&e[0], 0, 10U * sizeof(real_T));
  b_emxInit_real_T(&work, 1);
  mm = work->size[0];
  work->size[0] = n;
  emxEnsureCapacity((emxArray__common *)work, mm, (int32_T)sizeof(real_T));
  kase = n - 1;
  for (mm = 0; mm <= kase; mm++) {
    work->data[mm] = 0.0;
  }

  if (A->size[0] == 0) {
  } else {
    if (8 <= n) {
      nrt = 8;
    } else {
      nrt = n;
    }

    if (n < 1) {
      kase = 1;
    } else {
      kase = n;
    }

    nct = kase - 1;
    if (nct <= 10) {
    } else {
      nct = 10;
    }

    if (nct >= nrt) {
      mm = nct;
    } else {
      mm = nrt;
    }

    for (q = 0; q + 1 <= mm; q++) {
      qp1 = q + 2;
      kase = q + n * q;
      iter = n - q;
      if (q + 1 <= nct) {
        ztest0 = c_eml_xnrm2(iter, b_A, kase + 1);
        if (ztest0 > 0.0) {
          if (b_A->data[kase] < 0.0) {
            ztest0 = -ztest0;
          }

          s_data[q] = ztest0;
          eml_xscal(iter, eml_div(1.0, s_data[q]), b_A, kase + 1);
          b_A->data[kase]++;
          s_data[q] = -s_data[q];
        } else {
          s_data[q] = 0.0;
        }
      }

      for (qs = qp1; qs < 11; qs++) {
        mm1 = q + n * (qs - 1);
        if ((q + 1 <= nct) && (s_data[q] != 0.0)) {
          ztest0 = eml_xdotc(iter, b_A, kase + 1, b_A, mm1 + 1);
          ztest0 = -eml_div(ztest0, b_A->data[q + b_A->size[0] * q]);
          eml_xaxpy(iter, ztest0, kase + 1, b_A, mm1 + 1);
        }

        e[qs - 1] = b_A->data[mm1];
      }

      if (q + 1 <= nrt) {
        ztest0 = d_eml_xnrm2(9 - q, e, qp1);
        if (ztest0 == 0.0) {
          e[q] = 0.0;
        } else {
          if (e[qp1 - 1] < 0.0) {
            ztest0 = -ztest0;
          }

          e[q] = ztest0;
          b_eml_xscal(9 - q, eml_div(1.0, e[q]), e, qp1);
          e[qp1 - 1]++;
        }

        e[q] = -e[q];
        if ((qp1 <= n) && (e[q] != 0.0)) {
          for (kase = qp1; kase <= n; kase++) {
            work->data[kase - 1] = 0.0;
          }

          for (qs = qp1; qs < 11; qs++) {
            b_eml_xaxpy(iter - 1, e[qs - 1], b_A, qp1 + n * (qs - 1), work, qp1);
          }

          for (qs = qp1; qs < 11; qs++) {
            c_eml_xaxpy(iter - 1, eml_div(-e[qs - 1], e[qp1 - 1]), work, qp1,
                        b_A, qp1 + n * (qs - 1));
          }
        }
      }
    }

    m = n + 1;
    if (10 <= m) {
      m = 10;
    }

    if (nct < 10) {
      s_data[nct] = b_A->data[nct + b_A->size[0] * nct];
    }

    if (n < m) {
      s_data[m - 1] = 0.0;
    }

    if (nrt + 1 < m) {
      e[nrt] = b_A->data[nrt + b_A->size[0] * (m - 1)];
    }

    e[m - 1] = 0.0;
    for (q = 0; q + 1 <= m; q++) {
      if (s_data[q] != 0.0) {
        ztest = fabs(s_data[q]);
        ztest0 = eml_div(s_data[q], ztest);
        s_data[q] = ztest;
        if (q + 1 < m) {
          e[q] = eml_div(e[q], ztest0);
        }
      }

      if ((q + 1 < m) && (e[q] != 0.0)) {
        ztest = fabs(e[q]);
        ztest0 = eml_div(ztest, e[q]);
        e[q] = ztest;
        s_data[q + 1] *= ztest0;
      }
    }

    mm = m;
    iter = 0;
    tiny = eml_div(2.2250738585072014E-308, 2.2204460492503131E-16);
    snorm = 0.0;
    for (kase = 0; kase + 1 <= m; kase++) {
      ztest0 = fabs(s_data[kase]);
      ztest = fabs(e[kase]);
      if ((ztest0 >= ztest) || rtIsNaN(ztest)) {
        ztest = ztest0;
      }

      if ((snorm >= ztest) || rtIsNaN(ztest)) {
      } else {
        snorm = ztest;
      }
    }

    while ((m > 0) && (!(iter >= 75))) {
      q = m - 1;
      exitg3 = FALSE;
      while (!((exitg3 == 1U) || (q == 0))) {
        ztest0 = fabs(e[q - 1]);
        if ((ztest0 <= 2.2204460492503131E-16 * (fabs(s_data[q - 1]) + fabs
              (s_data[q]))) || (ztest0 <= tiny) || ((iter > 20) && (ztest0 <=
              2.2204460492503131E-16 * snorm))) {
          e[q - 1] = 0.0;
          exitg3 = TRUE;
        } else {
          q--;
        }
      }

      if (q == m - 1) {
        kase = 4;
      } else {
        qs = m;
        kase = m;
        exitg2 = FALSE;
        while ((exitg2 == 0U) && (kase >= q)) {
          qs = kase;
          if (kase == q) {
            exitg2 = TRUE;
          } else {
            ztest0 = 0.0;
            if (kase < m) {
              ztest0 = fabs(e[kase - 1]);
            }

            if (kase > q + 1) {
              ztest0 += fabs(e[kase - 2]);
            }

            ztest = fabs(s_data[kase - 1]);
            if ((ztest <= 2.2204460492503131E-16 * ztest0) || (ztest <= tiny)) {
              s_data[kase - 1] = 0.0;
              exitg2 = TRUE;
            } else {
              kase--;
            }
          }
        }

        if (qs == q) {
          kase = 3;
        } else if (qs == m) {
          kase = 1;
        } else {
          kase = 2;
          q = qs;
        }
      }

      switch (kase) {
       case 1:
        ztest = e[m - 2];
        e[m - 2] = 0.0;
        for (qs = m - 1; qs >= q + 1; qs--) {
          ztest0 = s_data[qs - 1];
          eml_xrotg(&ztest0, &ztest, &sm, &sn);
          s_data[qs - 1] = ztest0;
          if (qs > q + 1) {
            kase = qs - 2;
            ztest = -sn * e[kase];
            e[kase] *= sm;
          }
        }
        break;

       case 2:
        kase = q - 1;
        ztest = e[kase];
        e[kase] = 0.0;
        while (q + 1 <= m) {
          eml_xrotg(&s_data[q], &ztest, &sm, &sn);
          ztest = -sn * e[q];
          e[q] *= sm;
          q++;
        }
        break;

       case 3:
        mm1 = m - 2;
        varargin_1[0] = fabs(s_data[m - 1]);
        varargin_1[1] = fabs(s_data[mm1]);
        varargin_1[2] = fabs(e[mm1]);
        varargin_1[3] = fabs(s_data[q]);
        varargin_1[4] = fabs(e[q]);
        kase = 1;
        sn = varargin_1[0];
        if (rtIsNaN(varargin_1[0])) {
          qs = 2;
          exitg1 = FALSE;
          while ((exitg1 == 0U) && (qs < 6)) {
            kase = qs;
            if (!rtIsNaN(varargin_1[qs - 1])) {
              sn = varargin_1[qs - 1];
              exitg1 = TRUE;
            } else {
              qs++;
            }
          }
        }

        if (kase < 5) {
          while (kase + 1 < 6) {
            if (varargin_1[kase] > sn) {
              sn = varargin_1[kase];
            }

            kase++;
          }
        }

        sm = eml_div(s_data[m - 1], sn);
        ztest0 = eml_div(s_data[mm1], sn);
        ztest = eml_div(e[mm1], sn);
        sqds = eml_div(s_data[q], sn);
        b = eml_div((ztest0 + sm) * (ztest0 - sm) + ztest * ztest, 2.0);
        ztest0 = sm * ztest;
        ztest0 *= ztest0;
        ztest = 0.0;
        if ((b != 0.0) || (ztest0 != 0.0)) {
          ztest = sqrt(b * b + ztest0);
          if (b < 0.0) {
            ztest = -ztest;
          }

          ztest = eml_div(ztest0, b + ztest);
        }

        ztest += (sqds + sm) * (sqds - sm);
        ztest0 = sqds * eml_div(e[q], sn);
        for (qs = q; qs + 1 <= mm1 + 1; qs++) {
          kase = qs + 1;
          eml_xrotg(&ztest, &ztest0, &sm, &sn);
          if (qs + 1 > q + 1) {
            e[qs - 1] = ztest;
          }

          ztest0 = sm * s_data[qs];
          ztest = sn * e[qs];
          e[qs] = sm * e[qs] - sn * s_data[qs];
          b = s_data[kase];
          s_data[kase] *= sm;
          s_data[qs] = ztest0 + ztest;
          ztest0 = sn * b;
          eml_xrotg(&s_data[qs], &ztest0, &sm, &sn);
          ztest = sm * e[qs] + sn * s_data[kase];
          s_data[kase] = -sn * e[qs] + sm * s_data[kase];
          ztest0 = sn * e[kase];
          e[kase] *= sm;
        }

        e[m - 2] = ztest;
        iter++;
        break;

       default:
        if (s_data[q] < 0.0) {
          s_data[q] = -s_data[q];
        }

        qp1 = q + 1;
        while ((q + 1 < mm) && (s_data[q] < s_data[qp1])) {
          ztest = s_data[q];
          s_data[q] = s_data[qp1];
          s_data[qp1] = ztest;
          q = qp1;
          qp1++;
        }

        iter = 0;
        m--;
        break;
      }
    }
  }

  emxFree_real_T(&work);
  emxFree_real_T(&b_A);
  S_size[0] = minnp;
  for (qs = 0; qs + 1 <= minnp; qs++) {
    S_data[qs] = s_data[qs];
  }
}

static void eml_xrotg(real_T *a, real_T *b, real_T *c, real_T *s)
{
  real_T roe;
  real_T absa;
  real_T absb;
  real_T scale;
  real_T ads;
  real_T bds;
  roe = *b;
  absa = fabs(*a);
  absb = fabs(*b);
  if (absa > absb) {
    roe = *a;
  }

  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    ads = 0.0;
    scale = 0.0;
  } else {
    ads = absa / scale;
    bds = absb / scale;
    ads = scale * sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      ads = -ads;
    }

    *c = *a / ads;
    *s = *b / ads;
    if (absa > absb) {
      scale = *s;
    } else if (*c != 0.0) {
      scale = 1.0 / *c;
    } else {
      scale = 1.0;
    }
  }

  *a = ads;
  *b = scale;
}

static void eml_xscal(int32_T n, real_T a, emxArray_real_T *x, int32_T ix0)
{
  emxArray_real_T *b_x;
  int32_T i19;
  int32_T k;
  emxInit_real_T(&b_x, 2);
  i19 = (ix0 + n) - 1;
  for (k = ix0; k <= i19; k++) {
    x->data[k - 1] *= a;
  }

  emxFree_real_T(&b_x);
}

real_T cond(const emxArray_real_T *A)
{
  real_T c;
  int32_T U_size[1];
  real_T U_data[10];
  if (A->size[0] == 0) {
    c = 0.0;
  } else {
    eml_xgesvd(A, U_data, U_size);
    if (U_data[U_size[0] - 1] == 0.0) {
      c = rtInf;
    } else {
      c = U_data[0] / U_data[U_size[0] - 1];
    }
  }

  return c;
}

/* End of code generation (cond.c) */
