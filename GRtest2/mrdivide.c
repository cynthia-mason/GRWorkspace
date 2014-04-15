/*
 * mrdivide.c
 *
 * Code generation for function 'mrdivide'
 *
 * C source code generated on: Fri Apr 11 11:48:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen2.h"
#include "mrdivide.h"
#include "Function_PCG_Wood_clean_codegen2_emxutil.h"
#include "colon.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void eml_lusolve(const emxArray_real_T *A, emxArray_real_T *B);
static real_T eml_matlab_zlarfg(int32_T n, real_T *alpha1, emxArray_real_T *x,
  int32_T ix0);
static void eml_qrsolve(const emxArray_real_T *A, emxArray_real_T *B,
  emxArray_real_T *Y);
static real_T eml_xnrm2(int32_T n, const emxArray_real_T *x, int32_T ix0);
static void eml_xscal(int32_T n, real_T a, emxArray_real_T *x, int32_T ix0);
static void eml_xswap(int32_T n, emxArray_real_T *x, int32_T ix0, int32_T incx,
                      int32_T iy0, int32_T incy);
static real_T rt_hypotd_snf(real_T u0, real_T u1);

/* Function Definitions */
static void eml_lusolve(const emxArray_real_T *A, emxArray_real_T *B)
{
  emxArray_real_T *X;
  emxArray_real_T *b_A;
  int32_T n;
  int32_T i13;
  int32_T jy;
  emxArray_int32_T *ipiv;
  int32_T u0;
  int32_T j;
  int32_T mmj;
  int32_T jj;
  int32_T jp1j;
  int32_T c;
  int32_T ix;
  real_T smax;
  int32_T jA;
  real_T s;
  int32_T i;
  emxInit_real_T(&X, 2);
  emxInit_real_T(&b_A, 2);
  n = A->size[1];
  i13 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)b_A, i13, (int32_T)sizeof(real_T));
  jy = A->size[0] * A->size[1] - 1;
  for (i13 = 0; i13 <= jy; i13++) {
    b_A->data[i13] = A->data[i13];
  }

  b_emxInit_int32_T(&ipiv, 2);
  eml_signed_integer_colon(n, ipiv);
  if (n < 1) {
  } else {
    u0 = n - 1;
    if (u0 <= n) {
    } else {
      u0 = n;
    }

    for (j = 1; j <= u0; j++) {
      mmj = n - j;
      jj = (j - 1) * (n + 1);
      jp1j = jj + 2;
      c = mmj + 1;
      if (c < 1) {
        jy = -1;
      } else {
        jy = 0;
        if (c > 1) {
          ix = jj;
          smax = fabs(b_A->data[jj]);
          for (jA = 1; jA + 1 <= c; jA++) {
            ix++;
            s = fabs(b_A->data[ix]);
            if (s > smax) {
              jy = jA;
              smax = s;
            }
          }
        }
      }

      if (b_A->data[jj + jy] != 0.0) {
        if (jy != 0) {
          ipiv->data[j - 1] = j + jy;
          eml_xswap(n, b_A, j, n, j + jy, n);
        }

        i13 = (jp1j + mmj) - 1;
        for (i = jp1j; i <= i13; i++) {
          b_A->data[i - 1] /= b_A->data[jj];
        }
      }

      c = n - j;
      jA = (jj + n) + 1;
      jy = jj + n;
      for (jj = 1; jj <= c; jj++) {
        smax = b_A->data[jy];
        if (b_A->data[jy] != 0.0) {
          ix = jp1j;
          i13 = mmj + jA;
          for (i = jA; i + 1 <= i13; i++) {
            b_A->data[i] += b_A->data[ix - 1] * -smax;
            ix++;
          }
        }

        jy += n;
        jA += n;
      }
    }
  }

  jy = B->size[1];
  for (i = 0; i + 1 <= n; i++) {
    if (ipiv->data[i] != i + 1) {
      for (j = 0; j + 1 <= jy; j++) {
        smax = B->data[i + B->size[0] * j];
        B->data[i + B->size[0] * j] = B->data[(ipiv->data[i] + B->size[0] * j) -
          1];
        B->data[(ipiv->data[i] + B->size[0] * j) - 1] = smax;
      }
    }
  }

  emxFree_int32_T(&ipiv);
  if ((jy == 0) || ((B->size[0] == 0) || (B->size[1] == 0))) {
  } else {
    for (j = 1; j <= jy; j++) {
      c = n * (j - 1) - 1;
      for (jA = 1; jA <= n; jA++) {
        jj = n * (jA - 1);
        if (B->data[jA + c] != 0.0) {
          for (i = jA + 1; i <= n; i++) {
            B->data[i + c] -= B->data[jA + c] * b_A->data[(i + jj) - 1];
          }
        }
      }
    }
  }

  if ((jy == 0) || ((B->size[0] == 0) || (B->size[1] == 0))) {
  } else {
    for (j = 1; j <= jy; j++) {
      c = n * (j - 1) - 1;
      for (jA = n; jA > 0; jA--) {
        jj = n * (jA - 1) - 1;
        if (B->data[jA + c] != 0.0) {
          B->data[jA + c] /= b_A->data[jA + jj];
          for (i = 1; i <= jA - 1; i++) {
            B->data[i + c] -= B->data[jA + c] * b_A->data[i + jj];
          }
        }
      }
    }
  }

  emxFree_real_T(&b_A);
  emxFree_real_T(&X);
}

static real_T eml_matlab_zlarfg(int32_T n, real_T *alpha1, emxArray_real_T *x,
  int32_T ix0)
{
  real_T tau;
  int32_T nm1;
  real_T xnorm;
  int32_T knt;
  real_T d1;
  tau = 0.0;
  if (n <= 0) {
  } else {
    nm1 = n - 1;
    xnorm = eml_xnrm2(nm1, x, ix0);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(fabs(*alpha1), xnorm);
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          eml_xscal(nm1, 9.9792015476736E+291, x, ix0);
          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

        xnorm = eml_xnrm2(nm1, x, ix0);
        xnorm = rt_hypotd_snf(fabs(*alpha1), xnorm);
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        tau = (xnorm - *alpha1) / xnorm;
        d1 = 1.0 / (*alpha1 - xnorm);
        eml_xscal(nm1, d1, x, ix0);
        for (nm1 = 1; nm1 <= knt; nm1++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        tau = (xnorm - *alpha1) / xnorm;
        d1 = 1.0 / (*alpha1 - xnorm);
        eml_xscal(nm1, d1, x, ix0);
        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

static void eml_qrsolve(const emxArray_real_T *A, emxArray_real_T *B,
  emxArray_real_T *Y)
{
  emxArray_real_T *b_A;
  int32_T nb;
  real_T smax;
  real_T s;
  int32_T mn;
  int32_T i11;
  int32_T pvt;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  int32_T m;
  int32_T n;
  int32_T b_mn;
  emxArray_real_T *work;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int32_T k;
  int32_T j;
  int32_T i;
  int32_T i_i;
  int32_T nmi;
  int32_T mmi;
  int32_T mmip1;
  int32_T ix;
  int32_T i_ip1;
  int32_T lastc;
  boolean_T exitg3;
  int32_T exitg2;
  real_T y;
  real_T t;
  boolean_T exitg1;
  uint32_T unnamed_idx_0;
  uint32_T unnamed_idx_1;
  emxInit_real_T(&b_A, 2);
  nb = B->size[1] - 1;
  smax = (real_T)A->size[0];
  s = (real_T)A->size[1];
  if (smax <= s) {
    s = smax;
  }

  mn = (int32_T)s;
  i11 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)b_A, i11, (int32_T)sizeof(real_T));
  pvt = A->size[0] * A->size[1] - 1;
  for (i11 = 0; i11 <= pvt; i11++) {
    b_A->data[i11] = A->data[i11];
  }

  b_emxInit_real_T(&tau, 1);
  b_emxInit_int32_T(&jpvt, 2);
  m = A->size[0];
  n = A->size[1];
  if (m <= n) {
    b_mn = m;
  } else {
    b_mn = n;
  }

  i11 = tau->size[0];
  tau->size[0] = b_mn;
  emxEnsureCapacity((emxArray__common *)tau, i11, (int32_T)sizeof(real_T));
  eml_signed_integer_colon(n, jpvt);
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
  } else {
    b_emxInit_real_T(&work, 1);
    i11 = work->size[0];
    work->size[0] = n;
    emxEnsureCapacity((emxArray__common *)work, i11, (int32_T)sizeof(real_T));
    pvt = n - 1;
    for (i11 = 0; i11 <= pvt; i11++) {
      work->data[i11] = 0.0;
    }

    b_emxInit_real_T(&vn1, 1);
    b_emxInit_real_T(&vn2, 1);
    i11 = vn1->size[0];
    vn1->size[0] = n;
    emxEnsureCapacity((emxArray__common *)vn1, i11, (int32_T)sizeof(real_T));
    i11 = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity((emxArray__common *)vn2, i11, (int32_T)sizeof(real_T));
    k = 1;
    for (j = 0; j + 1 <= n; j++) {
      vn1->data[j] = eml_xnrm2(m, A, k);
      vn2->data[j] = vn1->data[j];
      k += m;
    }

    for (i = 0; i + 1 <= b_mn; i++) {
      i_i = i + i * m;
      nmi = n - i;
      mmi = (m - i) - 1;
      mmip1 = 1 + mmi;
      if (nmi < 1) {
        pvt = -1;
      } else {
        pvt = 0;
        if (nmi > 1) {
          ix = i;
          smax = fabs(vn1->data[i]);
          for (k = 2; k <= nmi; k++) {
            ix++;
            s = fabs(vn1->data[ix]);
            if (s > smax) {
              pvt = k - 1;
              smax = s;
            }
          }
        }
      }

      pvt += i;
      if (pvt + 1 != i + 1) {
        eml_xswap(m, b_A, m * pvt + 1, 1, m * i + 1, 1);
        k = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i];
        jpvt->data[i] = k;
        vn1->data[pvt] = vn1->data[i];
        vn2->data[pvt] = vn2->data[i];
      }

      if (i + 1 < m) {
        s = b_A->data[i_i];
        tau->data[i] = eml_matlab_zlarfg(mmip1, &s, b_A, i_i + 2);
      } else {
        smax = b_A->data[i_i];
        s = b_A->data[i_i];
        b_A->data[i_i] = smax;
        tau->data[i] = 0.0;
      }

      b_A->data[i_i] = s;
      if (i + 1 < n) {
        s = b_A->data[i_i];
        b_A->data[i_i] = 1.0;
        i_ip1 = (i + (i + 1) * m) + 1;
        if (tau->data[i] != 0.0) {
          pvt = i_i + mmip1;
          while ((mmip1 > 0) && (b_A->data[pvt - 1] == 0.0)) {
            mmip1--;
            pvt--;
          }

          lastc = nmi - 1;
          exitg3 = FALSE;
          while ((exitg3 == 0U) && (lastc > 0)) {
            pvt = i_ip1 + (lastc - 1) * m;
            nmi = pvt;
            do {
              exitg2 = 0;
              if (nmi <= (pvt + mmip1) - 1) {
                if (b_A->data[nmi - 1] != 0.0) {
                  exitg2 = 1;
                } else {
                  nmi++;
                }
              } else {
                lastc--;
                exitg2 = 2;
              }
            } while (exitg2 == 0U);

            if (exitg2 == 1U) {
              exitg3 = TRUE;
            }
          }
        } else {
          mmip1 = 0;
          lastc = 0;
        }

        if (mmip1 > 0) {
          if (lastc == 0) {
          } else {
            pvt = lastc - 1;
            for (j = 1; j <= pvt + 1; j++) {
              work->data[j - 1] = 0.0;
            }

            j = 0;
            i11 = i_ip1 + m * pvt;
            pvt = i_ip1;
            while ((m > 0) && (pvt <= i11)) {
              ix = i_i;
              smax = 0.0;
              k = pvt + mmip1;
              for (nmi = pvt; nmi <= k - 1; nmi++) {
                smax += b_A->data[nmi - 1] * b_A->data[ix];
                ix++;
              }

              work->data[j] += smax;
              j++;
              pvt += m;
            }
          }

          if (-tau->data[i] == 0.0) {
          } else {
            pvt = i_ip1 - 1;
            k = 0;
            for (j = 1; j <= lastc; j++) {
              if (work->data[k] != 0.0) {
                smax = work->data[k] * -tau->data[i];
                ix = i_i;
                i11 = mmip1 + pvt;
                for (nmi = pvt; nmi + 1 <= i11; nmi++) {
                  b_A->data[nmi] += b_A->data[ix] * smax;
                  ix++;
                }
              }

              k++;
              pvt += m;
            }
          }
        }

        b_A->data[i_i] = s;
      }

      for (j = i + 1; j + 1 <= n; j++) {
        if (vn1->data[j] != 0.0) {
          s = fabs(b_A->data[i + b_A->size[0] * j]) / vn1->data[j];
          y = s * s;
          s = 1.0 - s * s;
          if (1.0 - y < 0.0) {
            s = 0.0;
          }

          smax = vn1->data[j] / vn2->data[j];
          if (s * (smax * smax) <= 1.4901161193847656E-8) {
            if (i + 1 < m) {
              pvt = (i + m * j) + 1;
              y = 0.0;
              if (mmi < 1) {
              } else if (mmi == 1) {
                y = fabs(b_A->data[pvt]);
              } else {
                smax = 2.2250738585072014E-308;
                k = pvt + mmi;
                while (pvt + 1 <= k) {
                  s = fabs(b_A->data[pvt]);
                  if (s > smax) {
                    t = smax / s;
                    y = 1.0 + y * t * t;
                    smax = s;
                  } else {
                    t = s / smax;
                    y += t * t;
                  }

                  pvt++;
                }

                y = smax * sqrt(y);
              }

              vn1->data[j] = y;
              vn2->data[j] = vn1->data[j];
            } else {
              vn1->data[j] = 0.0;
              vn2->data[j] = 0.0;
            }
          } else {
            vn1->data[j] *= sqrt(s);
          }
        }
      }
    }

    emxFree_real_T(&vn2);
    emxFree_real_T(&vn1);
    emxFree_real_T(&work);
  }

  t = 0.0;
  if (mn > 0) {
    k = 0;
    exitg1 = FALSE;
    while ((exitg1 == 0U) && (k <= mn - 1)) {
      smax = (real_T)A->size[0];
      s = (real_T)A->size[1];
      if (smax >= s) {
        s = smax;
      }

      if (fabs(b_A->data[k + b_A->size[0] * k]) <= s * fabs(b_A->data[0]) *
          2.2204460492503131E-16) {
        exitg1 = TRUE;
      } else {
        t++;
        k++;
      }
    }
  }

  unnamed_idx_0 = (uint32_T)A->size[1];
  unnamed_idx_1 = (uint32_T)B->size[1];
  i11 = Y->size[0] * Y->size[1];
  Y->size[0] = (int32_T)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)Y, i11, (int32_T)sizeof(real_T));
  i11 = Y->size[0] * Y->size[1];
  Y->size[1] = (int32_T)unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)Y, i11, (int32_T)sizeof(real_T));
  pvt = (int32_T)unnamed_idx_0 * (int32_T)unnamed_idx_1 - 1;
  for (i11 = 0; i11 <= pvt; i11++) {
    Y->data[i11] = 0.0;
  }

  for (j = 0; j <= mn - 1; j++) {
    if (tau->data[j] != 0.0) {
      for (k = 0; k <= nb; k++) {
        smax = B->data[j + B->size[0] * k];
        i11 = A->size[0] + (int32_T)(1.0 - ((1.0 + (real_T)j) + 1.0));
        for (i = 0; i <= i11 - 1; i++) {
          unnamed_idx_0 = ((uint32_T)j + (uint32_T)i) + 2U;
          smax += b_A->data[((int32_T)unnamed_idx_0 + b_A->size[0] * j) - 1] *
            B->data[((int32_T)unnamed_idx_0 + B->size[0] * k) - 1];
        }

        smax *= tau->data[j];
        if (smax != 0.0) {
          B->data[j + B->size[0] * k] -= smax;
          i11 = A->size[0] + (int32_T)(1.0 - ((1.0 + (real_T)j) + 1.0));
          for (i = 0; i <= i11 - 1; i++) {
            unnamed_idx_0 = ((uint32_T)j + (uint32_T)i) + 2U;
            B->data[((int32_T)unnamed_idx_0 + B->size[0] * k) - 1] -= b_A->data
              [((int32_T)unnamed_idx_0 + b_A->size[0] * j) - 1] * smax;
          }
        }
      }
    }
  }

  emxFree_real_T(&tau);
  for (k = 0; k <= nb; k++) {
    for (i = 0; i <= (int32_T)t - 1; i++) {
      Y->data[(jpvt->data[(int32_T)(1.0 + (real_T)i) - 1] + Y->size[0] * k) - 1]
        = B->data[((int32_T)(1.0 + (real_T)i) + B->size[0] * k) - 1];
    }

    for (j = 0; j <= (int32_T)-(1.0 + (-1.0 - t)) - 1; j++) {
      smax = t + -(real_T)j;
      Y->data[(jpvt->data[(int32_T)smax - 1] + Y->size[0] * k) - 1] /= b_A->
        data[((int32_T)smax + b_A->size[0] * ((int32_T)smax - 1)) - 1];
      for (i = 0; i <= (int32_T)smax - 2; i++) {
        Y->data[(jpvt->data[i] + Y->size[0] * k) - 1] -= Y->data[(jpvt->data
          [(int32_T)smax - 1] + Y->size[0] * k) - 1] * b_A->data[i + b_A->size[0]
          * ((int32_T)smax - 1)];
      }
    }
  }

  emxFree_int32_T(&jpvt);
  emxFree_real_T(&b_A);
}

static real_T eml_xnrm2(int32_T n, const emxArray_real_T *x, int32_T ix0)
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

static void eml_xscal(int32_T n, real_T a, emxArray_real_T *x, int32_T ix0)
{
  emxArray_real_T *b_x;
  int32_T i14;
  int32_T k;
  emxInit_real_T(&b_x, 2);
  i14 = (ix0 + n) - 1;
  for (k = ix0; k <= i14; k++) {
    x->data[k - 1] *= a;
  }

  emxFree_real_T(&b_x);
}

static void eml_xswap(int32_T n, emxArray_real_T *x, int32_T ix0, int32_T incx,
                      int32_T iy0, int32_T incy)
{
  emxArray_real_T *b_x;
  int32_T ix;
  int32_T iy;
  int32_T k;
  real_T temp;
  emxInit_real_T(&b_x, 2);
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 1; k <= n; k++) {
    temp = x->data[ix];
    x->data[ix] = x->data[iy];
    x->data[iy] = temp;
    ix += incx;
    iy += incy;
  }

  emxFree_real_T(&b_x);
}

static real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = a * sqrt(y * y + 1.0);
  } else if (rtIsNaN(y)) {
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

void mrdivide(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *y)
{
  emxArray_real_T *b_A;
  int32_T i9;
  int32_T loop_ub;
  int32_T b_loop_ub;
  int32_T i10;
  emxArray_real_T *b_B;
  emxArray_real_T *c_B;
  uint32_T unnamed_idx_0;
  uint32_T unnamed_idx_1;
  emxInit_real_T(&b_A, 2);
  i9 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = B->size[1];
  b_A->size[1] = B->size[0];
  emxEnsureCapacity((emxArray__common *)b_A, i9, (int32_T)sizeof(real_T));
  loop_ub = B->size[0] - 1;
  for (i9 = 0; i9 <= loop_ub; i9++) {
    b_loop_ub = B->size[1] - 1;
    for (i10 = 0; i10 <= b_loop_ub; i10++) {
      b_A->data[i10 + b_A->size[0] * i9] = B->data[i9 + B->size[0] * i10];
    }
  }

  emxInit_real_T(&b_B, 2);
  i9 = b_B->size[0] * b_B->size[1];
  b_B->size[0] = A->size[1];
  b_B->size[1] = A->size[0];
  emxEnsureCapacity((emxArray__common *)b_B, i9, (int32_T)sizeof(real_T));
  loop_ub = A->size[0] - 1;
  for (i9 = 0; i9 <= loop_ub; i9++) {
    b_loop_ub = A->size[1] - 1;
    for (i10 = 0; i10 <= b_loop_ub; i10++) {
      b_B->data[i10 + b_B->size[0] * i9] = A->data[i9 + A->size[0] * i10];
    }
  }

  emxInit_real_T(&c_B, 2);
  if ((b_A->size[0] == 0) || (b_A->size[1] == 0) || ((b_B->size[0] == 0) ||
       (b_B->size[1] == 0))) {
    unnamed_idx_0 = (uint32_T)b_A->size[1];
    unnamed_idx_1 = (uint32_T)b_B->size[1];
    i9 = b_B->size[0] * b_B->size[1];
    b_B->size[0] = (int32_T)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)b_B, i9, (int32_T)sizeof(real_T));
    i9 = b_B->size[0] * b_B->size[1];
    b_B->size[1] = (int32_T)unnamed_idx_1;
    emxEnsureCapacity((emxArray__common *)b_B, i9, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)unnamed_idx_0 * (int32_T)unnamed_idx_1 - 1;
    for (i9 = 0; i9 <= loop_ub; i9++) {
      b_B->data[i9] = 0.0;
    }
  } else if (b_A->size[0] == b_A->size[1]) {
    eml_lusolve(b_A, b_B);
  } else {
    i9 = c_B->size[0] * c_B->size[1];
    c_B->size[0] = b_B->size[0];
    c_B->size[1] = b_B->size[1];
    emxEnsureCapacity((emxArray__common *)c_B, i9, (int32_T)sizeof(real_T));
    loop_ub = b_B->size[0] * b_B->size[1] - 1;
    for (i9 = 0; i9 <= loop_ub; i9++) {
      c_B->data[i9] = b_B->data[i9];
    }

    eml_qrsolve(b_A, c_B, b_B);
  }

  emxFree_real_T(&c_B);
  emxFree_real_T(&b_A);
  i9 = y->size[0] * y->size[1];
  y->size[0] = b_B->size[1];
  y->size[1] = b_B->size[0];
  emxEnsureCapacity((emxArray__common *)y, i9, (int32_T)sizeof(real_T));
  loop_ub = b_B->size[0] - 1;
  for (i9 = 0; i9 <= loop_ub; i9++) {
    b_loop_ub = b_B->size[1] - 1;
    for (i10 = 0; i10 <= b_loop_ub; i10++) {
      y->data[i10 + y->size[0] * i9] = b_B->data[i9 + b_B->size[0] * i10];
    }
  }

  emxFree_real_T(&b_B);
}

/* End of code generation (mrdivide.c) */
