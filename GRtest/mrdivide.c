/*
 * mrdivide.c
 *
 * Code generation for function 'mrdivide'
 *
 * C source code generated on: Mon Apr 07 11:43:50 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Function_PCG_Wood_clean_codegen.h"
#include "mrdivide.h"
#include "colon.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static real_T b_eml_matlab_zlarfg(real_T *alpha1, const real_T *x);
static real_T b_eml_xnrm2(int32_T n, const real_T x_data[289], const int32_T
  x_size[2], int32_T ix0);
static int32_T eml_ixamax(int32_T n, const real_T x_data[17], const int32_T
  x_size[1], int32_T ix0);
static void eml_lusolve(const real_T A_data[289], const int32_T A_size[2],
  real_T B_data[289], int32_T B_size[2]);
static void eml_matlab_zlarf(int32_T m, int32_T n, int32_T iv0, real_T tau,
  real_T C_data[289], int32_T C_size[2], int32_T ic0, int32_T ldc, real_T
  work_data[17], int32_T work_size[1]);
static real_T eml_matlab_zlarfg(int32_T n, real_T *alpha1, real_T x_data[289],
  int32_T x_size[2], int32_T ix0);
static void eml_qrsolve(const real_T A_data[289], const int32_T A_size[2],
  real_T B_data[289], int32_T B_size[2], real_T Y_data[289], int32_T Y_size[2]);
static real_T eml_xnrm2(int32_T n, const real_T x_data[289], const int32_T
  x_size[2], int32_T ix0);
static void eml_xswap(int32_T n, real_T x_data[289], int32_T x_size[2], int32_T
                      ix0, int32_T incx, int32_T iy0, int32_T incy);
static real_T rt_hypotd_snf(real_T u0, real_T u1);

/* Function Definitions */
static real_T b_eml_matlab_zlarfg(real_T *alpha1, const real_T *x)
{
  return 0.0;
}

static real_T b_eml_xnrm2(int32_T n, const real_T x_data[289], const int32_T
  x_size[2], int32_T ix0)
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
    y = fabs(x_data[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x_data[k - 1]);
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

static int32_T eml_ixamax(int32_T n, const real_T x_data[17], const int32_T
  x_size[1], int32_T ix0)
{
  int32_T idxmax;
  int32_T ix;
  real_T smax;
  int32_T k;
  real_T s;
  if (n < 1) {
    idxmax = 0;
  } else {
    idxmax = 1;
    if (n > 1) {
      ix = ix0 - 1;
      smax = fabs(x_data[ix0 - 1]);
      for (k = 2; k <= n; k++) {
        ix++;
        s = fabs(x_data[ix]);
        if (s > smax) {
          idxmax = k;
          smax = s;
        }
      }
    }
  }

  return idxmax;
}

static void eml_lusolve(const real_T A_data[289], const int32_T A_size[2],
  real_T B_data[289], int32_T B_size[2])
{
  int32_T b_A_size[2];
  int32_T jy;
  int32_T i13;
  real_T b_A_data[289];
  int32_T u1;
  int32_T ipiv_size[2];
  int32_T ipiv_data[17];
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
  b_A_size[0] = A_size[0];
  b_A_size[1] = A_size[1];
  jy = A_size[0] * A_size[1] - 1;
  for (i13 = 0; i13 <= jy; i13++) {
    b_A_data[i13] = A_data[i13];
  }

  jy = A_size[1];
  u1 = A_size[1];
  if (jy <= u1) {
    u1 = jy;
  }

  eml_signed_integer_colon(u1, ipiv_data, ipiv_size);
  if (A_size[1] < 1) {
  } else {
    jy = A_size[1] - 1;
    u1 = A_size[1];
    if (jy <= u1) {
      u1 = jy;
    }

    for (j = 1; j <= u1; j++) {
      mmj = A_size[1] - j;
      jj = (j - 1) * (A_size[1] + 1);
      jp1j = jj + 2;
      c = mmj + 1;
      if (c < 1) {
        jy = -1;
      } else {
        jy = 0;
        if (c > 1) {
          ix = jj;
          smax = fabs(b_A_data[jj]);
          for (jA = 1; jA + 1 <= c; jA++) {
            ix++;
            s = fabs(b_A_data[ix]);
            if (s > smax) {
              jy = jA;
              smax = s;
            }
          }
        }
      }

      if (b_A_data[jj + jy] != 0.0) {
        if (jy != 0) {
          ipiv_data[j - 1] = j + jy;
          eml_xswap(A_size[1], b_A_data, b_A_size, j, A_size[1], j + jy, A_size
                    [1]);
        }

        i13 = (jp1j + mmj) - 1;
        for (i = jp1j; i <= i13; i++) {
          b_A_data[i - 1] /= b_A_data[jj];
        }
      }

      c = A_size[1] - j;
      jA = (jj + A_size[1]) + 1;
      jy = jj + A_size[1];
      for (jj = 1; jj <= c; jj++) {
        smax = b_A_data[jy];
        if (b_A_data[jy] != 0.0) {
          ix = jp1j;
          i13 = mmj + jA;
          for (i = jA; i + 1 <= i13; i++) {
            b_A_data[i] += b_A_data[ix - 1] * -smax;
            ix++;
          }
        }

        jy += A_size[1];
        jA += A_size[1];
      }
    }
  }

  jy = B_size[1];
  for (i = 0; i + 1 <= A_size[1]; i++) {
    if (ipiv_data[i] != i + 1) {
      for (j = 0; j + 1 <= jy; j++) {
        smax = B_data[i + B_size[0] * j];
        B_data[i + B_size[0] * j] = B_data[(ipiv_data[i] + B_size[0] * j) - 1];
        B_data[(ipiv_data[i] + B_size[0] * j) - 1] = smax;
      }
    }
  }

  if ((jy == 0) || ((B_size[0] == 0) || (B_size[1] == 0))) {
  } else {
    for (j = 1; j <= jy; j++) {
      c = A_size[1] * (j - 1) - 1;
      for (jA = 1; jA <= A_size[1]; jA++) {
        jj = A_size[1] * (jA - 1);
        if (B_data[jA + c] != 0.0) {
          for (i = jA + 1; i <= A_size[1]; i++) {
            B_data[i + c] -= B_data[jA + c] * b_A_data[(i + jj) - 1];
          }
        }
      }
    }
  }

  if ((jy == 0) || ((B_size[0] == 0) || (B_size[1] == 0))) {
  } else {
    for (j = 1; j <= jy; j++) {
      c = A_size[1] * (j - 1) - 1;
      for (jA = A_size[1]; jA > 0; jA--) {
        jj = A_size[1] * (jA - 1) - 1;
        if (B_data[jA + c] != 0.0) {
          B_data[jA + c] /= b_A_data[jA + jj];
          for (i = 1; i <= jA - 1; i++) {
            B_data[i + c] -= B_data[jA + c] * b_A_data[i + jj];
          }
        }
      }
    }
  }
}

static void eml_matlab_zlarf(int32_T m, int32_T n, int32_T iv0, real_T tau,
  real_T C_data[289], int32_T C_size[2], int32_T ic0, int32_T ldc, real_T
  work_data[17], int32_T work_size[1])
{
  int32_T lastv;
  int32_T i;
  int32_T lastc;
  boolean_T exitg2;
  int32_T ia;
  int32_T exitg1;
  int32_T iy;
  int32_T i15;
  int32_T ix;
  real_T c;
  int32_T jy;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C_data[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n;
    exitg2 = FALSE;
    while ((exitg2 == 0U) && (lastc > 0)) {
      i = ic0 + (lastc - 1) * ldc;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C_data[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0U);

      if (exitg1 == 1U) {
        exitg2 = TRUE;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc == 0) {
    } else {
      i = lastc - 1;
      for (iy = 1; iy <= i + 1; iy++) {
        work_data[iy - 1] = 0.0;
      }

      iy = 0;
      i15 = ic0 + ldc * i;
      i = ic0;
      while ((ldc > 0) && (i <= i15)) {
        ix = iv0;
        c = 0.0;
        jy = i + lastv;
        for (ia = i; ia <= jy - 1; ia++) {
          c += C_data[ia - 1] * C_data[ix - 1];
          ix++;
        }

        work_data[iy] += c;
        iy++;
        i += ldc;
      }
    }

    if (-tau == 0.0) {
    } else {
      i = ic0 - 1;
      jy = 0;
      for (ia = 1; ia <= lastc; ia++) {
        if (work_data[jy] != 0.0) {
          c = work_data[jy] * -tau;
          ix = iv0;
          i15 = lastv + i;
          for (iy = i; iy + 1 <= i15; iy++) {
            C_data[iy] += C_data[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += ldc;
      }
    }
  }
}

static real_T eml_matlab_zlarfg(int32_T n, real_T *alpha1, real_T x_data[289],
  int32_T x_size[2], int32_T ix0)
{
  real_T tau;
  int32_T nm1;
  real_T xnorm;
  int32_T knt;
  int32_T i14;
  int32_T k;
  tau = 0.0;
  if (n <= 0) {
  } else {
    nm1 = n - 1;
    xnorm = eml_xnrm2(nm1, x_data, x_size, ix0);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(fabs(*alpha1), xnorm);
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          i14 = (ix0 + nm1) - 1;
          for (k = ix0; k <= i14; k++) {
            x_data[k - 1] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

        xnorm = eml_xnrm2(nm1, x_data, x_size, ix0);
        xnorm = rt_hypotd_snf(fabs(*alpha1), xnorm);
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        tau = (xnorm - *alpha1) / xnorm;
        *alpha1 = 1.0 / (*alpha1 - xnorm);
        i14 = (ix0 + nm1) - 1;
        for (k = ix0; k <= i14; k++) {
          x_data[k - 1] *= *alpha1;
        }

        for (k = 1; k <= knt; k++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        tau = (xnorm - *alpha1) / xnorm;
        *alpha1 = 1.0 / (*alpha1 - xnorm);
        i14 = (ix0 + nm1) - 1;
        for (k = ix0; k <= i14; k++) {
          x_data[k - 1] *= *alpha1;
        }

        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

static void eml_qrsolve(const real_T A_data[289], const int32_T A_size[2],
  real_T B_data[289], int32_T B_size[2], real_T Y_data[289], int32_T Y_size[2])
{
  real_T atmp;
  real_T temp1;
  int32_T mn;
  int32_T b_A_size[2];
  int32_T itemp;
  int32_T pvt;
  real_T b_A_data[289];
  int32_T b_mn;
  real_T tau_data[17];
  int32_T jpvt_size[2];
  int32_T jpvt_data[17];
  int32_T work_size[1];
  real_T work_data[17];
  real_T vn1_data[17];
  int32_T vn1_size[1];
  real_T vn2_data[17];
  int32_T k;
  int32_T i_i;
  int32_T i;
  int32_T nmi;
  int32_T mmi;
  int32_T mmip1;
  real_T rankR;
  boolean_T exitg1;
  atmp = (real_T)A_size[0];
  temp1 = (real_T)A_size[1];
  if (atmp <= temp1) {
    temp1 = atmp;
  }

  mn = (int32_T)temp1;
  b_A_size[0] = A_size[0];
  b_A_size[1] = A_size[1];
  itemp = A_size[0] * A_size[1] - 1;
  for (pvt = 0; pvt <= itemp; pvt++) {
    b_A_data[pvt] = A_data[pvt];
  }

  itemp = A_size[0];
  b_mn = A_size[1];
  if (itemp <= b_mn) {
    b_mn = itemp;
  }

  eml_signed_integer_colon(A_size[1], jpvt_data, jpvt_size);
  if ((A_size[0] == 0) || (A_size[1] == 0)) {
  } else {
    work_size[0] = (int8_T)A_size[1];
    itemp = (int8_T)A_size[1] - 1;
    for (pvt = 0; pvt <= itemp; pvt++) {
      work_data[pvt] = 0.0;
    }

    vn1_size[0] = (int8_T)A_size[1];
    k = 1;
    for (i_i = 0; i_i + 1 <= A_size[1]; i_i++) {
      vn1_data[i_i] = eml_xnrm2(A_size[0], A_data, A_size, k);
      vn2_data[i_i] = vn1_data[i_i];
      k += A_size[0];
    }

    for (i = 1; i <= b_mn; i++) {
      k = i - 1;
      i_i = (i + k * A_size[0]) - 1;
      nmi = A_size[1] - i;
      mmi = A_size[0] - i;
      mmip1 = 1 + mmi;
      itemp = eml_ixamax(1 + nmi, vn1_data, vn1_size, i) - 1;
      pvt = k + itemp;
      if (pvt + 1 != i) {
        eml_xswap(A_size[0], b_A_data, b_A_size, A_size[0] * pvt + 1, 1, A_size
                  [0] * k + 1, 1);
        itemp = jpvt_data[pvt];
        jpvt_data[pvt] = jpvt_data[i - 1];
        jpvt_data[i - 1] = itemp;
        vn1_data[pvt] = vn1_data[i - 1];
        vn2_data[pvt] = vn2_data[i - 1];
      }

      if (i < A_size[0]) {
        atmp = b_A_data[i_i];
        tau_data[i - 1] = eml_matlab_zlarfg(mmip1, &atmp, b_A_data, b_A_size,
          i_i + 2);
      } else {
        atmp = b_A_data[i_i];
        tau_data[i - 1] = b_eml_matlab_zlarfg(&atmp, &b_A_data[i_i]);
      }

      b_A_data[i_i] = atmp;
      if (i < A_size[1]) {
        atmp = b_A_data[i_i];
        b_A_data[i_i] = 1.0;
        eml_matlab_zlarf(mmip1, nmi, i_i + 1, tau_data[i - 1], b_A_data,
                         b_A_size, i + i * A_size[0], A_size[0], work_data,
                         work_size);
        b_A_data[i_i] = atmp;
      }

      for (i_i = i; i_i + 1 <= A_size[1]; i_i++) {
        if (vn1_data[i_i] != 0.0) {
          temp1 = fabs(b_A_data[(i + b_A_size[0] * i_i) - 1]) / vn1_data[i_i];
          atmp = temp1 * temp1;
          temp1 = 1.0 - temp1 * temp1;
          if (1.0 - atmp < 0.0) {
            temp1 = 0.0;
          }

          atmp = vn1_data[i_i] / vn2_data[i_i];
          if (temp1 * (atmp * atmp) <= 1.4901161193847656E-8) {
            if (i < A_size[0]) {
              vn1_data[i_i] = b_eml_xnrm2(mmi, b_A_data, b_A_size, (i + A_size[0]
                * i_i) + 1);
              vn2_data[i_i] = vn1_data[i_i];
            } else {
              vn1_data[i_i] = 0.0;
              vn2_data[i_i] = 0.0;
            }
          } else {
            vn1_data[i_i] *= sqrt(temp1);
          }
        }
      }
    }
  }

  rankR = 0.0;
  if (mn > 0) {
    k = 0;
    exitg1 = FALSE;
    while ((exitg1 == 0U) && (k <= mn - 1)) {
      atmp = (real_T)A_size[0];
      temp1 = (real_T)A_size[1];
      if (atmp >= temp1) {
        temp1 = atmp;
      }

      if (fabs(b_A_data[k + b_A_size[0] * k]) <= temp1 * fabs(b_A_data[0]) *
          2.2204460492503131E-16) {
        exitg1 = TRUE;
      } else {
        rankR++;
        k++;
      }
    }
  }

  Y_size[0] = (int8_T)A_size[1];
  Y_size[1] = (int8_T)B_size[1];
  itemp = (int8_T)A_size[1] * (int8_T)B_size[1] - 1;
  for (pvt = 0; pvt <= itemp; pvt++) {
    Y_data[pvt] = 0.0;
  }

  for (i_i = 0; i_i <= mn - 1; i_i++) {
    if (tau_data[i_i] != 0.0) {
      for (k = 0; k <= B_size[1] - 1; k++) {
        atmp = B_data[i_i + B_size[0] * k];
        pvt = A_size[0] - i_i;
        for (i = 0; i <= pvt - 2; i++) {
          itemp = (i_i + i) + 1;
          atmp += b_A_data[itemp + b_A_size[0] * i_i] * B_data[itemp + B_size[0]
            * k];
        }

        atmp *= tau_data[i_i];
        if (atmp != 0.0) {
          B_data[i_i + B_size[0] * k] -= atmp;
          pvt = A_size[0] - i_i;
          for (i = 0; i <= pvt - 2; i++) {
            itemp = (i_i + i) + 1;
            B_data[itemp + B_size[0] * k] -= b_A_data[itemp + b_A_size[0] * i_i]
              * atmp;
          }
        }
      }
    }
  }

  for (k = 0; k <= B_size[1] - 1; k++) {
    for (i = 0; i <= (int32_T)rankR - 1; i++) {
      Y_data[(jpvt_data[(int32_T)(1.0 + (real_T)i) - 1] + Y_size[0] * k) - 1] =
        B_data[((int32_T)(1.0 + (real_T)i) + B_size[0] * k) - 1];
    }

    for (i_i = 0; i_i <= (int32_T)-(1.0 + (-1.0 - rankR)) - 1; i_i++) {
      atmp = rankR + -(real_T)i_i;
      Y_data[(jpvt_data[(int32_T)atmp - 1] + Y_size[0] * k) - 1] /= b_A_data
        [((int32_T)atmp + b_A_size[0] * ((int32_T)atmp - 1)) - 1];
      for (i = 0; i <= (int32_T)atmp - 2; i++) {
        Y_data[(jpvt_data[i] + Y_size[0] * k) - 1] -= Y_data[(jpvt_data[(int32_T)
          atmp - 1] + Y_size[0] * k) - 1] * b_A_data[i + b_A_size[0] * ((int32_T)
          atmp - 1)];
      }
    }
  }
}

static real_T eml_xnrm2(int32_T n, const real_T x_data[289], const int32_T
  x_size[2], int32_T ix0)
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
    y = fabs(x_data[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x_data[k - 1]);
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

static void eml_xswap(int32_T n, real_T x_data[289], int32_T x_size[2], int32_T
                      ix0, int32_T incx, int32_T iy0, int32_T incy)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  real_T temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 1; k <= n; k++) {
    temp = x_data[ix];
    x_data[ix] = x_data[iy];
    x_data[iy] = temp;
    ix += incx;
    iy += incy;
  }
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

void mrdivide(const real_T A_data[289], const int32_T A_size[2], const real_T
              B_data[289], const int32_T B_size[2], real_T y_data[289], int32_T
              y_size[2])
{
  int32_T b_A_size[2];
  int32_T loop_ub;
  int32_T i10;
  int32_T b_loop_ub;
  int32_T i11;
  real_T b_A_data[289];
  int32_T b_B_size[2];
  real_T b_B_data[289];
  int8_T unnamed_idx_1;
  real_T c_B_data[289];
  int32_T c_B_size[2];
  b_A_size[0] = B_size[1];
  b_A_size[1] = B_size[0];
  loop_ub = B_size[0] - 1;
  for (i10 = 0; i10 <= loop_ub; i10++) {
    b_loop_ub = B_size[1] - 1;
    for (i11 = 0; i11 <= b_loop_ub; i11++) {
      b_A_data[i11 + b_A_size[0] * i10] = B_data[i10 + B_size[0] * i11];
    }
  }

  b_B_size[0] = A_size[1];
  b_B_size[1] = A_size[0];
  loop_ub = A_size[0] - 1;
  for (i10 = 0; i10 <= loop_ub; i10++) {
    b_loop_ub = A_size[1] - 1;
    for (i11 = 0; i11 <= b_loop_ub; i11++) {
      b_B_data[i11 + b_B_size[0] * i10] = A_data[i10 + A_size[0] * i11];
    }
  }

  if ((b_A_size[0] == 0) || (b_A_size[1] == 0) || ((b_B_size[0] == 0) ||
       (b_B_size[1] == 0))) {
    unnamed_idx_1 = (int8_T)b_B_size[1];
    b_B_size[0] = (int8_T)b_A_size[1];
    b_B_size[1] = unnamed_idx_1;
    loop_ub = (int8_T)b_A_size[1] * unnamed_idx_1 - 1;
    for (i10 = 0; i10 <= loop_ub; i10++) {
      b_B_data[i10] = 0.0;
    }
  } else if (b_A_size[0] == b_A_size[1]) {
    eml_lusolve(b_A_data, b_A_size, b_B_data, b_B_size);
  } else {
    c_B_size[0] = b_B_size[0];
    c_B_size[1] = b_B_size[1];
    loop_ub = b_B_size[0] * b_B_size[1] - 1;
    for (i10 = 0; i10 <= loop_ub; i10++) {
      c_B_data[i10] = b_B_data[i10];
    }

    eml_qrsolve(b_A_data, b_A_size, c_B_data, c_B_size, b_B_data, b_B_size);
  }

  y_size[0] = b_B_size[1];
  y_size[1] = b_B_size[0];
  loop_ub = b_B_size[0] - 1;
  for (i10 = 0; i10 <= loop_ub; i10++) {
    b_loop_ub = b_B_size[1] - 1;
    for (i11 = 0; i11 <= b_loop_ub; i11++) {
      y_data[i11 + y_size[0] * i10] = b_B_data[i10 + b_B_size[0] * i11];
    }
  }
}

/* End of code generation (mrdivide.c) */
