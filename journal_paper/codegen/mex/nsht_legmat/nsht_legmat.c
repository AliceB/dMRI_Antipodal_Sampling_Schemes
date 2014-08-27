/*
 * nsht_legmat.c
 *
 * Code generation for function 'nsht_legmat'
 *
 * C source code generated on: Mon Aug 11 15:44:13 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nsht_legmat.h"
#include "nsht_legmat_emxutil.h"
#include "nsht_legmat_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 19, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRSInfo b_emlrtRSI = { 31, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRSInfo c_emlrtRSI = { 39, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRSInfo d_emlrtRSI = { 41, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRSInfo e_emlrtRSI = { 66, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRSInfo f_emlrtRSI = { 16, "error",
  "C:/Program Files/matlab/r2013b/toolbox/eml/lib/matlab/lang/error.m" };

static emlrtRSInfo g_emlrtRSI = { 102, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRSInfo h_emlrtRSI = { 14, "sqrt",
  "C:/Program Files/matlab/r2013b/toolbox/eml/lib/matlab/elfun/sqrt.m" };

static emlrtRSInfo i_emlrtRSI = { 20, "eml_error",
  "C:/Program Files/matlab/r2013b/toolbox/eml/lib/matlab/eml/eml_error.m" };

static emlrtRSInfo j_emlrtRSI = { 37, "mpower",
  "C:/Program Files/matlab/r2013b/toolbox/eml/lib/matlab/ops/mpower.m" };

static emlrtRSInfo k_emlrtRSI = { 42, "power",
  "C:/Program Files/matlab/r2013b/toolbox/eml/lib/matlab/ops/power.m" };

static emlrtRSInfo l_emlrtRSI = { 56, "power",
  "C:/Program Files/matlab/r2013b/toolbox/eml/lib/matlab/ops/power.m" };

static emlrtRSInfo m_emlrtRSI = { 85, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRSInfo n_emlrtRSI = { 87, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtMCInfo emlrtMCI = { 16, 1, "error",
  "C:/Program Files/matlab/r2013b/toolbox/eml/lib/matlab/lang/error.m" };

static emlrtRTEInfo emlrtRTEI = { 1, 19, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtBCInfo emlrtBCI = { -1, -1, 37, 10, "thetas", "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 0 };

static emlrtRTEInfo c_emlrtRTEI = { 45, 5, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRTEInfo d_emlrtRTEI = { 100, 1, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtRTEInfo e_emlrtRTEI = { 126, 5, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m" };

static emlrtDCInfo emlrtDCI = { 26, 11, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 1 };

static emlrtDCInfo b_emlrtDCI = { 26, 11, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 4 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 61, 11, "P", "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 0 };

static emlrtDCInfo c_emlrtDCI = { 61, 11, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 1 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 61, 19, "P", "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 0 };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 63, 12, "Sc", "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 0 };

static emlrtDCInfo d_emlrtDCI = { 63, 12, "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 1 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 63, 20, "Sc", "nsht_legmat",
  "C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m", 0 };

static emlrtRTEInfo f_emlrtRTEI = { 20, 5, "eml_error",
  "C:/Program Files/matlab/r2013b/toolbox/eml/lib/matlab/eml/eml_error.m" };

/* Function Declarations */
static void b_eml_error(void);
static void eml_error(void);
static void error(const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */
static void b_eml_error(void)
{
  emlrtPushRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
  emlrtErrorWithMessageIdR2012b(emlrtRootTLSGlobal, &f_emlrtRTEI,
    "Coder:toolbox:power_domainError", 0);
  emlrtPopRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
}

static void eml_error(void)
{
  static char_T cv2[4][1] = { { 's' }, { 'q' }, { 'r' }, { 't' } };

  emlrtPushRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
  emlrtErrorWithMessageIdR2012b(emlrtRootTLSGlobal, &f_emlrtRTEI,
    "Coder:toolbox:ElFunDomainError", 3, 4, 4, cv2);
  emlrtPopRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
}

static void error(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 0, NULL, 1, &pArray, "error", TRUE,
                        location);
}

void nsht_legmat(const emxArray_real_T *thetas, real_T L, real_T m,
                 emxArray_real_T *P, emxArray_real_T *Sc)
{
  const mxArray *y;
  static const int32_T iv0[2] = { 1, 14 };

  const mxArray *m0;
  char_T cv0[14];
  int32_T i;
  static const char_T cv1[14] = { 'R', 'e', 'q', 'u', 'i', 'r', 'e', ' ', 'm',
    ' ', '<', ' ', 'L', '.' };

  real_T factor_remain2;
  int32_T i0;
  uint32_T unnamed_idx_0;
  uint32_T unnamed_idx_1;
  real_T Km;
  real_T power_10_first;
  int32_T ii;
  real_T b_y;
  real_T dlm1;
  real_T temp;
  real_T dlm;
  real_T dlm_1;
  int32_T ell;
  real_T b_ell;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  real_T c_y;
  real_T d_y;

  /*  nsht_legmat - Compute Legendre matrix */
  /*  */
  /*  Compute the Legendre matrix.  Usage is given by */
  /*  */
  /*   [P Sc] = nsht_mat(thetas, L, m) */
  /*  */
  /*  where thetas is the vector of theta samples, L is the harmonic band-limit */
  /*  and m is the order considered.  The computed matrix P is ordered  */
  /*  P(ell, theta) where ell \in [m, L-1] and theta \in thetas. */
  /*  */
  /*  Notes: */
  /*   - Kostelec recusrsion is implemented to compute the scaled legendre */
  /*   coefficients for each theta \in thetas */
  /*  Check arguments.  */
  if (m >= L) {
    emlrtPushRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
    y = NULL;
    m0 = mxCreateCharArray(2, iv0);
    for (i = 0; i < 14; i++) {
      cv0[i] = cv1[i];
    }

    emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 14, m0, cv0);
    emlrtAssign(&y, m0);
    error(y, &emlrtMCI);
    emlrtPopRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPopRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
  }

  /* if ( length(thetas) ~= L)  */
  /*   error('Require L theta samples.'); */
  /* end */
  factor_remain2 = L - m;
  factor_remain2 = emlrtNonNegativeCheckFastR2012b(factor_remain2, &b_emlrtDCI,
    emlrtRootTLSGlobal);
  i = (int32_T)emlrtIntegerCheckFastR2012b(factor_remain2, &emlrtDCI,
    emlrtRootTLSGlobal);
  i0 = P->size[0] * P->size[1];
  P->size[0] = i;
  emxEnsureCapacity((emxArray__common *)P, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  factor_remain2 = L - m;
  factor_remain2 = emlrtNonNegativeCheckFastR2012b(factor_remain2, &b_emlrtDCI,
    emlrtRootTLSGlobal);
  emlrtIntegerCheckFastR2012b(factor_remain2, &emlrtDCI, emlrtRootTLSGlobal);
  i = thetas->size[1];
  i0 = P->size[0] * P->size[1];
  P->size[1] = i;
  emxEnsureCapacity((emxArray__common *)P, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  factor_remain2 = L - m;
  factor_remain2 = emlrtNonNegativeCheckFastR2012b(factor_remain2, &b_emlrtDCI,
    emlrtRootTLSGlobal);
  i = (int32_T)emlrtIntegerCheckFastR2012b(factor_remain2, &emlrtDCI,
    emlrtRootTLSGlobal) * thetas->size[1];
  for (i0 = 0; i0 < i; i0++) {
    P->data[i0] = 0.0;
  }

  unnamed_idx_0 = (uint32_T)(int32_T)(L - m);
  unnamed_idx_1 = (uint32_T)thetas->size[1];
  i0 = Sc->size[0] * Sc->size[1];
  Sc->size[0] = (int32_T)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)Sc, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  i0 = Sc->size[0] * Sc->size[1];
  Sc->size[1] = (int32_T)unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)Sc, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  i = (int32_T)unnamed_idx_0 * (int32_T)unnamed_idx_1;
  for (i0 = 0; i0 < i; i0++) {
    Sc->data[i0] = 0.0;
  }

  /*  scaling matrix contains exponents of 10 */
  emlrtPushRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);
  Km = 1.0;
  power_10_first = 0.0;
  emlrtForLoopVectorCheckR2012b(1.0, 1.0, m, mxDOUBLE_CLASS, (int32_T)m,
    &d_emlrtRTEI, emlrtRootTLSGlobal);
  ii = 0;
  while (ii <= (int32_T)m - 1) {
    emlrtPushRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
    b_y = ((2.0 * m - (1.0 + (real_T)ii)) + 1.0) / ((m - (1.0 + (real_T)ii)) +
      1.0);
    if (b_y < 0.0) {
      emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
      eml_error();
      emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
    }

    Km = Km * muDoubleScalarSqrt(b_y) / 2.0;
    emlrtPopRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
    if (Km > 100.0) {
      while (Km > 100.0) {
        power_10_first += 2.0;
        Km *= 0.01;
        emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar,
          emlrtRootTLSGlobal);
      }
    }

    if (Km < 0.01) {
      while (Km < 0.01) {
        power_10_first -= 2.0;
        Km *= 100.0;
        emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar,
          emlrtRootTLSGlobal);
      }
    }

    ii++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
  }

  emlrtPopRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);

  /*  find scaling factor */
  ii = 0;
  while (ii <= thetas->size[1] - 1) {
    i0 = thetas->size[1];
    i = (int32_T)(1.0 + (real_T)ii);
    emlrtDynamicBoundsCheckFastR2012b(i, 1, i0, &emlrtBCI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
    factor_remain2 = 1.0;
    dlm1 = 0.0;
    emlrtForLoopVectorCheckR2012b(1.0, 1.0, m, mxDOUBLE_CLASS, (int32_T)m,
      &e_emlrtRTEI, emlrtRootTLSGlobal);
    i = 0;
    while (i <= (int32_T)m - 1) {
      factor_remain2 *= muDoubleScalarSin(thetas->data[ii]);
      if (muDoubleScalarAbs(factor_remain2) < 0.01) {
        while (muDoubleScalarAbs(factor_remain2) < 0.01) {
          dlm1 -= 2.0;
          factor_remain2 *= 100.0;
          emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar,
            emlrtRootTLSGlobal);
        }
      }

      i++;
      emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
    }

    emlrtPopRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
    temp = power_10_first + dlm1;
    emlrtPushRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
    if (muDoubleScalarFloor(m) != m) {
      emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
      b_eml_error();
      emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
    }

    emlrtPopRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPopRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
    b_y = (2.0 * m + 1.0) / 2.0;
    if (b_y < 0.0) {
      emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
      eml_error();
      emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
    }

    dlm = Km * factor_remain2 * muDoubleScalarPower(-1.0, m) *
      muDoubleScalarSqrt(b_y) / 2.5066282746310002;
    emlrtPopRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
    dlm_1 = 0.0;
    i0 = (int32_T)((L - 1.0) + (1.0 - m));
    emlrtForLoopVectorCheckR2012b(m, 1.0, L - 1.0, mxDOUBLE_CLASS, i0,
      &c_emlrtRTEI, emlrtRootTLSGlobal);
    ell = 0;
    while (ell <= i0 - 1) {
      b_ell = m + (real_T)ell;
      if (muDoubleScalarAbs(dlm) < 0.01) {
        while (muDoubleScalarAbs(dlm) < 0.01) {
          temp -= 2.0;
          dlm *= 100.0;
          dlm_1 *= 100.0;
          emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar,
            emlrtRootTLSGlobal);
        }
      }

      if (muDoubleScalarAbs(dlm) > 100.0) {
        while (muDoubleScalarAbs(dlm) > 100.0) {
          temp += 2.0;
          dlm /= 100.0;
          dlm_1 /= 100.0;
          emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar,
            emlrtRootTLSGlobal);
        }
      }

      i = P->size[0];
      factor_remain2 = (b_ell - m) + 1.0;
      i1 = (int32_T)emlrtIntegerCheckFastR2012b(factor_remain2, &c_emlrtDCI,
        emlrtRootTLSGlobal);
      i2 = P->size[1];
      i3 = 1 + ii;
      P->data[(emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &b_emlrtBCI,
                emlrtRootTLSGlobal) + P->size[0] *
               (emlrtDynamicBoundsCheckFastR2012b(i3, 1, i2, &c_emlrtBCI,
                 emlrtRootTLSGlobal) - 1)) - 1] = dlm;
      i = Sc->size[0];
      factor_remain2 = (b_ell - m) + 1.0;
      i1 = (int32_T)emlrtIntegerCheckFastR2012b(factor_remain2, &d_emlrtDCI,
        emlrtRootTLSGlobal);
      i2 = Sc->size[1];
      i3 = 1 + ii;
      Sc->data[(emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &d_emlrtBCI,
                 emlrtRootTLSGlobal) + Sc->size[0] *
                (emlrtDynamicBoundsCheckFastR2012b(i3, 1, i2, &e_emlrtBCI,
                  emlrtRootTLSGlobal) - 1)) - 1] = temp;
      emlrtPushRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
      if (b_ell == 0.0) {
        emlrtPushRtStackR2012b(&m_emlrtRSI, emlrtRootTLSGlobal);
        dlm1 = m * m;
        if (1.0 - dlm1 < 0.0) {
          emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
          eml_error();
          emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        }

        dlm1 = 1.7320508075688772 * (1.0 / muDoubleScalarSqrt(1.0 - dlm1)) *
          (muDoubleScalarCos(thetas->data[ii]) * dlm);
        emlrtPopRtStackR2012b(&m_emlrtRSI, emlrtRootTLSGlobal);
      } else {
        emlrtPushRtStackR2012b(&n_emlrtRSI, emlrtRootTLSGlobal);
        b_y = (2.0 * b_ell + 3.0) / (2.0 * b_ell + 1.0);
        if (b_y < 0.0) {
          emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
          eml_error();
          emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        }

        factor_remain2 = ((b_ell + 1.0) * (b_ell + 1.0) - m * m) * ((b_ell + 1.0)
          * (b_ell + 1.0));
        if (factor_remain2 < 0.0) {
          emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
          eml_error();
          emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        }

        dlm1 = (2.0 * b_ell + 3.0) / (2.0 * b_ell - 1.0);
        if (dlm1 < 0.0) {
          emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
          eml_error();
          emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        }

        c_y = (b_ell * b_ell - m * m) * (b_ell * b_ell);
        if (c_y < 0.0) {
          emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
          eml_error();
          emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        }

        d_y = ((b_ell + 1.0) * (b_ell + 1.0) - m * m) * ((b_ell + 1.0) * (b_ell
          + 1.0));
        if (d_y < 0.0) {
          emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
          eml_error();
          emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        }

        dlm1 = muDoubleScalarSqrt(b_y) * ((b_ell + 1.0) * (2.0 * b_ell + 1.0) /
          muDoubleScalarSqrt(factor_remain2)) * (muDoubleScalarCos(thetas->
          data[ii]) * dlm) - muDoubleScalarSqrt(dlm1) * (muDoubleScalarSqrt(c_y)
          / muDoubleScalarSqrt(d_y)) * (b_ell + 1.0) * dlm_1 / b_ell;
        emlrtPopRtStackR2012b(&n_emlrtRSI, emlrtRootTLSGlobal);
      }

      emlrtPopRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
      dlm_1 = dlm;
      dlm = dlm1;
      ell++;
      emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
    }

    ii++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
  }
}

/* End of code generation (nsht_legmat.c) */
