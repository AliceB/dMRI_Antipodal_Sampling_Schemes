/*
 * nsht_legmat_initialize.c
 *
 * Code generation for function 'nsht_legmat_initialize'
 *
 * C source code generated on: Mon Aug 11 15:44:13 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nsht_legmat.h"
#include "nsht_legmat_initialize.h"
#include "nsht_legmat_data.h"

/* Function Definitions */
void nsht_legmat_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (nsht_legmat_initialize.c) */
