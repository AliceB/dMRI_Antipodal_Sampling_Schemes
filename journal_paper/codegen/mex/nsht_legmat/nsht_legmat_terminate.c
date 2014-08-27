/*
 * nsht_legmat_terminate.c
 *
 * Code generation for function 'nsht_legmat_terminate'
 *
 * C source code generated on: Mon Aug 11 15:44:13 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nsht_legmat.h"
#include "nsht_legmat_terminate.h"

/* Function Definitions */
void nsht_legmat_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void nsht_legmat_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (nsht_legmat_terminate.c) */
