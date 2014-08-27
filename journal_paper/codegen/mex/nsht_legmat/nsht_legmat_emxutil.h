/*
 * nsht_legmat_emxutil.h
 *
 * Code generation for function 'nsht_legmat_emxutil'
 *
 * C source code generated on: Mon Aug 11 15:44:13 2014
 *
 */

#ifndef __NSHT_LEGMAT_EMXUTIL_H__
#define __NSHT_LEGMAT_EMXUTIL_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "nsht_legmat_types.h"

/* Function Declarations */
extern void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize, const emlrtRTEInfo *srcLocation);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, const emlrtRTEInfo *srcLocation, boolean_T doPush);
#endif
/* End of code generation (nsht_legmat_emxutil.h) */
