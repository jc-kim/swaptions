#ifndef __HJM_SECURITIES__
#define __HJM_SECURITIES__

#include "HJM_type.h"

int HJM_Yield_to_Forward(FTYPE *pdForward, int iN, FTYPE *pdYield);
int HJM_Drifts(FTYPE *pdTotalDrift, FTYPE **ppdDrifts, int iN, int iFactors, FTYPE dYears, FTYPE **ppdFactors);
FTYPE dMax( FTYPE dA, FTYPE dB );

#endif //__HJM_SECURITIES__
