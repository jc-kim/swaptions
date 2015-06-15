//HJM.cpp
//Routine to setup HJM framework.
//Authors: Mark Broadie, Jatin Dewanwala, Columbia University 
//Collaborator: Mikhail Smelyanskiy, Intel
//Based on hjm_simn.xls created by Mark Broadie

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "nr_routines.h"
#include "HJM.h"
#include "HJM_type.h"

int HJM_Yield_to_Forward(FTYPE *pdForward, int iN, FTYPE *pdYield);
int HJM_Drifts(FTYPE *pdTotalDrift, FTYPE **ppdDrifts, int iN, int iFactors, FTYPE dYears, FTYPE **ppdFactors);

int HJM_Yield_to_Forward (FTYPE *pdForward,	//Forward curve to be outputted
    int iN,				//Number of time-steps
    FTYPE *pdYield)		//Input yield curve 
{
  //This function computes forward rates from supplied yield rates.
  int i;

  //forward curve computation
  pdForward[0] = pdYield[0];
  for(i = 1; i < iN; i++) pdForward[i] = (i + 1) * pdYield[i] - i * pdYield[i - 1];	//as per formula

  return 1;
}

int HJM_Drifts(FTYPE *pdTotalDrift,	//Output vector that stores the total drift correction for each maturity
    FTYPE **ppdDrifts,		//Output matrix that stores drift correction for each factor for each maturity
    int iN, 
    int iFactors,
    FTYPE dYears,
    FTYPE **ppdFactors)		//Input factor volatilities
{
  //This function computes drift corrections required for each factor for each maturity based on given factor volatilities

  int i, j, l;
  FTYPE ddelt = (FTYPE)(dYears/iN);
  FTYPE dSumVol;

  //computation of factor drifts for shortest maturity
  for (i = 0; i < iFactors; i++)
    ppdDrifts[i][0] = 0.5*ddelt*(ppdFactors[i][0])*(ppdFactors[i][0]);

  //computation of factor drifts for other maturities
  for (i = 0; i < iFactors;i++)
    for (j = 1; j < iN - 1; j++)
    {
      ppdDrifts[i][j] = 0;
      for(l = 0; l < j; l++) ppdDrifts[i][j] -= ppdDrifts[i][l];

      dSumVol=0;
      for(l = 0; l <= j; l++) dSumVol += ppdFactors[i][l];

      ppdDrifts[i][j] += 0.5 * ddelt * dSumVol * dSumVol;
    }

  //computation of total drifts for all maturities
  for(i = 0; i < iN - 1; i++)
  {
    pdTotalDrift[i]=0;
    for(j = 0; j < iFactors; j++) pdTotalDrift[i] += ppdDrifts[j][i];
  }

  return 1;
}

int Discount_Factors_Blocking(FTYPE *pdDiscountFactors, 
    int iN, 
    FTYPE dYears, 
    FTYPE *pdRatePath,
    int BLOCKSIZE) {
  int i, j, b;
  FTYPE ddelt = (FTYPE)(dYears/iN); //HJM time-step length

  FTYPE *pdexpRes = dvector((iN - 1) * BLOCKSIZE);

  //precompute the exponientials
  for (j = 0; j < (iN - 1) * BLOCKSIZE; j++) pdexpRes[j] = -pdRatePath[j] * ddelt;
  for (j = 0; j < (iN - 1) * BLOCKSIZE; j++) pdexpRes[j] = exp(pdexpRes[j]);

  //initializing the discount factor vector
  for (i = 0; i < iN * BLOCKSIZE; i++) pdDiscountFactors[i] = 1.0;

  for (i = 1; i < iN; i++){
    for (b = 0; b < BLOCKSIZE; b++){
      for (j = 0; j < i; j++){
        pdDiscountFactors[i * BLOCKSIZE + b] *= pdexpRes[j * BLOCKSIZE + b];
      }
    }
  } 

  free_dvector(pdexpRes);
  return 1;
}
