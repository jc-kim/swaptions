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

int HJM_SimPath_Forward(FTYPE **ppdHJMPath, int iN, int iFactors, FTYPE dYears, FTYPE *pdForward, FTYPE *pdTotalDrift, FTYPE **ppdFactors, long *lRndSeed);
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

int HJM_SimPath_Forward(FTYPE **ppdHJMPath,	//Matrix that stores generated HJM path (Output)
    int iN,					//Number of time-steps
    int iFactors,			//Number of factors in the HJM framework
    FTYPE dYears,			//Number of years
    FTYPE *pdForward,		//t=0 Forward curve
    FTYPE *pdTotalDrift,	//Vector containing total drift corrections for different maturities
    FTYPE **ppdFactors,	//Factor volatilities
    long *lRndSeed)			//Random number seed
{	
  //This function computes and stores an HJM Path for given inputs
  int i,j,l;
  FTYPE ddelt; //length of time steps
  FTYPE dTotalShock; //total shock by which the forward curve is hit at (t, T-t)
  FTYPE *pdZ; //vector to store random normals

  ddelt = (FTYPE)(dYears/iN);

  pdZ = dvector(0, iFactors - 1); //assigning memory

  for(i = 0; i < iN; i++)
    for(j = 0; j < iN; j++)
      ppdHJMPath[i][j]=0; //initializing HJMPath to zero

  //t=0 forward curve stored iN first row of ppdHJMPath
  for(i = 0; i < iN; i++) ppdHJMPath[0][i] = pdForward[i]; 

  //Generation of HJM Path
  for(j = 1; j < iN; j++)
  {
		//shocks to hit various factors for forward curve at t
    for (l = 0; l < iFactors; l++) pdZ[l]= CumNormalInv(RanUnif(lRndSeed));

    for (l = 0; l < iN - j; l++)
    {
      dTotalShock = 0;
      for (i = 0; i < iFactors; i++) dTotalShock += ppdFactors[i][l]* pdZ[i];

      ppdHJMPath[j][l] = ppdHJMPath[j - 1][l + 1] + (pdTotalDrift[l] * ddelt) + (sqrt(ddelt) * dTotalShock);
      //as per formula
    }
  }

  free_dvector(pdZ, 0, iFactors - 1);
  return 1;
}

int Discount_Factors_Blocking(FTYPE *pdDiscountFactors, 
    int iN, 
    FTYPE dYears, 
    FTYPE *pdRatePath,
    int BLOCKSIZE) {
  int i,j,b;
  FTYPE ddelt;			//HJM time-step length
  ddelt = (FTYPE)(dYears/iN);

  FTYPE *pdexpRes;
  pdexpRes = dvector(0,(iN-1)*BLOCKSIZE-1);

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

  free_dvector(pdexpRes, 0, (iN - 1) * BLOCKSIZE - 1);
  return 1;
}
