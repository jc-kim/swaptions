//HJM_Swaption_Blocking.cpp
//Routines to compute various security prices using HJM framework (via Simulation).
//Authors: Mark Broadie, Jatin Dewanwala
//Collaborator: Mikhail Smelyanskiy, Intel, Jike Chong (Berkeley)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr_routines.h"
#include "HJM_Securities.h"
#include "HJM.h"
#include "HJM_type.h"

int HJM_Swaption_Blocking(FTYPE *pdSwaptionPrice, //Output vector that will store simulation results in the form:
    //Swaption Price
    //Swaption Standard Error
    //Swaption Parameters
    FTYPE dStrike,
    FTYPE dCompounding,     //Compounding convention used for quoting the strike (0 => continuous,
    //0.5 => semi-annual, 1 => annual).
    FTYPE dMaturity,	      //Maturity of the swaption (time to expiration)
    FTYPE dTenor,	      //Tenor of the swap
    FTYPE dPaymentInterval, //frequency of swap payments e.g. dPaymentInterval = 0.5 implies a swap payment every half
    //year
    //HJM Framework Parameters (please refer HJM.cpp for explanation of variables and functions)
    int iN,
    int iFactors,
    FTYPE dYears,
    FTYPE *pdYield,
    FTYPE *ppdFactors,
    //Simulation Parameters
    long iRndSeed,
    long lTrials,
    int BLOCKSIZE, int tid,
		FTYPE *ppdHJMPath, FTYPE *pdForward, FTYPE *ppdDrifts, FTYPE *pdTotalDrift,
		FTYPE *pdDiscountingRatePath, FTYPE *pdPayoffDiscountFactors, FTYPE *pdSwapRatePath, FTYPE *pdSwapDiscountFactors, FTYPE *pdSwapPayoffs,
		FTYPE *pdZ, FTYPE *randZ,
		FTYPE *pdexpRes)
{
  int iSuccess = 0;
  int i;
  int b; //block looping variable
  long l; //looping variables

  FTYPE ddelt = (FTYPE)(dYears/iN);
  int iFreqRatio = (int)(dPaymentInterval/ddelt + 0.5);
  FTYPE dStrikeCont;
	dStrikeCont = dCompounding == 0 ? dStrike: (1 / dCompounding) * log(1 + dStrike * dCompounding);  

	int iSwapVectorLength = (int)(iN - dMaturity / ddelt + 0.5);

  int iSwapStartTimeIndex;
  int iSwapTimePoints;
  FTYPE dSwapVectorYears;

  FTYPE dSwaptionPayoff;
  FTYPE dDiscSwaptionPayoff;
  FTYPE dFixedLegValue;

  // Accumulators
  FTYPE dSumSimSwaptionPrice; 
  FTYPE dSumSquareSimSwaptionPrice;

  // Final returned results
  FTYPE dSimSwaptionMeanPrice;
  FTYPE dSimSwaptionStdError;

  iSwapStartTimeIndex = (int)(dMaturity/ddelt + 0.5);	//Swap starts at swaption maturity
  iSwapTimePoints = (int)(dTenor/ddelt + 0.5);			//Total HJM time points corresponding to the swap's tenor
  dSwapVectorYears = (FTYPE)(iSwapVectorLength*ddelt);

  for (i = 0; i < iSwapVectorLength; i++) pdSwapPayoffs[i] = 0.0;
  for (i = iFreqRatio; i <= iSwapTimePoints; i += iFreqRatio) {
    if(i != iSwapTimePoints)
      pdSwapPayoffs[i] = exp(dStrikeCont * dPaymentInterval) - 1;
    else
      pdSwapPayoffs[i] = exp(dStrikeCont * dPaymentInterval);
  }

  //generating forward curve at t=0 from supplied yield curve
  iSuccess = HJM_Yield_to_Forward(pdForward, iN, pdYield);
  if(iSuccess != 1) return iSuccess;

  //computation of drifts from factor volatilities
  iSuccess = HJM_Drifts(pdTotalDrift, ppdDrifts, iN, iFactors, dYears, ppdFactors);
  if(iSuccess != 1) return iSuccess;

  dSumSimSwaptionPrice = 0.0;
  dSumSquareSimSwaptionPrice = 0.0;

  //Simulations begin:
  for(l = 0; l < lTrials; l += BLOCKSIZE) {
    //For each trial a new HJM Path is generated
    /* GC: 51% of the time goes here */
    iSuccess = HJM_SimPath_Forward_Blocking(ppdHJMPath, iN, iFactors, dYears, pdForward, pdTotalDrift, ppdFactors, &iRndSeed, BLOCKSIZE, pdZ, randZ);
    if(iSuccess != 1) return iSuccess;

    //now we compute the discount factor vector
    for(i = 0; i < iN; i++){
      for(b = 0; b < BLOCKSIZE; b++){
        pdDiscountingRatePath[BLOCKSIZE*i + b] = ppdHJMPath[i * iN * BLOCKSIZE + b];
      }
    }

    /* 15% of the time goes here */
    iSuccess = Discount_Factors_Blocking(pdPayoffDiscountFactors, iN, dYears, pdDiscountingRatePath, BLOCKSIZE, pdexpRes);
    if(iSuccess != 1) return iSuccess;

    //now we compute discount factors along the swap path
    for(i = 0; i < iSwapVectorLength; i++){
      for(b = 0; b < BLOCKSIZE; b++){
        pdSwapRatePath[i * BLOCKSIZE + b] = ppdHJMPath[iSwapStartTimeIndex * iN * BLOCKSIZE + i * BLOCKSIZE + b];
      }
    }

    iSuccess = Discount_Factors_Blocking(pdSwapDiscountFactors, iSwapVectorLength, dSwapVectorYears, pdSwapRatePath, BLOCKSIZE, pdexpRes);
    if(iSuccess != 1) return iSuccess;

    // Simulation
    for(b = 0; b < BLOCKSIZE; b++){
      dFixedLegValue = 0.0;

      for(i = 0; i < iSwapVectorLength; i++) dFixedLegValue += pdSwapPayoffs[i] * pdSwapDiscountFactors[i*BLOCKSIZE + b];
      dSwaptionPayoff = dMax(dFixedLegValue - 1.0, 0);

      dDiscSwaptionPayoff = dSwaptionPayoff * pdPayoffDiscountFactors[iSwapStartTimeIndex*BLOCKSIZE + b];
      // end simulation

      // accumulate into the aggregating variables
      dSumSimSwaptionPrice += dDiscSwaptionPayoff;
      dSumSquareSimSwaptionPrice += dDiscSwaptionPayoff*dDiscSwaptionPayoff;
    }
  }

  // Simulation Results Stored
  dSimSwaptionMeanPrice = dSumSimSwaptionPrice / lTrials;
  dSimSwaptionStdError = sqrt((dSumSquareSimSwaptionPrice-dSumSimSwaptionPrice*dSumSimSwaptionPrice/lTrials)/
      (lTrials-1.0))/sqrt((FTYPE)lTrials);

  //results returned
  pdSwaptionPrice[0] = dSimSwaptionMeanPrice;
  pdSwaptionPrice[1] = dSimSwaptionStdError;

  return 1;
}
