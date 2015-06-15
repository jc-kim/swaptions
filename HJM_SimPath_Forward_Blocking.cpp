#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "HJM_type.h"
#include "HJM.h"
#include "nr_routines.h"

void serialB(FTYPE *pdZ, FTYPE *randZ, int BLOCKSIZE, int iN, int iFactors)
{
  for(int l = 0; l < iFactors; l++){
    for(int b = 0; b < BLOCKSIZE; b++){
      for (int j = 1; j < iN; j++){
        pdZ[l * iN * BLOCKSIZE + BLOCKSIZE*j + b]= CumNormalInv(randZ[l * iN * BLOCKSIZE + BLOCKSIZE * j + b]);  /* 18% of the total executition time */
      }
    }
  }
}

int HJM_SimPath_Forward_Blocking(FTYPE *ppdHJMPath,	//Matrix that stores generated HJM path (Output)
    int iN,					//Number of time-steps
    int iFactors,			//Number of factors in the HJM framework
    FTYPE dYears,			//Number of years
    FTYPE *pdForward,		//t=0 Forward curve
    FTYPE *pdTotalDrift,	//Vector containing total drift corrections for different maturities
    FTYPE *ppdFactors,	//Factor volatilities
    long *lRndSeed,			//Random number seed
    int BLOCKSIZE,
		FTYPE *pdZ, FTYPE *randZ)
{	
  //This function computes and stores an HJM Path for given inputs

  int iSuccess = 0;
  int i, j, l; //looping variables
  int b;
  FTYPE dTotalShock; //total shock by which the forward curve is hit at (t, T-t)
  FTYPE ddelt, sqrt_ddelt; //length of time steps	

  ddelt = (FTYPE)(dYears/iN);
  sqrt_ddelt = sqrt(ddelt);

  // t=0 forward curve stored iN first row of ppdHJMPath
  // At time step 0: insert expected drift 
  // rest reset to 0
  for(b = 0; b < BLOCKSIZE; b++){
    for(j = 0; j < iN; j++){
      ppdHJMPath[BLOCKSIZE*j + b] = pdForward[j]; 

      //initializing HJMPath to zero
      for(i = 1; i < iN; i++) ppdHJMPath[i * iN * BLOCKSIZE + BLOCKSIZE*j + b] = 0;
    }
  }

  // sequentially generating random numbers
  //compute random number in exact same sequence
  /* 10% of the total executition time */
  for(b = 0; b < BLOCKSIZE; b++){
    for(j = 1; j < iN; j++){
      for(l = 0; l < iFactors; l++){
        randZ[l * iN * BLOCKSIZE + BLOCKSIZE * j + b] = RanUnif(lRndSeed);
      }
    }
  }

  // shocks to hit various factors for forward curve at t
  /* 18% of the total executition time */
  serialB(pdZ, randZ, BLOCKSIZE, iN, iFactors);

  // Generation of HJM Path1
  for(b = 0; b < BLOCKSIZE; b++){ // b is the blocks
    for(j = 1; j < iN; j++) { // j is the timestep
      for(l = 0; l < iN - j; l++){ // l is the future steps
        dTotalShock = 0;

        // i steps through the stochastic factors
        for (i = 0; i < iFactors; i++) dTotalShock += ppdFactors[i * (iN - 1) + l] * pdZ[i * iN * BLOCKSIZE + BLOCKSIZE * j + b];

        ppdHJMPath[j * iN * BLOCKSIZE + BLOCKSIZE * l + b] = ppdHJMPath[(j - 1) * iN * BLOCKSIZE + BLOCKSIZE * (l + 1) + b] + pdTotalDrift[l] * ddelt + sqrt_ddelt * dTotalShock;
        //as per formula
      }
    }
  }

  return 1;
}
