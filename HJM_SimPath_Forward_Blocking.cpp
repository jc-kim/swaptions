#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "HJM_type.h"
#include "HJM.h"
#include "nr_routines.h"

void serialB(FTYPE **pdZ, FTYPE **randZ, int BLOCKSIZE, int iN, int iFactors)
{
  for(int l = 0; l < iFactors; l++){
    for(int b = 0; b < BLOCKSIZE; b++){
      for (int j = 1; j < iN; j++){
        pdZ[l][BLOCKSIZE*j + b]= CumNormalInv(randZ[l][BLOCKSIZE*j + b]);  /* 18% of the total executition time */
      }
    }
  }
}

int HJM_SimPath_Forward_Blocking(FTYPE **ppdHJMPath,	//Matrix that stores generated HJM path (Output)
    int iN,					//Number of time-steps
    int iFactors,			//Number of factors in the HJM framework
    FTYPE dYears,			//Number of years
    FTYPE *pdForward,		//t=0 Forward curve
    FTYPE *pdTotalDrift,	//Vector containing total drift corrections for different maturities
    FTYPE **ppdFactors,	//Factor volatilities
    long *lRndSeed,			//Random number seed
    int BLOCKSIZE)
{	
  //This function computes and stores an HJM Path for given inputs

  int iSuccess = 0;
  int i, j, l; //looping variables
  int b;
  FTYPE **pdZ; //vector to store random normals
  FTYPE **randZ; //vector to store random normals
  FTYPE dTotalShock; //total shock by which the forward curve is hit at (t, T-t)
  FTYPE ddelt, sqrt_ddelt; //length of time steps	

  ddelt = (FTYPE)(dYears/iN);
  sqrt_ddelt = sqrt(ddelt);

  pdZ   = dmatrix(iFactors, iN * BLOCKSIZE);
  randZ = dmatrix(iFactors, iN * BLOCKSIZE);

  // t=0 forward curve stored iN first row of ppdHJMPath
  // At time step 0: insert expected drift 
  // rest reset to 0
  for(b = 0; b < BLOCKSIZE; b++){
    for(j = 0; j < iN; j++){
      ppdHJMPath[0][BLOCKSIZE*j + b] = pdForward[j]; 

      //initializing HJMPath to zero
      for(i = 1; i < iN; i++) ppdHJMPath[i][BLOCKSIZE*j + b] = 0;
    }
  }

  // sequentially generating random numbers
  //compute random number in exact same sequence
  /* 10% of the total executition time */
  for(b = 0; b < BLOCKSIZE; b++){
    for(j = 1; j < iN; j++){
      for(l = 0; l < iFactors; l++){
        randZ[l][BLOCKSIZE*j + b] = RanUnif(lRndSeed);
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
        for (i = 0; i < iFactors; i++) dTotalShock += ppdFactors[i][l] * pdZ[i][BLOCKSIZE * j + b];

        ppdHJMPath[j][BLOCKSIZE * l + b] = ppdHJMPath[j - 1][BLOCKSIZE * (l + 1) + b] + pdTotalDrift[l] * ddelt + sqrt_ddelt * dTotalShock;
        //as per formula
      }
    }
  }

  free_dmatrix(pdZ);
  free_dmatrix(randZ);
  return 1;
}
