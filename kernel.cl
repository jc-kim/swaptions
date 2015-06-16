#define FTYPE double

FTYPE dMax(FTYPE dA, FTYPE dB) {
  return dA > dB ? dA : dB;
}

FTYPE RanUnif(long *s)
{
  // uniform random number generator
  long   ix, k1;
  FTYPE dRes;

  ix = *s;
  *s = ix+1;
  ix *= 1513517L;
  ix %= 2147483647L;
  k1 = ix/127773L;
  ix = 16807L*( ix - k1*127773L ) - k1 * 2836L;
  if (ix < 0) ix = ix + 2147483647L;
  //*s   = ix;
  dRes = (ix * 4.656612875e-10);
  return (dRes);
}

__constant FTYPE a[4] = {
  2.50662823884,
  -18.61500062529,
  41.39119773534,
  -25.44106049637
};

__constant FTYPE b[4] = {
  -8.47351093090,
  23.08336743743,
  -21.06224101826,
  3.13082909833
};

__constant FTYPE c[9] = {
  0.3374754822726147,
  0.9761690190917186,
  0.1607979714918209,
  0.0276438810333863,
  0.0038405729373609,
  0.0003951896511919,
  0.0000321767881768,
  0.0000002888167364,
  0.0000003960315187
};

FTYPE CumNormalInv(FTYPE u)
{
  FTYPE x, r;

  x = u - 0.5;
  if(fabs (x) < 0.42)
  { 
    r = x * x;
    r = x * ((( a[3]*r + a[2]) * r + a[1]) * r + a[0])/
      ((((b[3] * r+ b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
    return (r);
  }

  r = u;
  if(x > 0.0) r = 1.0 - u;
  r = log(-log(r));
  r = c[0] + r * (c[1] + r * 
      (c[2] + r * (c[3] + r * 
                   (c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r*c[8])))))));

  if(x < 0.0) r = -r;

  return (r);
}

void HJM_Yield_to_Forward (__local FTYPE *pdForward, int iN, __global FTYPE *pdYield) {
  int i;
  pdForward[0] = pdYield[0];
  for(i = 1; i < iN; i++) pdForward[i] = (i + 1) * pdYield[i] - i * pdYield[i - 1];
}

void HJM_Drifts(__local FTYPE *pdTotalDrift,
    __local FTYPE *ppdDrifts,
    int iN, 
    int iFactors,
    FTYPE dYears,
    __global FTYPE *ppdFactors) {
  int i, j, l;
  FTYPE ddelt = (FTYPE)(dYears/iN);
  FTYPE dSumVol;

  for (i = 0; i < iFactors; i++)
    ppdDrifts[i * (iN - 1)] = 0.5 * ddelt * (ppdFactors[i * (iN - 1)]) * (ppdFactors[i * (iN - 1)]);

  for (i = 0; i < iFactors;i++) {
    for (j = 1; j < iN - 1; j++) {
      ppdDrifts[i * (iN - 1) + j] = 0;
      for(l = 0; l < j; l++) ppdDrifts[i * (iN - 1) + j] -= ppdDrifts[i * (iN - 1) + l];

      dSumVol=0;
      for(l = 0; l <= j; l++) dSumVol += ppdFactors[i * (iN - 1) + l];

      ppdDrifts[i * (iN - 1) + j] += 0.5 * ddelt * dSumVol * dSumVol;
    }
  }

  for(i = 0; i < iN - 1; i++) {
    pdTotalDrift[i]=0;
    for(j = 0; j < iFactors; j++) pdTotalDrift[i] += ppdDrifts[j * (iN - 1) + i];
  }
}

void serialB(__local FTYPE *pdZ, __local FTYPE *randZ, int BLOCKSIZE, int iN, int iFactors)
{
  int l, b, j;
  for(l = 0; l < iFactors; l++){
    for(b = 0; b < BLOCKSIZE; b++){
      for (j = 1; j < iN; j++){
        pdZ[l * iN * BLOCKSIZE + BLOCKSIZE*j + b]= CumNormalInv(randZ[l * iN * BLOCKSIZE + BLOCKSIZE * j + b]);  /* 18% of the total executition time */
      }
    }
  }
}

void HJM_SimPath_Forward_Blocking(__local FTYPE *ppdHJMPath,
    int iN,
    int iFactors,
    FTYPE dYears,
    __local FTYPE *pdForward,
    __local FTYPE *pdTotalDrift,
    __global FTYPE *ppdFactors,
    long *lRndSeed,
    int BLOCKSIZE,
		__local FTYPE *pdZ, __local FTYPE *randZ) {	
  int i, j, l;
  int b;
  FTYPE dTotalShock;
  FTYPE ddelt, sqrt_ddelt;

  ddelt = (FTYPE)(dYears/iN);
  sqrt_ddelt = sqrt(ddelt);

  for(b = 0; b < BLOCKSIZE; b++){
    for(j = 0; j < iN; j++){
      ppdHJMPath[BLOCKSIZE*j + b] = pdForward[j]; 

      for(i = 1; i < iN; i++) ppdHJMPath[i * iN * BLOCKSIZE + BLOCKSIZE*j + b] = 0;
    }
  }
  for(b = 0; b < BLOCKSIZE; b++){
    for(j = 1; j < iN; j++){
      for(l = 0; l < iFactors; l++){
        randZ[l * iN * BLOCKSIZE + BLOCKSIZE * j + b] = RanUnif(lRndSeed);
      }
    }
  }
  serialB(pdZ, randZ, BLOCKSIZE, iN, iFactors);

  for(b = 0; b < BLOCKSIZE; b++){
    for(j = 1; j < iN; j++) {
      for(l = 0; l < iN - j; l++){
        dTotalShock = 0;

        for (i = 0; i < iFactors; i++) dTotalShock += ppdFactors[i * (iN - 1) + l] * pdZ[i * iN * BLOCKSIZE + BLOCKSIZE * j + b];

        ppdHJMPath[j * iN * BLOCKSIZE + BLOCKSIZE * l + b] = ppdHJMPath[(j - 1) * iN * BLOCKSIZE + BLOCKSIZE * (l + 1) + b] + pdTotalDrift[l] * ddelt + sqrt_ddelt * dTotalShock;
      }
    }
  }
}

void Discount_Factors_Blocking(__local FTYPE *pdDiscountFactors, 
    int iN, 
    FTYPE dYears, 
    __local FTYPE *pdRatePath,
    int BLOCKSIZE,
		__local FTYPE *pdexpRes) {
  int i, j, b;
  FTYPE ddelt = (FTYPE)(dYears/iN); //HJM time-step length

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
}

__kernel void HJM_Swaption_Blocking(__global FTYPE *pdSwaptionPrice,
    const FTYPE dStrike,
    const FTYPE dCompounding,
    const FTYPE dMaturity,
    const FTYPE dTenor,
    const FTYPE dPaymentInterval,
    const int iN,
    const int iFactors,
    const FTYPE dYears,
    __global FTYPE *pdYield,
    __global FTYPE *ppdFactors,
    const long lRndSeed,
    const long lTrials,
    int BLOCKSIZE, int tid,
    __local FTYPE *ppdHJMPath, __local FTYPE *pdForward, __local FTYPE *ppdDrifts, __local FTYPE *pdTotalDrift,
    __local FTYPE *pdDiscountingRatePath, __local FTYPE *pdPayoffDiscountFactors, __local FTYPE *pdSwapRatePath, __local FTYPE *pdSwapDiscountFactors, __local FTYPE *pdSwapPayoffs,
    __local FTYPE *pdZ, __local FTYPE *randZ,
    __local FTYPE *pdexpRes)
{
  int i, b;
  long l;

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
  HJM_Yield_to_Forward(pdForward, iN, pdYield);

  //computation of drifts from factor volatilities
  HJM_Drifts(pdTotalDrift, ppdDrifts, iN, iFactors, dYears, ppdFactors);

  dSumSimSwaptionPrice = 0.0;
  dSumSquareSimSwaptionPrice = 0.0;

  //Simulations begin:
  for(l = 0L; l < lTrials; l += (long)BLOCKSIZE) {
    HJM_SimPath_Forward_Blocking(ppdHJMPath, iN, iFactors, dYears, pdForward, pdTotalDrift, ppdFactors, &lRndSeed, BLOCKSIZE, pdZ, randZ);

    for(i = 0; i < iN; i++){
      for(b = 0; b < BLOCKSIZE; b++){
        pdDiscountingRatePath[BLOCKSIZE*i + b] = ppdHJMPath[i * iN * BLOCKSIZE + b];
      }
    }

    Discount_Factors_Blocking(pdPayoffDiscountFactors, iN, dYears, pdDiscountingRatePath, BLOCKSIZE, pdexpRes);

    for(i = 0; i < iSwapVectorLength; i++){
      for(b = 0; b < BLOCKSIZE; b++){
        pdSwapRatePath[i * BLOCKSIZE + b] = ppdHJMPath[iSwapStartTimeIndex * iN * BLOCKSIZE + i * BLOCKSIZE + b];
      }
    }

    Discount_Factors_Blocking(pdSwapDiscountFactors, iSwapVectorLength, dSwapVectorYears, pdSwapRatePath, BLOCKSIZE, pdexpRes);

    for(b = 0; b < BLOCKSIZE; b++){
      dFixedLegValue = 0.0;

      for(i = 0; i < iSwapVectorLength; i++) dFixedLegValue += pdSwapPayoffs[i] * pdSwapDiscountFactors[i*BLOCKSIZE + b];
      dSwaptionPayoff = dMax(dFixedLegValue - 1.0, 0);

      dDiscSwaptionPayoff = dSwaptionPayoff * pdPayoffDiscountFactors[iSwapStartTimeIndex*BLOCKSIZE + b];

      dSumSimSwaptionPrice += dDiscSwaptionPayoff;
      dSumSquareSimSwaptionPrice += dDiscSwaptionPayoff*dDiscSwaptionPayoff;
    }
  }

  dSimSwaptionMeanPrice = dSumSimSwaptionPrice / lTrials;
  dSimSwaptionStdError = sqrt((dSumSquareSimSwaptionPrice-dSumSimSwaptionPrice*dSumSimSwaptionPrice/lTrials)/
      (lTrials-1.0))/sqrt((FTYPE)lTrials);

  pdSwaptionPrice[0] = dSimSwaptionMeanPrice;
  pdSwaptionPrice[1] = dSimSwaptionStdError;
}
