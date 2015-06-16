//HJM_Securities.cpp
//Routines to compute various security prices using HJM framework (via Simulation).
//Authors: Mark Broadie, Jatin Dewanwala
//Collaborator: Mikhail Smelyanskiy, Jike Chong, Intel

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <iostream>

#include "nr_routines.h"
#include "HJM.h"
#include "HJM_Securities.h"
#include "HJM_type.h"

#ifdef ENABLE_THREADS
#include <pthread.h>
#define MAX_THREAD 1024
#endif //ENABLE_THREADS
#ifdef ENABLE_OPENCL
#include <CL/cl.h>
#endif

int NUM_TRIALS = DEFAULT_NUM_TRIALS;
int nThreads = 1;
int nSwaptions = 1;
int iN = 11;
FTYPE dYears = 5.5; 
int iFactors = 3; 
parm *swaptions;

#ifdef ENABLE_OPENCL
char* getKernelSrc() {
	FILE* fp = fopen("kernel.cl", "r");
	char *result;
	long file_size;
	fseek(fp, 0, SEEK_END);
	file_size = ftell(fp);
	rewind(fp);

	result = (char*)malloc(sizeof(char) * (file_size + 1));
	fread(result, sizeof(char), file_size, fp);
  result[file_size] = '\0';
	fclose(fp);

	return result;
}

void checkError(cl_int err) {
  if(err != CL_SUCCESS) {
    printf("Error occured! %d\n", err);
    exit(1);
  }
}

void checkError(cl_int err, int line_no) {
  if(err != CL_SUCCESS) {
    printf("Error occured in line %d!: %d\n", line_no, err);
  }
}
#endif

void* worker(void *arg){
  int tid = *((int *)arg);

  int chunksize = nSwaptions / nThreads;
  int beg = tid * chunksize;
  int end = (tid == nThreads - 1 ? nSwaptions : (tid + 1) * chunksize);

#ifdef ENABLE_OPENCL
#ifdef DEVICE_CPU
  cl_device_type device_type = CL_DEVICE_TYPE_CPU;
#else
  cl_device_type device_type = CL_DEVICE_TYPE_GPU;
#endif
  FTYPE* swaptionPrices = dmatrix(chunksize, 2);
  cl_platform_id platform;
  cl_device_id* devices;
  cl_uint device_count;
  cl_program program;
  cl_kernel kernel;
  cl_context context;
  cl_command_queue* queues;
  cl_int err;

  checkError(clGetPlatformIDs(1, &platform, NULL));
  checkError(clGetDeviceIDs(platform, device_type, NULL, NULL, &device_count));
  devices = (cl_device_id*)malloc(sizeof(cl_device_id) * device_count);
  checkError(clGetDeviceIDs(platform, device_type, device_count, devices, NULL));
  context = clCreateContext(NULL, device_count, devices, NULL, NULL, &err);
  checkError(err);

  char* kernelSource = getKernelSrc();
  size_t kernelLength = strlen(kernelSource);

  program = clCreateProgramWithSource(context, 1, (const char**)&kernelSource, &kernelLength, &err);
  checkError(err);
  err = clBuildProgram(program, device_count, devices, NULL, NULL, NULL);
  if(err != CL_SUCCESS) {
    printf("Kernel Program build failed: \n");
    size_t len;
    char *buffer;
    cl_build_status status;
    buffer = (char*)calloc(2048,sizeof(char));
    clGetProgramBuildInfo(program, *devices, CL_PROGRAM_BUILD_LOG, 2048*sizeof(char), buffer, &len);
    printf("\tLog: %s\n", buffer);
    clGetProgramBuildInfo(program, *devices, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, &len);
    printf("\tSTATUS: %s\n", status == CL_BUILD_SUCCESS ? "success" : (status == CL_BUILD_ERROR ? "error": (status == CL_BUILD_NONE ? "none": "progress")));
    clGetProgramBuildInfo(program, *devices, CL_PROGRAM_BUILD_OPTIONS, 2048*sizeof(char), buffer, &len);
    printf("\tOPTIONS: %s\n", buffer);
    exit(1);
  }
  kernel = clCreateKernel(program, "HJM_Swaption_Blocking", &err);
  checkError(err);

  queues = (cl_command_queue*)malloc(sizeof(cl_command_queue) * device_count);
  for(int i = 0; i < device_count; i++) {
    queues[i] = clCreateCommandQueue(context, devices[i], CL_QUEUE_PROFILING_ENABLE, &err);
    checkError(err);
  }

  int BLOCKSIZE = BLOCK_SIZE;
  int i, device_idx;
  for(i = beg, device_idx = 0; i < end; device_idx = (device_idx + 1) % device_count, i++) {
    FTYPE ddelt = (FTYPE)(swaptions[i].dYears / swaptions[i].iN);
    int iSwapVectorLength = (int)(swaptions[i].iN - swaptions[i].dMaturity / ddelt + 0.5);

    cl_command_queue queue = queues[device_idx];

    cl_mem pdSwaptionPriceBuf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(FTYPE) * 2, NULL, &err);
    checkError(err);
    cl_mem pdYieldBuf = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(FTYPE) * swaptions[i].iN, NULL, &err);
    checkError(err);
    cl_mem ppdFactorsBuf = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(FTYPE) * swaptions[i].iFactors * (swaptions[i].iN - 1), NULL, &err);
    checkError(err);
    long lRndSeed = 100;
    long trials = NUM_TRIALS;

    checkError(clSetKernelArg(kernel, 0, sizeof(cl_mem), &pdSwaptionPriceBuf), __LINE__);
    checkError(clSetKernelArg(kernel, 1, sizeof(FTYPE), &swaptions[i].dStrike), __LINE__);
    checkError(clSetKernelArg(kernel, 2, sizeof(FTYPE), &swaptions[i].dCompounding), __LINE__);
    checkError(clSetKernelArg(kernel, 3, sizeof(FTYPE), &swaptions[i].dMaturity), __LINE__);
    checkError(clSetKernelArg(kernel, 4, sizeof(FTYPE), &swaptions[i].dTenor), __LINE__);
    checkError(clSetKernelArg(kernel, 5, sizeof(FTYPE), &swaptions[i].dPaymentInterval), __LINE__);
    checkError(clSetKernelArg(kernel, 6, sizeof(int), &swaptions[i].iN), __LINE__);
    checkError(clSetKernelArg(kernel, 7, sizeof(int), &swaptions[i].iFactors), __LINE__);
    checkError(clSetKernelArg(kernel, 8, sizeof(FTYPE), &swaptions[i].dYears), __LINE__);
    checkError(clSetKernelArg(kernel, 9, sizeof(cl_mem), &pdYieldBuf), __LINE__);
    checkError(clSetKernelArg(kernel, 10, sizeof(cl_mem), &ppdFactorsBuf), __LINE__);
    checkError(clSetKernelArg(kernel, 11, sizeof(long), &lRndSeed), __LINE__);
    checkError(clSetKernelArg(kernel, 12, sizeof(long), &trials), __LINE__);
    checkError(clSetKernelArg(kernel, 13, sizeof(int), &BLOCKSIZE), __LINE__);
    checkError(clSetKernelArg(kernel, 14, sizeof(int), &tid), __LINE__);
    checkError(clSetKernelArg(kernel, 15, sizeof(FTYPE) * swaptions[i].iN * swaptions[i].iN * BLOCK_SIZE, NULL), __LINE__); // ppdHJMPath
    checkError(clSetKernelArg(kernel, 16, sizeof(FTYPE) * swaptions[i].iN, NULL), __LINE__); // pdForward
    checkError(clSetKernelArg(kernel, 17, sizeof(FTYPE) * swaptions[i].iFactors * (swaptions[i].iN - 1), NULL), __LINE__); // ppdDrifts
    checkError(clSetKernelArg(kernel, 18, sizeof(FTYPE) * (swaptions[i].iN - 1), NULL), __LINE__); // pdTotalDrift
    checkError(clSetKernelArg(kernel, 19, sizeof(FTYPE) * swaptions[i].iN * BLOCK_SIZE, NULL), __LINE__); // pdDiscountingRatePath
    checkError(clSetKernelArg(kernel, 20, sizeof(FTYPE) * swaptions[i].iN * BLOCK_SIZE, NULL), __LINE__); // pdPayoffDiscountFactors
    checkError(clSetKernelArg(kernel, 21, sizeof(FTYPE) * iSwapVectorLength * BLOCK_SIZE, NULL), __LINE__); // pdSwapRatePath
    checkError(clSetKernelArg(kernel, 22, sizeof(FTYPE) * iSwapVectorLength * BLOCK_SIZE, NULL), __LINE__); // pdSwapDiscountFactors
    checkError(clSetKernelArg(kernel, 23, sizeof(FTYPE) * iSwapVectorLength, NULL), __LINE__); // pdSwapPayoffs
    checkError(clSetKernelArg(kernel, 24, sizeof(FTYPE) * swaptions[i].iFactors * swaptions[i].iN * BLOCK_SIZE, NULL), __LINE__); // pdZ
    checkError(clSetKernelArg(kernel, 25, sizeof(FTYPE) * swaptions[i].iFactors * swaptions[i].iN * BLOCK_SIZE, NULL), __LINE__); // randZ
    checkError(clSetKernelArg(kernel, 26, sizeof(FTYPE) * (swaptions[i].iN - 1) * BLOCK_SIZE, NULL), __LINE__); // pdexpRes

    checkError(clEnqueueWriteBuffer(queue, pdYieldBuf, CL_FALSE, 0, sizeof(FTYPE) * swaptions[i].iN, swaptions[i].pdYield, NULL, NULL, NULL));
    checkError(clEnqueueWriteBuffer(queue, ppdFactorsBuf, CL_FALSE, 0, sizeof(FTYPE) * swaptions[i].iFactors * (swaptions[i].iN - 1), swaptions[i].ppdFactors, NULL, NULL, NULL));

    size_t group_size[1] = { 1 };
    size_t work_items[1] = { 1 };
    checkError(clEnqueueNDRangeKernel(queues[device_idx], kernel, 1, 0, group_size, work_items, 0, NULL, NULL));
    checkError(clEnqueueReadBuffer(queue, pdSwaptionPriceBuf, CL_FALSE, 0, sizeof(FTYPE) * 2, (void*)&swaptionPrices[(i - beg) * 2], NULL, NULL, NULL));

    clReleaseMemObject(ppdFactorsBuf);
    clReleaseMemObject(pdYieldBuf);
    clReleaseMemObject(pdSwaptionPriceBuf);
  }
  for(i = 0; i < device_count; i++) {
    clFinish(queues[i]);
  }

  for(i = beg; i < end; i++) {
    swaptions[i].dSimSwaptionMeanPrice = swaptionPrices[(i - beg) * 2];
    swaptions[i].dSimSwaptionStdError = swaptionPrices[(i - beg) * 2 + 1];
  }
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseContext(context);
  for(i = 0; i < device_count; i++) {
    clReleaseCommandQueue(queues[i]);
    clReleaseDevice(devices[i]);
  }
  free(swaptionPrices);
#else
  FTYPE pdSwaptionPrice[2];
  for(int i = beg; i < end; i++) {
    FTYPE ddelt = (FTYPE)(swaptions[i].dYears / swaptions[i].iN);
    int iSwapVectorLength = (int)(swaptions[i].iN - swaptions[i].dMaturity / ddelt + 0.5);
    FTYPE *ppdHJMPath = dmatrix(iN, iN * BLOCK_SIZE),
          *pdForward = dvector(iN),
          *ppdDrifts = dmatrix(iFactors, iN - 1),
          *pdTotalDrift = dvector(iN - 1);
    FTYPE *pdDiscountingRatePath = dvector(iN * BLOCK_SIZE),
          *pdPayoffDiscountFactors = dvector(iN * BLOCK_SIZE),
          *pdSwapRatePath = dvector(iSwapVectorLength * BLOCK_SIZE),
          *pdSwapDiscountFactors = dvector(iSwapVectorLength * BLOCK_SIZE),
          *pdSwapPayoffs = dvector(iSwapVectorLength);
    FTYPE *pdZ = dmatrix(iFactors, iN * BLOCK_SIZE),
          *randZ = dmatrix(iFactors, iN * BLOCK_SIZE);
    FTYPE *pdexpRes = dvector((iN - 1) * BLOCK_SIZE);

    int iSuccess = HJM_Swaption_Blocking(pdSwaptionPrice, swaptions[i].dStrike,
        swaptions[i].dCompounding, swaptions[i].dMaturity,
        swaptions[i].dTenor, swaptions[i].dPaymentInterval,
        swaptions[i].iN, swaptions[i].iFactors, swaptions[i].dYears,
        swaptions[i].pdYield, swaptions[i].ppdFactors,
        100, NUM_TRIALS, BLOCK_SIZE, 0,
        ppdHJMPath, pdForward, ppdDrifts, pdTotalDrift,
        pdDiscountingRatePath, pdPayoffDiscountFactors, pdSwapRatePath, pdSwapDiscountFactors, pdSwapPayoffs,
        pdZ, randZ, pdexpRes);
		assert(iSuccess == 1);

		free(ppdHJMPath);
		free(pdForward);
		free(ppdDrifts);
		free(pdTotalDrift);
		free(pdDiscountingRatePath);
		free(pdPayoffDiscountFactors);
		free(pdSwapRatePath);
		free(pdSwapDiscountFactors);
		free(pdSwapPayoffs);
		free(pdZ);
		free(randZ);
		free(pdexpRes);

    swaptions[i].dSimSwaptionMeanPrice = pdSwaptionPrice[0];
    swaptions[i].dSimSwaptionStdError = pdSwaptionPrice[1];
  }
#endif
  return NULL;
}

void parseOpt(int argc, char** argv) {
  if(argc == 1)
  {
    fprintf(stderr," usage: \n\t-ns [number of swaptions (should be > number of threads]\n\t-sm [number of simulations]\n\t-nt [number of threads]\n"); 
    exit(1);
  }

  for (int j=1; j<argc; j++) {
    if (!strcmp("-sm", argv[j])) {NUM_TRIALS = atol(argv[++j]);}
    else if (!strcmp("-nt", argv[j])) {nThreads = atoi(argv[++j]);} 
    else if (!strcmp("-ns", argv[j])) {nSwaptions = atoi(argv[++j]);} 
    else {
      fprintf(stderr," usage: \n\t-ns [number of swaptions (should be > number of threads]\n\t-sm [number of simulations]\n\t-nt [number of threads]\n"); 
    }
  }

  if(nSwaptions < nThreads) {
    nSwaptions = nThreads; 
  }
}

FTYPE* getFactors() {
  // initialize input dataset
  FTYPE* factors = dmatrix(iFactors, iN - 1);
  //the three rows store vol data for the three factors
  factors[0]= .01;
  factors[1]= .01;
  factors[2]= .01;
  factors[3]= .01;
  factors[4]= .01;
  factors[5]= .01;
  factors[6]= .01;
  factors[7]= .01;
  factors[8]= .01;
  factors[9]= .01;

  factors[(iN - 1)]= .009048;
  factors[(iN - 1) + 1]= .008187;
  factors[(iN - 1) + 2]= .007408;
  factors[(iN - 1) + 3]= .006703;
  factors[(iN - 1) + 4]= .006065;
  factors[(iN - 1) + 5]= .005488;
  factors[(iN - 1) + 6]= .004966;
  factors[(iN - 1) + 7]= .004493;
  factors[(iN - 1) + 8]= .004066;
  factors[(iN - 1) + 9]= .003679;

  factors[(iN - 1) * 2]= .001000;
  factors[(iN - 1) * 2 + 1]= .000750;
  factors[(iN - 1) * 2 + 2]= .000500;
  factors[(iN - 1) * 2 + 3]= .000250;
  factors[(iN - 1) * 2 + 4]= .000000;
  factors[(iN - 1) * 2 + 5]= -.000250;
  factors[(iN - 1) * 2 + 6]= -.000500;
  factors[(iN - 1) * 2 + 7]= -.000750;
  factors[(iN - 1) * 2 + 8]= -.001000;
  factors[(iN - 1) * 2 + 9]= -.001250;

	return factors;
}

void initSwaption(FTYPE* factors) {
  // setting up multiple swaptions
	int i, j, k;
  swaptions = (parm *)malloc(sizeof(parm)*nSwaptions);

  for (i = 0; i < nSwaptions; i++) {
    swaptions[i].Id = i;
    swaptions[i].iN = iN;
    swaptions[i].iFactors = iFactors;
    swaptions[i].dYears = dYears;

    swaptions[i].dStrike = (double)i / (double)nSwaptions;
    swaptions[i].dCompounding = 0;
    swaptions[i].dMaturity = 1;
    swaptions[i].dTenor = 2.0;
    swaptions[i].dPaymentInterval = 1.0;

    swaptions[i].pdYield = dvector(iN);
    swaptions[i].pdYield[0] = 0.1;
    for(j = 1; j <= swaptions[i].iN; j++) swaptions[i].pdYield[j] = swaptions[i].pdYield[j - 1] + 0.005;

    swaptions[i].ppdFactors = dmatrix(swaptions[i].iFactors, swaptions[i].iN - 1);
    for(k = 0; k < swaptions[i].iFactors; k++)
      for(j = 0; j < swaptions[i].iN - 1; j++)
        swaptions[i].ppdFactors[k * (iN - 1) + j] = factors[k * (iN - 1) + j];
  }

}
//Please note: Whenever we type-cast to (int), we add 0.5 to ensure that the value is rounded to the correct number. 
//For instance, if X/Y = 0.999 then (int) (X/Y) will equal 0 and not 1 (as (int) rounds down).
//Adding 0.5 ensures that this does not happen. Therefore we use (int) (X/Y + 0.5); instead of (int) (X/Y);
int main(int argc, char *argv[])
{
  int i, j;
  FTYPE *factors = NULL;

#ifdef PARSEC_VERSION
  printf("PARSEC Benchmark Suite\n");
  fflush(NULL);
#endif //PARSEC_VERSION

	parseOpt(argc, argv);
  printf("Number of Simulations: %d,  Number of threads: %d Number of swaptions: %d\n", NUM_TRIALS, nThreads, nSwaptions);
#ifdef ENABLE_THREADS
  pthread_t* threads;
  pthread_attr_t pthread_custom_attr;

  if ((nThreads < 1) || (nThreads > MAX_THREAD)) {
    fprintf(stderr,"Number of threads must be between 1 and %d.\n", MAX_THREAD);
    exit(1);
  }
  threads = (pthread_t*)malloc(nThreads * sizeof(pthread_t));
  pthread_attr_init(&pthread_custom_attr);

  if ((nThreads < 1) || (nThreads > MAX_THREAD)) {
    fprintf(stderr,"Number of threads must be between 1 and %d.\n", MAX_THREAD);
    exit(1);
  }
#else
  if (nThreads != 1) {
    fprintf(stderr,"Number of threads must be 1 (serial version)\n");
    exit(1);
  }
#endif

	factors = getFactors();
	initSwaption(factors);

  // Calling the Swaption Pricing Routine
#ifdef ENABLE_THREADS
  int threadIDs[nThreads];
  for (i = 0; i < nThreads; i++) {
    threadIDs[i] = i;
    pthread_create(&threads[i], &pthread_custom_attr, worker, &threadIDs[i]);
  }
  for (i = 0; i < nThreads; i++) {
    pthread_join(threads[i], NULL);
  }

  free(threads);
#else
  int threadID = 0;
  worker(&threadID);
#endif //ENABLE_THREADS

  for (i = 0; i < nSwaptions; i++) {
    fprintf(stderr,"Swaption%d: [SwaptionPrice: %.10lf StdError: %.10lf] \n", 
        i, swaptions[i].dSimSwaptionMeanPrice, swaptions[i].dSimSwaptionStdError);
  }

  for (i = 0; i < nSwaptions; i++) {
    free_dvector(swaptions[i].pdYield);
    free_dmatrix(swaptions[i].ppdFactors);
  }

  free(swaptions);

  return 0;
}
