#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>

#include "nr_routines.h"
#include "HJM_type.h"

void nrerror(std::string error_text )
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text.c_str());
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

int *ivector(long nl, long nh)
{
  /* allocate an int vector with subscript range v[nl..nh] */
  int *v;

  v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+1;
}

void free_ivector(int *v, long nl, long nh)
{
  /* free an int vector allocated with ivector() */
  free((char *) (v+nl-1));
}

FTYPE *dvector( long nl, long nh )
{
  // allocate a FTYPE vector with subscript range v[nl..nh]
  FTYPE *v;

  v=(FTYPE *)malloc((size_t)((nh - nl + 2) * sizeof(FTYPE)));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl+1;
}

void free_dvector( FTYPE *v, long nl, long nh )
{
  // free a FTYPE vector allocated with dvector()
  free((char*) (v+nl-1));
}

FTYPE **dmatrix( long nrl, long nrh, long ncl, long nch )
{
  // allocate a FTYPE matrix with subscript range m[nrl..nrh][ncl..nch]
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  FTYPE **m;

  // allocate pointers to rows
  m=(FTYPE **) malloc((size_t)((nrow+1)*sizeof(FTYPE*)));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m += 1;
  m -= nrl;

  // allocate rows and set pointers to them
  m[nrl]=(FTYPE *)malloc((size_t)((nrow*ncol + 1)*sizeof(FTYPE)));
  if(!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i = nrl + 1; i <= nrh; i++)
		m[i] = m[i - 1] + ncol;

  // return pointer to array of pointers to rows
  return m;
}

void free_dmatrix( FTYPE **m, long nrl, long nrh, long ncl, long nch )
{
  // free a FTYPE matrix allocated by dmatrix()
  free((char*)(m[nrl] + ncl - 1));
  free((char*)(m + nrl - 1));
}
