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

FTYPE *dvector(long length) {
  FTYPE *v = (FTYPE*)malloc(sizeof(FTYPE) * length);
  if(!v) nrerror("allocation failure in dvector()");
  return v;
}

void free_dvector(FTYPE *v)
{
  // free a FTYPE vector allocated with dvector()
  free((char*)v);
}

FTYPE** dmatrix(long rows, long cols) {
  int i;
  FTYPE** m = (FTYPE **)malloc(sizeof(FTYPE *) * rows);
  m[0] = (FTYPE *)malloc(sizeof(FTYPE) * rows * cols);

  for(i = 1; i < rows; i++) {
    m[i] = m[i - 1] + cols;
  }

  return m;
}

void free_dmatrix(FTYPE **m)
{
  free((char*)m[0]);
  free((char*)m);
}
