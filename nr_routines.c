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
  FTYPE *v = (FTYPE*)calloc(1, sizeof(FTYPE) * length);
  if(!v) nrerror("allocation failure in dvector()");
  return v;
}

void free_dvector(FTYPE *v)
{
  free(v);
  v = NULL;
}

FTYPE* dmatrix(long rows, long cols) {
  FTYPE* m = (FTYPE *)calloc(1, sizeof(FTYPE) * rows * cols);
  if(!m) nrerror("allocation failure in dmatrix()");
  return m;
}

void free_dmatrix(FTYPE *m)
{
  free(m);
  m = NULL;
}
