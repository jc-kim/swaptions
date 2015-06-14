#include "HJM_type.h"
#include <string>

void     nrerror(std::string error_text);
int      *ivector(long nl, long nh);
void     free_ivector(int *v, long nl, long nh);
FTYPE   *dvector( long nl, long nh );
void     free_dvector( FTYPE *v, long nl, long nh );
FTYPE   **dmatrix( long nrl, long nrh, long ncl, long nch );
void     free_dmatrix( FTYPE **m, long nrl, long nrh, long ncl, long nch );
