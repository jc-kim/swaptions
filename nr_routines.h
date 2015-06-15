#include "HJM_type.h"
#include <string>

void nrerror(std::string error_text);
FTYPE* dvector(long length);
void free_dvector(FTYPE *v);
FTYPE** dmatrix(long rows, long cols);
void free_dmatrix(FTYPE **m);
