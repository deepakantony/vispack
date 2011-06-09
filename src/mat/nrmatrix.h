#include "nrutil.h"
// you will need to put a declaration for each c routine you plan to use
void gaussj(float **a, int n, float **b, int m);
void svdcmp(float **a, int m, int n, float w[], float **v);
float pythag(float a, float b);
void ludcmp(float **a, int n, int *indx, float *d);
void qrdcmp(float **a, int n, float *c, float *d, int *sing);
