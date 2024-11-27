#ifndef LA2D_H
#define LA2D_H

#include <stdlib.h>
#include "ode.h"

#define aij(a,lda, i,j) (a)[(j)*(lda)+(i)] // col-major

#define alloc(n,x) (x)= (__typeof__(x)) malloc((n)*sizeof(*(x)))
#define dealloc(x) if (x) free(x),(x)=NULL

int eigen2d(int n, MY_FLOAT a[], int lda, MY_FLOAT vaps[], MY_FLOAT veps[], int ldv);
void solve2d(int n, int ns, MY_FLOAT a[], int lda, MY_FLOAT b[], int ldb, double log10tol, int info[]);

#endif // LA2D_H
