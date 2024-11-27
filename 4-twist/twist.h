#ifndef TWIST_H
#define TWIST_H

#include <stdlib.h>
#include "ode.h"

#ifndef alloc
#define alloc(n,x) (x)= (__typeof__(x)) malloc((n)*sizeof(*(x)))
#endif
#ifndef dealloc
#define dealloc(x) if (x) free(x),(x)=NULL
#endif


int compute_twist(int n, MY_JET F[], MY_JET twists[]);

int read_nofo4twist(FILE *file, int index[], MY_JET F);

void save_twist(FILE *file, int mxdeg, MY_JET F, double log10tol);


#endif // TWIST_H
