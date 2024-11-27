#ifndef NOFO_H
#define NOFO_H

#include "la2d.h"

void nofo_init(int n);
void nofo_clear(void);

int compute_nofo(int n, int mxdeg, MY_JET G[], MY_JET changes[],
                 double scaling, FILE *finfo, double log10tol, int info[1]);



void save_nofo(FILE *file, int mxdeg, int i, MY_JET F, double log10tol);
void read_nofo(FILE *file, MY_JET F); // unused

#endif // NOFO_H
