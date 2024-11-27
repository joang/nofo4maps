#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "newton.h"
#include "la2d.h"

const int N=2;
#define LOG10TOL -14
double ham;

static FILE *file=NULL;

#define print_usage() printf("%s <fout>\n", __FILE__)
int main(int argc, char *argv[])
{
#define fnameout argv[1]
  int i, j;
  double x[N], df[N*N], T, err;

  if (argc > 1 && strcmp(argv[1], "-h")==0) {print_usage(); exit(1);}

  // initial guesses
  ham = 0.125e0;
  x[0]= 3.0140065033328262e-01;
  x[1]= 0.0;
  
  /*
   * Newton scheme to compute fixed point.
   *   Input: x initial guess,
   *          log10tolerance
   *   Return: ==0 if success
   *   Outputs: x solution,
   *            T period T,
   *            df monodromy matrix,
   *            err final error
   */
  i = newton(x, &T, df, &err, LOG10TOL, stdout);
  
  if (i < 0) {return 1;}


  if (argc > 1) {
    file = fopen(fnameout, "w");
    if (file==NULL) {printf("%s cannot be opened\n", fnameout); exit(1);}
  } else { file = stdout; }
  
  fprintf(file, "# phasedim = %d\n", N);
  fprintf(file, "# ham = % .15e\n", ham);
  fprintf(file, "# period = % .15e\n", T);
  for (i = 0; i < N; ++i)   fprintf(file, "% .15e\n", x[i]);
  
  for (i = 0; i < N; ++i) for (j = 0; j < N; ++j) fprintf(file, "% .15e\n", aij(df,N, i,j));


  if (argc > 1) {
      fclose(file); file=NULL;
      printf("Saved in %s\n", fnameout);
    }
  return 0;
}
