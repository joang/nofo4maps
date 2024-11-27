#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "twist.h"
#include "util_jet.h"

/**
 * @brief local and global variables
 *    n: dim phase space
 *    mxdeg: max. degree
 *    log10tol: log10 tolerance
 *    scaing: scaling factor
 */
static int n, mxdeg=1;
static double log10tol, scaling;
static MY_JET *F=NULL;
static clock_t cputime;
static FILE *file=NULL;

/**
 * @brief read_input_nofo
 *    Reads the input from previous normal form computation
 *    More precisely:
 *      n: number of coordinates
 *      mxdeg: max. degree
 *      F: vector of normal forms (size n)
 *
 * @param fname
 * @return
 */
static int read_input_nofo(const char *fname)
{
  int i, l, k, nnz, *idx;
  double re, im;

  file = fopen(fname, "r");
  if (file==NULL) { printf("%s cannot be opened\n", fname); exit(1); }


  // read nofo dimension
  l = fscanf(file, "%*[^=]%*c%d ", &n);
  if (l != 1) {printf("%s:%d I/O error\n", __FILE__, __LINE__); exit(2);}

  n/= 2;
  alloc(n, idx);

  // read max. degree
  l = fscanf(file, "%*[^=]%*c%d ", &mxdeg);
  if (l != 1) {printf("%s fails for mxdeg\n", fname); exit(2);}

  // read log10tol
  l = fscanf(file, "%*[^(]%*c%lf ", &log10tol);
  if (l != 1) {printf("%s:%d I/O error\n", __FILE__, __LINE__); exit(2);}

  // read scaling factor
  l = fscanf(file, "%*[^=]%*c%lf ", &scaling);
  if (l != 1) {printf("%s:%d I/O error\n", __FILE__, __LINE__); exit(2);}

  if (F==NULL) {
      alloc(n, F);
      for (i = 0; i < n; ++i) {InitJet(F[i]);}
    }
  idx[0] = 1;
  for (i = 0; i < n; ++i) {read_nofo4twist(file, idx, F[i]);}

  fclose(file); file=NULL;
  dealloc(idx);
  return 0;
}

#define print_usage() printf("%s <fin-nofo> <fout-twist>\n", __FILE__)
int main(int argc, char *argv[])
{
#define fnamein_nofo argv[1]
#define fnameout_twist argv[2]
  int i;
  double smax;

  if (argc < 2 || (argc > 1 && strcmp(argv[1], "-h")==0)) {print_usage(); exit(1);}

  InitUpJet();
  // read inputs
  read_input_nofo(fnamein_nofo);


  cputime = clock();

  // computation of twist
  compute_twist(n, F, F);

  cputime = clock() - cputime;


  // save the result
  file = fopen(fnameout_twist, "w");
  if (file==NULL) {printf("%s cannot be opened\n", fnameout_twist); exit(1);}


  for (i = 0; i < n; ++i)
    {
      smax = compute_validity_range(mxdeg, F[i], 2.0, log10tol);
      fprintf(file, "# smax(%g) = % .15e\n", log10tol, smax);
      save_twist(file, mxdeg, F[i], log10tol);
      fprintf(file, "\n");
    }
  fprintf(file, "# cputime = % .10e sec\n", ((double) cputime) / CLOCKS_PER_SEC);

  fclose(file); file=NULL;
  printf("Saved in %s\n", fnameout_twist);



  // deallocate some memory
  for (i = 0; i < n; ++i) {ClearJet(F[i]);} dealloc(F);
  ClearUpJet();
  return 0;
}
