#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "la2d.h"
#include "nofo.h"
#include "util_jet.h"

/**
 * @brief local and global variables
 *    ham: Hamiltonian value
 *    n: dim phase space
 *    jx0: Jet at the fixed point
 *    jperiod: Jet at the period
 */
double ham;
static int n;
static int mxdeg=1;
static MY_JET *jx0;
//static MY_JET jperiod; // not used
static clock_t cputime;
static FILE *file=NULL;

/**
 * @brief read_input_jet_pfix
 *    Reads the input from the file jet at the fixed point computed at the step 1&2
 *    It updates to the global variables:
 *      mxdeg: max. degree
 *      jx0: jet at fixed point
 *      jperiod: jet of the period of the periodic orbit
 *               (since we are in spatial section)
 *
 * @param[in] fname filename
 * @return 0 if success
 */
static int read_input_jet_pfix(const char *fname)
{
  int i, l;
  unsigned int k;
  double re;

  file = fopen(fname, "r");
  if (file==NULL) { printf("%s cannot be opened\n", fname); exit(1); }

  // ndim
  l = fscanf(file, "%*[^=]%*c%d ", &n);
  if (l != 1) {printf("%s fails for n\n", fname); exit(2);}

  // read max. degree
  l = fscanf(file, "%*[^=]%*c%d ", &mxdeg);
  if (l != 1) {printf("%s fails for mxdeg\n", fname); exit(2);}

  // allocate vector of jets and read their coefficients
  alloc(n, jx0);
  for (i = 0; i < n; ++i)
    {
      InitJet(jx0[i]);
      for (k = 0; k < _JET_COEFFICIENTS_COUNT_TOTAL_; ++k){
          l = fscanf(file, "%lf ", &re);
          if (l != 1) {printf("%s:%d failed reading %s\n", __FILE__, __LINE__, fname); exit(2);}

          MY_JET_DATA(jx0[i], k) = CMPLX(re, 0.0e0);
        }
    }

  fscanf(file, "%*[^=]%*c ");

  // allocate jet of period and read its coefficients
  /* // not needed
  InitJet(jperiod);
  for (k = 0; k < _JET_COEFFICIENTS_COUNT_TOTAL_; ++k){
      l = fscanf(file, "%lf ", &re);
      if (l != 1) {printf("%s:%d failed reading %s\n", __FILE__, __LINE__, fname); exit(2);}

      MY_JET_DATA(jperiod, k) = CMPLX(re, 0.0e0);
    }
  */

  fclose(file); file=NULL;
  return 0;
}


int is_unavoidable_resonant(int i, int im[])
{
  return (i==0 && im[0]-im[1]==1) || (i==1 && im[0]-im[1]==-1);
}

#define print_usage() printf("%s <fin-jet-pfix> <fout>\n", __FILE__)
int main(int argc, char *argv[])
{
#define fnamein_jet_pfix argv[1]
#define fnameout argv[2]
#define LOG10TOL -15.0
  int i;
  MY_JET *G=NULL, *changes=NULL;
  double scaling;

  if (argc < 3 || (argc > 1 && strcmp(argv[1], "-h")==0)) {print_usage(); exit(1);}

  InitUpJet();
  // read inputs
  read_input_jet_pfix(fnamein_jet_pfix);


  // allocate memory
  nofo_init(n);
  alloc(n, changes); for (i = 0; i < n; ++i) {InitJet(changes[i]);}
  alloc(n, G); for (i = 0; i < n; ++i) {InitJet(G[i]);}

  // init the output
  for (i = 0; i < n; ++i) {AssignJetToJet(G[i], jx0[i]);}

  cputime = clock();

  // computation of normal form with no scaling
  compute_nofo(n, mxdeg, G, changes, 1.0, NULL, LOG10TOL, &i);

  // computation of a scaling factor
  for (scaling = 1.0, i = 0; i < n; ++i)
    {scaling = fmin(scaling, compute_validity_range(mxdeg, G[i], 1.0, LOG10TOL));}
  printf("scaling(%g) = % .15e\n", LOG10TOL, scaling);

  // recompute the normal form now with scaling
  for (i = 0; i < n; ++i) {AssignJetToJet(G[i], jx0[i]);}
  compute_nofo(n, mxdeg, G, changes, scaling, stdout, LOG10TOL, &i);

  cputime = clock() - cputime;



  // save the results
  file = fopen(fnameout, "w");
  if (file==NULL) {printf("%s cannot be opened\n", fnameout); exit(1);}

  fprintf(file, "# n = %d\n", n);
  fprintf(file, "# mxdeg = %d\n", GetJetVarDegree());
  fprintf(file, "# scal(%g) = % .15e\n", LOG10TOL, scaling);
  for (i = 0; i < n; ++i) {
      save_nofo(file, mxdeg, i, G[i], LOG10TOL);
      fprintf(file, "\n");
    }
  fprintf(file, "# changes =\n");
  for (i = 0; i < n; ++i) {
      OutputJet2File(file, "% .15e % .15e\n", changes[i]);
      fprintf(file, "\n");
    }

  fprintf(file, "# cputime = % .10e sec\n", ((double) cputime) / CLOCKS_PER_SEC);
  fclose(file); file=NULL;
  printf("Saved in %s\n", fnameout);



  // deallocate some memory
  for (i = 0; i < n; ++i) {ClearJet(changes[i]);} dealloc(changes);
  for (i = 0; i < n; ++i) {ClearJet(G[i]);} dealloc(G);

  for (i = 0; i < n; ++i) {ClearJet(jx0[i]);} dealloc(jx0);
//  ClearJet(jperiod);

  ClearUpJet();
  nofo_clear();
  return 0;
}
