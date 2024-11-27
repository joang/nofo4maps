#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "poinca.h"

/**
 * @brief local and global variables
 *    ham: Hamiltonian value
 *    n: dim phase space
 *    jx0: Jet at the fixed point
 *    jperiod: Jet at the period
 */
double ham;
static int n;
static MY_JET *jx0, jperiod;
static clock_t cputime;
static FILE *file=NULL;

#define alloc(n,x) (x)= (__typeof__(x)) malloc((n)*sizeof(*(x)))
#define dealloc(x) if (x) free(x),(x)=NULL


/**
 * @brief read_input
 *    Reads the input from the file fixed point compute at the step 1
 *    It updates to the global variables:
 *      n: dimension of the phase space
 *      ham: Hamiltonian value
 *      period: period of the periodic orbit (since we are in spatial section)
 *      x: fixed point
 *
 * @param[in] fname filename
 * @return 0 if success
 */
static int read_input(const char *fname)
{
  int i, l;

  file = fopen(fname, "r");
  if (file==NULL) { printf("%s cannot be opened\n", fname); exit(1); }

  // read phase space
  l = fscanf(file, "%*[^=]%*c%d ", &n);
  if (l != 1) {printf("%s fails for n\n", fname); exit(2);}

  // read hamiltonian value
  l = fscanf(file, "%*[^=]%*c%lf ", &ham);
  if (l != 1) {printf("%s:%d failed reading %s\n", __FILE__, __LINE__, fname); exit(2);}

  /* Init jet of periods,
   * set it to zero, and
   * read period to save it in the 0th term
   */
  InitJet(jperiod); MY_JET_FUN(set_si)(jperiod, 0);
  l = fscanf(file, "%*[^=]%*c%lf ", &MY_JET_DATA(jperiod, 0));
  if (l != 1) {printf("%s:%d failed reading %s\n", __FILE__, __LINE__, fname); exit(2);}

  /* allocate vector of jets,
   * initialize them,
   * set them to zero,
   * read the fixed point coordinate and save them in the 0th term, and
   * set the symbols value to zero
   */
  alloc(n, jx0);
  for (i = 0; i < n; ++i)
    {
      InitJet(jx0[i]); MY_JET_FUN(set_si)(jx0[i], 0);

      l = fscanf(file, "%lf ", &MY_JET_DATA(jx0[i], 0));
      if (l != 1) {printf("%s:%d failed reading %s\n", __FILE__, __LINE__, fname); exit(2);}

      MakeMyFloatC(MY_JET_DATA(jx0[i], 1+i), "1.0", 1.0);
    }

  fclose(file); file=NULL;
  return 0;
}

#define print_usage() printf("%s <fin> <fout>\n", __FILE__)
int main(int argc, char *argv[])
{
#define fnamein argv[1]
#define fnameout argv[2]
  int i;

  if (argc < 2 || (argc > 1 && strcmp(argv[1], "-h")==0)) {print_usage(); exit(1);}

  InitUpJet();
  // read input
  read_input(fnamein);

  // compute the Poincare map with jets
  cputime = clock();
  poinca_jets(jx0, jperiod);
  cputime = clock() - cputime;

  // save the result into a file
  if (argc > 1) {
    file = fopen(fnameout, "w");
    if (file==NULL) {printf("%s cannot be opened\n", fnameout); exit(1);}
  } else { file = stdout; }


  fprintf(file, "# n = %d\n", n);
  fprintf(file, "# deg = %d\n", GetJetVarDegree());
  for (i = 0; i < n; ++i) {
      OutputJet2File(file, "% .15e\n", jx0[i]);
      fprintf(file, "\n");
    }

  fprintf(file, "\n# period = \n");
  OutputJet2File(file, "% .15e\n", jperiod);

  fprintf(file, "# cputime = % .10e sec\n", ((double) cputime) / CLOCKS_PER_SEC);


  if (argc > 1) {
      fclose(file); file=NULL;
      printf("Saved in %s\n", fnameout);
    }

  // deallocate some memory
  for (i = 0; i < n; ++i) {ClearJet(jx0[i]);} dealloc(jx0);
  ClearJet(jperiod);
  ClearUpJet();
  return 0;
}
