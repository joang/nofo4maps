#include <stdio.h>
#include <stdlib.h>
#include "nofo.h"

static MY_JET *c=NULL, *G=NULL;
static int ndim=0;
#define nsym _NUMBER_OF_MAX_SYMBOLS_
#define mdeg _MAX_DEGREE_OF_JET_VARS_
#define OFMT "% .15e % .15e " // output printing format
#define IFMT "%lf %lf "       // input printing format


extern int is_unavoidable_resonant(int i, int im[]);

void nofo_init(int n)
{
  int i;
  if (ndim) {
      fprintf(stderr, "%s:%d nofo already initialized\n", __FILE__, __LINE__);
      exit(40);
    }
  ndim = n;

  InitUpJet();
  alloc(n, c); for (i = 0; i < n; ++i) {InitJet(c[i]);}
  alloc(n, G); for (i = 0; i < n; ++i) {InitJet(G[i]);}
}

void nofo_clear(void)
{
  int i;

  for (i = 0; i < ndim; ++i) {ClearJet(c[i]);} dealloc(c);
  for (i = 0; i < ndim; ++i) {ClearJet(G[i]);} dealloc(G);
  ClearUpJet();
  ndim=0;
}

static int is_resonant_monomial(MY_FLOAT out[1], int n, int i, MY_FLOAT vaps[], int im[], double log10tol)
{
  int j;
  MY_FLOAT tmp;
  double dtmp;

  InitMyFloat(tmp);

  MakeMyFloatC(*out, "1.0", 1.0);
  for (j = 0; j < n; ++j) if (im[j])
    {
      ExponentiateMyFloatA(tmp, vaps[j], im[j]);
      MultiplyMyFloatA(tmp, *out, tmp);

      AssignMyFloat(*out, tmp);
    }
  SubtractMyFloatA(*out, *out, vaps[i]);

  fabsMyFloatA(tmp, *out);
  MyFloatToDouble(dtmp, tmp);
  dtmp = log10(dtmp);

  if (dtmp < log10tol) {j=1;}
  else {j=0;}

  ClearMyFloat(tmp);
  return j;
}




/**
 * @brief compute_nofo
 *    Computes normal for from a given formal power series
 *    It outputs the normal form, and the changes
 *
 * @param n
 * @param mxdeg
 * @param F
 * @param changes
 * @param scaling
 * @param finfo
 * @param log10tol
 * @param info
 * @return
 */
int compute_nofo(int n, int mxdeg, MY_JET F[], MY_JET changes[],
                 double scaling, FILE *finfo, double log10tol, int info[1])
{
  int i, j, l, *monomial_counts=NULL;
  MY_FLOAT zero, lambda, cwork[n];
  MY_FLOAT vaps[n], veps[n*n];
  int k, nch, nv, ord;
  int im[nsym];
  typeof(F[0]) hf, hc;

  *info=0;

  // init memory
  InitMyFloat(zero); MakeMyFloatC(zero, "0.0", 0.0);
  InitMyFloat(lambda);
  for (i = 0; i < n; ++i)   {InitMyFloat(vaps[i]);}
  for (i = 0; i < n*n; ++i) {InitMyFloat(veps[i]);}
  for (i = 0; i < n; ++i)   {InitMyFloat(cwork[i]);}

  // get the number of monomials for each homogenous polynomial
  monomial_counts = MY_JET_FUN(monomial_counts)();

  // 0th order: F0 <- G - x0
  for (i = 0; i < n; ++i)
    {
      if (changes) {AssignFloatToJet(changes[i], MY_JET_DATA(F[i], 0));}
      AssignMyFloat(MY_JET_DATA(F[i], 0), zero);
    }

  // 1st order: F1 <- c1^{-1} F0 c1  with c1(s) = V s matrix of veps */
  { /* 2-by-2 eigenvalue problem */
    for (i = 0; i < n; ++i) for (j = 0; j < n; ++j)
      {AssignMyFloat(aij(veps,n, i,j), MY_JET_DATA(F[i], 1+j));}

    eigen2d(n, veps, n, vaps, veps, n);

    for (i = 0; i < n; ++i)
      {
        AssignFloatToJet(c[i], zero);
        for (j = 0; j < n; ++j)
          {
            MultiplyMyFloatD(aij(veps,n, i,j), aij(veps,n, i,j), scaling);
            AssignMyFloat(MY_JET_DATA(c[i], 1+j), aij(veps,n, i,j));

            if (changes) { AssignMyFloat(MY_JET_DATA(changes[i], 1+j), MY_JET_DATA(c[i], 1+j)); }
          }
      }

    for (i = 0; i < n; ++i) {MY_JET_FUN(compo)(G[i], F[i], c);}

    // the inverse by c requires to solve (in x) the linear system VEPS x = F
    for (l = 0; l <= mxdeg; ++l)
      { // we loop on each homogenous polynomial

        hf = F[0]+l;
        nch = monomial_counts[hf->deg]; // get the number of monomials of order l
        for (k = 0; k < nch; ++k)
          {
            for (i = 0; i < n; ++i)
              { // save the coefficient on the auxiliary vector cwork
                hf = G[i]+l;
                AssignMyFloat(cwork[i], hf->coef[k]);
              }

            // linear solving
            solve2d(n, 1, veps, n, cwork, n, log10tol, &i);
            if (i) {fprintf(stderr, "%s:%d solve2d(%d) error\n", __FILE__, __LINE__, i); exit(20);}

            for (i = 0; i < n; ++i)
              { // assign the solution back on the coefficient
                hf = F[i]+l;
                AssignMyFloat(hf->coef[k], cwork[i]);
              }
          }
      }
  }
  // we ensure that entries outside the diagonal are exactly zero
  for (i = 0; i < n; ++i) for (j = 0; j < n; ++j) if (i != j)
    {AssignMyFloat(MY_JET_DATA(F[i], 1+j), zero);}

  // now we apply the procedure on all the other coefficient
  for (l = 2; l <= mxdeg; ++l)
    {
      for (i = 0; i < n; ++i)
        { // define the change of coordinate
          AssignFloatToJet(c[i], zero);
          MakeMyFloatC(MY_JET_DATA(c[i], 1+i), "1.0", 1.0);

          hc = c[i]+l;
          hf = F[i]+l;

          nv = hc->nsymb;
          ord = hc->deg;
          nch = monomial_counts[ord];
          im[0]=ord; for (j = 1; j < nv; j++) im[j]=0;
          for (k = 0; k < nch; ++k)
            {
              // ressonance treatment
              if (is_unavoidable_resonant(i, im))
                {
                  if (finfo) {
                      fprintf(finfo, "im="); for (j = 0; j < nv; j++) {fprintf(finfo, "%d ", im[j]);}
                      fprintf(finfo, " UNavoidable resonance %d %d %d\n", l, i, k);
                    }
                }
              else if (is_resonant_monomial(&lambda, n, i, vaps, im, log10tol))
                {
                  if (finfo) {
                      fprintf(finfo, "im="); for (j = 0; j < nv; j++) {fprintf(finfo, "%d ", im[j]);}
                      fprintf(finfo, " resonance %d %d %d\n", l, i, k);
                    }
                }
              else
                {
                  DivideMyFloatA(hc->coef[k],hf->coef[k], lambda);

                  if (changes) {AssignMyFloat(changes[i][l].coef[k], hc->coef[k]);}
                }

              MY_JET_FUN(genidx)(im,nv);
            }
        }

      // composition and reverse transformation
      for (i = 0; i < n; ++i) {MY_JET_FUN(compo)(G[i], F[i], c);}
      MY_JET_FUN(algT)(F, G, c);

      // we set to zero what it was supposed to be zero
      // to prevent sporious error propagation
      nch = monomial_counts[l];
      for (i = 0; i < n; ++i) {
          hf = F[i]+l;

          nv = hc->nsymb;
          ord = hc->deg;
          nch = monomial_counts[ord];
          im[0]=ord; for (j = 1; j < nv; j++) {im[j]=0;}
          for (k = 0; k < nch; ++k)
            {
              if (is_unavoidable_resonant(i, im) == 0 &&
                  is_resonant_monomial(&lambda, n, i, vaps, im, log10tol) == 0)
                {
                  AssignMyFloat(hf->coef[k], zero);
                }

              MY_JET_FUN(genidx)(im,nv);
            }
        }

      if (finfo) {
          fprintf(finfo, "Order %d done\n", l);
        }
    }

  // deallocate memory
  for (i = 0; i < n; ++i) {ClearMyFloat(cwork[i]);}
  for (i = 0; i < n*n; ++i) {ClearMyFloat(veps[i]);}
  for (i = 0; i < n; ++i) {ClearMyFloat(vaps[i]);}
  ClearMyFloat(lambda);
  ClearMyFloat(zero);
  return 0;
}

void save_nofo(FILE *file, int mxdeg, int i, MY_JET F, double log10tol)
{
  MY_FLOAT z;
  int nnz, j, l, im[nsym];
  int k, nch, nv, ord;
  typeof(F) hf;
  int *monomial_counts;
  double dtmp;

  InitMyFloat(z);

  if (mxdeg > mdeg) mxdeg = mdeg;

  // first count how many non zero (nnz) elements F has
  monomial_counts = MY_JET_FUN(monomial_counts)();
  for (nnz=0, l = 0; l <= mxdeg; ++l)
    {
      hf = F+l;

      nv = hf->nsymb; ord = hf->deg; nch = monomial_counts[ord];
      im[0]=ord; for (j = 1; j < nv; j++) {im[j]=0;}
      for (k = 0; k < nch; ++k)
        {
          if (is_unavoidable_resonant(i, im)) {
              fabsMyFloatA(z, hf->coef[k]);
              MyFloatToDouble(dtmp, z);
              if (log10(dtmp) > log10tol) {nnz++;}
            }

          MY_JET_FUN(genidx)(im,nv);
        }
    }
  fprintf(file, "# nnz = %d\n", nnz);

  // print
  for (l = 0; l <= mxdeg; ++l)
    {
      hf = F+l;

      nv = hf->nsymb;
      ord = hf->deg;
      nch = monomial_counts[ord];
      im[0]=ord; for (j = 1; j < nv; j++) {im[j]=0;}
      for (k = 0; k < nch; ++k)
        {
          if (is_unavoidable_resonant(i, im)) {
              fabsMyFloatA(z, hf->coef[k]);
              MyFloatToDouble(dtmp, z);
              if (log10(dtmp) > log10tol) {
                  nnz++;
                  AssignMyFloat(z, hf->coef[k]);

                  OutputMyFloat3(file, OFMT, z);
                  for (j = 0; j < nv; j++) fprintf(file, " %d", im[j]);
                  fprintf(file, "\n");
                }
            }

          MY_JET_FUN(genidx)(im,nv);
        }
    }
}

void read_nofo(FILE *file, MY_JET F)
{
  int nnz, i, k, l, idx[nsym];
  double re, im;

  l = fscanf(file, "%*[^=]%*c%d ", &nnz);
  if (l != 1) {printf("%s:%d I/O error for nnnz\n", __FILE__, __LINE__); exit(2);}

  MY_JET_FUN(set_si)(F, 0);
  for (k = 0; k < nnz; ++k)
    {
      l = fscanf(file, IFMT, &re, &im);
      if (l != 2) {printf("%s:%d I/O error\n", __FILE__, __LINE__); exit(2);}

      for (i = 0; i < nsym; ++i)
        {
          l = fscanf(file, "%d ", &idx[i]);
          if (l != 1) {printf("%s:%d I/O error\n", __FILE__, __LINE__); exit(2);}
        }
      AssignMyFloat(*MY_JET_FUN(get_coef)(F, idx), CMPLX(re, im));
    }
}

