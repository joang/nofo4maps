#include "twist.h"

#define nsym _NUMBER_OF_MAX_SYMBOLS_
#define mdeg _MAX_DEGREE_OF_JET_VARS_
#define OFMT "% .15e % .15e " // output printing format
#define IFMT "%lf %lf "       // input printing format

int read_nofo4twist(FILE *file, int index[], MY_JET F)
{
  int nnz, i, k, l, idx[2*nsym];
  double re, im;

  l = fscanf(file, "%*[^=]%*c%d ", &nnz);
  if (l != 1) {printf("%s:%d I/O error for nnnz\n", __FILE__, __LINE__); exit(2);}

  MY_JET_FUN(set_si)(F, 0);
  for (k = 0; k < nnz; ++k)
    {
      l = fscanf(file, IFMT, &re, &im);
      if (l != 2) {printf("%s:%d I/O error\n", __FILE__, __LINE__); exit(2);}

      for (i = 0; i < 2*nsym; ++i)
        {
          l = fscanf(file, "%d ", &idx[i]);
          if (l != 1) {printf("%s:%d I/O error\n", __FILE__, __LINE__); exit(2);}
        }
      for (i = 0; i < nsym; ++i) idx[i] = idx[index[i]];
      AssignMyFloat(*MY_JET_FUN(get_coef)(F, idx), CMPLX(re, im));
    }
  return 0;
}

int compute_twist(int n, MY_JET F[], MY_JET twist[])
{
  int i;
  MY_JET jtmp;
  MY_FLOAT a;

  InitJet(jtmp);
  InitMyFloat(a);

  for (i = 0; i < n; ++i)
    {
//      OutputJet2File(stdout, "% .10e % .10e\n ", F[i]);
      AssignMyFloat(a, MY_JET_DATA(F[i],0));
      MY_JET_FUN(div2_myfloat)(twist[i], F[i], a);

      MY_JET_FUN(set_log)(jtmp, twist[i]);
      MY_JET_FUN(div2_myfloat)(twist[i], jtmp, _Complex_I);

      AssignMyFloat(MY_JET_DATA(twist[i],0), atan2(cimag(a), creal(a)));
//      OutputJet2File(stdout, "% .10e % .10e\n ", twist[i]);
    }

  ClearMyFloat(a);
  ClearJet(jtmp);
}

void save_twist(FILE *file, int mxdeg, MY_JET F, double log10tol)
{
  MY_FLOAT z;
  int nnz, j, l, im[nsym];
  int k, nch, nv, ord;
  typeof(F) hf;
  int *monomial_counts;
  double dtmp;

  if (mxdeg > mdeg) mxdeg = mdeg;

  InitMyFloat(z);

  // first count how many non zero (nnz) elements F has
  monomial_counts = MY_JET_FUN(monomial_counts)();
  for (nnz=0, l = 0; l <= mxdeg; ++l)
    {
      hf = F+l;

      nv = hf->nsymb; ord = hf->deg; nch = monomial_counts[ord];
      for (j = 0; j < nv; j++) {im[j]=0;}
      im[0]=ord;
      for (k = 0; k < nch; ++k)
        {
          fabsMyFloatA(z, hf->coef[k]);
          MyFloatToDouble(dtmp, z);
          if (log10(dtmp) > log10tol) {nnz++;}

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
      for (j = 0; j < nv; j++) {im[j]=0;}
      im[0]=ord;
      for (k = 0; k < nch; ++k)
        {
          fabsMyFloatA(z, hf->coef[k]);
          MyFloatToDouble(dtmp, z);
          if (log10(dtmp) > log10tol) {
              nnz++;
              AssignMyFloat(z, hf->coef[k]);

              OutputMyFloat3(file, OFMT, z);
              for (j = 0; j < nv; j++) fprintf(file, " %d", im[j]);
              fprintf(file, "\n");
            }

          MY_JET_FUN(genidx)(im,nv);
        }
    }
}


