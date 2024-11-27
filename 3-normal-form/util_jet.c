#include "util_jet.h"

#define mdeg _MAX_DEGREE_OF_JET_VARS_

double compute_validity_range(int mxdeg, MY_JET in, double alph, double log10tol)
{
  int j, k;
  MY_FLOAT tmp;
  double dtmp, val;
  int offset, *monomial_counts, *monomial_offsets;

  monomial_counts=MY_JET_FUN(monomial_counts)();
  monomial_offsets=MY_JET_FUN(monomial_offsets)();


  if (mxdeg > mdeg) mxdeg = mdeg;

  InitMyFloat(tmp);
  for (val=0.0, k = mxdeg-1; k <= mxdeg; ++k)
    {
      offset = monomial_offsets[k];

      for (j = 0; j < monomial_counts[k]; ++j)
        {
          fabsMyFloatA(tmp, MY_JET_DATA(in, offset+j));
          MyFloatToDouble(dtmp, tmp);
          if (dtmp > 0.0e0) {
              dtmp = pow(10, (log10tol - log10(dtmp))/(alph*k));
              if (dtmp > val) val = dtmp;
            }
        }
    }
  ClearMyFloat(tmp);
  return val;
}


