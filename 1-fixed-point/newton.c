#include "newton.h"
#include "la2d.h"
#include "ode.h"

extern const int N;
extern MY_FLOAT ham;

#define REAL MY_FLOAT

#define NS _NUMBER_OF_STATE_VARS_
#define NJ _NUMBER_OF_JET_VARS_
#define NSYMB _NUMBER_OF_MAX_SYMBOLS_
#define log10abs -16
#define log10rel -16


/**
 * @brief poincare
 *  It computes the Poincare map of spatial section
 *
 * @param[in] x 2D vector
 * @param[out] f 2D vector
 * @param[out] df if no NULL, monodromy matrix (in col-major)
 * @return Flying time (after 2 crossing sections)
 */
REAL poincare(REAL x[N], REAL f[N], REAL df[N*N])
{
  static int init=0;
  static MY_FLOAT orbit[NS], t, ht;
  static MY_JET jet_orbit[NJ];
  int i, ord;
  MY_JET **jet_px;

  if (init==0) {
      InitUpJet();

      InitMyFloat(t);
      InitMyFloat(ht);
      for (i = 0; i < NS; ++i) {InitMyFloat(orbit[i]);}
      for (i = 0; i < NJ; ++i) {InitJet(jet_orbit[i]);}

      init=1;
    }


  // initial time
  MakeMyFloatC(t, "0.0", 0.0); // initial time

  /* (y,py) to (x,y,px,py) where:
   *    x=0,
   *    px is deduced from the energy level
   */
  MakeMyFloatC(orbit[0], "0.0", 0.0);
  AssignMyFloat(orbit[1], x[0]);
  MakeMyFloatC(orbit[2], "0.0", 0.0);
  AssignMyFloat(orbit[3], x[1]);

  for (i = 0; i < NJ; ++i)
    {AssignFloatToJet(jet_orbit[i], orbit[i]);}
  MakeMyFloatC(MY_JET_DATA(jet_orbit[1], 1+0), "1.0", 1.0);
  MakeMyFloatC(MY_JET_DATA(jet_orbit[3], 1+1), "1.0", 1.0);

  pxval(t, orbit, NULL, jet_orbit, &jet_px);

  MultiplyMyFloatSI(orbit[2], ham, 2);
  AddJetFloatA(jet_px[0][0], jet_px[0][0], orbit[2]);
  if (MyFloatA_LT_B(MY_JET_DATA(jet_px[0][0],0), t)) return 0; // out of domain

  sqrtJetA(jet_orbit[2], jet_px[0][0]);

  AssignJetToFloat(orbit[2], jet_orbit[2]);

  // Integration until first cross section
  do {
      taylor_step_ode(&t, orbit, +1, 2, log10abs, log10rel, NULL, NULL,&ord, jet_orbit);
    } while (orbit[0] > 0);
  // Integration until second cross section
  do {
      taylor_step_ode(&t, orbit, +1, 2, log10abs, log10rel, NULL, NULL,&ord, jet_orbit);
    } while (orbit[0] < 0);

  // Newton to compute flying time
  while (fabs(orbit[0]) > 1e-14)
    {
      DivideMyFloatA(ht, orbit[0], orbit[2]);
      NegateMyFloatA(ht, ht);
      taylor_step_ode(&t, orbit, +1, 0, log10abs, log10rel, NULL, &ht,&ord, jet_orbit);
    }

  // Assignment of the value
  if (f) {
      AssignMyFloat(f[0], orbit[1]);
      AssignMyFloat(f[1], orbit[3]);
    }

  if (df == NULL) return t;

  // Projection correction to the section, 1D projection
  MY_FLOAT **torbit;
  torbit = taylor_coefficients_ode_A(t, orbit, 1, 1, jet_orbit, NULL);

  ht = -MY_JET_DATA(jet_orbit[0], 1+0)/torbit[0][1];
  aij(df,N, 0,0) = MY_JET_DATA(jet_orbit[1], 1+0) + ht * torbit[1][1]; //[0][0]
  aij(df,N, 1,0) = MY_JET_DATA(jet_orbit[3], 1+0) + ht * torbit[3][1]; //[1][0]

  ht = -MY_JET_DATA(jet_orbit[0], 1+1)/torbit[0][1];
  aij(df,N, 0,1) = MY_JET_DATA(jet_orbit[1], 1+1) + ht * torbit[1][1]; //[0][1]
  aij(df,N, 1,1) = MY_JET_DATA(jet_orbit[3], 1+1) + ht * torbit[3][1]; //[1][1]

  return t;
}

int newton(REAL x[N], REAL T[1], REAL df[N*N], REAL errout[1], REAL log10tol, FILE *file)
{
  static int init_flag=0;
  int i, iter=0, max_iter=15;
  REAL y[N];
  double log10err, log10errbef=10;

  if (file) fprintf(file, "# T log10error log10correction\n");

  for (iter = 0; iter < max_iter; ++iter)
    {
      *T = poincare(x, y, df);

      for (log10err=0.0, i = 0; i < N; ++i) {
          y[i]-= x[i];
          log10err+= y[i]*y[i];
        }
      log10err = log10(log10err)/2;

      if (file) fprintf(file, "% .15e % 3.3f ", *T, log10err);
      if (log10err < log10tol)
        {
          if (file) fprintf(file, " ok!\n\n");
          goto end;
        }

      for (i = 0; i < N; ++i) {aij(df,N, i,i)-= 1.0e0;}
      solve2d(N,1, df,N, y,N, log10tol, &i);
      if (i!=0) {
          if (file) fprintf(file, "Unsolved linear system at step %d\n", iter);
          exit(1);
        }

      for (log10err=0.0, i = 0; i < N; ++i) {log10err+= y[i]*y[i];}
      log10err = log10(log10err)/2;

      if (file) fprintf(file, "% 3.3f\n", log10err);

      for (i = 0; i < N; ++i) {x[i]-= y[i];}

      if (log10err > 8 || log10errbef < log10err*1)
        {
          if (file) fprintf(file, "Newton fails at step %d ", iter);
          break;
        }
      log10errbef = log10err;
    }
  if (file) fprintf(file, "\nko!\n\n");
  i = -1;
end:
  return i;
}

