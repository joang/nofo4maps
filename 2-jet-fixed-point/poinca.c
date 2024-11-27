#include "poinca.h"

extern MY_FLOAT ham;

#define NS _NUMBER_OF_STATE_VARS_
#define NJ _NUMBER_OF_JET_VARS_
#define log10abs -16
#define log10rel -16

/**
 * @brief local global variables
 */
static int init=0;
static MY_FLOAT orbit[NS], ti;
static MY_JET jet_orbit[NJ];

/**
 * @brief projection
 *    Projects the infinitesimal terms into a section {x=0}
 *
 * @param[in] n: dimension of the system
 * @param[in] deg: max degree (set it by the model file)
 * @param[out] Tjet: infinitesimal corrections of the period jet
 * @param[in] x: 0th order
 * @param[in,out] y: infinitesimal projections
 * @param hjet:  auxiliary jet
 */
static void projection(int n, int deg, MY_JET Tjet,
                       MY_FLOAT x[], MY_JET y[], MY_JET hjet)
{
  int i, j, k, l;
  MY_FLOAT **xjet=NULL;
  MY_JET **yjet=NULL;
  extern MY_FLOAT ti;
  int *monomial_counts, offset;

  // obtain the number homogenous terms for each order
  monomial_counts=MY_JET_FUN(monomial_counts)();

  offset=monomial_counts[0]; // always 1 since the order=0 is the 0th term

  for (l = 1; l <= deg; ++l) // for each degre
    {
      // we compute the vector field with jets y to get yjet
      // note that we compute it up to order deg (no up to l)
      xjet = taylor_coefficients_ode_A(ti, x, deg, 0, y, &yjet);

      MY_JET_FUN(set_si)(hjet, 0);
      for (j = 0; j < monomial_counts[l]; ++j)
        {
          // project with the normal direction of the section {x=0}
          DivideMyFloatA(MY_JET_DATA(hjet, offset+j), MY_JET_DATA(y[0], offset+j),xjet[0][1]);
          NegateMyFloatA(MY_JET_DATA(hjet, offset+j), MY_JET_DATA(hjet, offset+j));

          // set the value for the infinitesimal time correction
          AssignMyFloat(MY_JET_DATA(Tjet, offset+j), MY_JET_DATA(hjet, offset+j));
        }

      // horner method with the updated infinitesimal time correction
      for (i = 0; i < n; ++i)
        {
          AssignJetToJet(y[i], yjet[i][deg]);
          for (k = deg-1; k >= 0; --k)
            {
              // hjet is just an auxiliary jet for this step
              MultiplyJetJetA(yjet[i][deg], y[i], hjet);
              AddJetJetA(y[i], yjet[i][deg], yjet[i][k]);
            }
        }

      // update the offset to move to the next order
      offset+= monomial_counts[l];
    }
}

/**
 * @brief poinca_jets
 *    Given initial jx and jperiod it computes
 *    the derivatives of the Poincare map
 *
 * @param[in,out] jx
 * @param[in,out] jperiod
 * @return 0 if success
 */
int poinca_jets(MY_JET jx[], MY_JET jperiod)
{
  int i;
  MY_FLOAT te;
  MY_JET **jet_px;

  // initializations
  if (init==0) {
      InitUpJet();

      InitMyFloat(t);
      for (i = 0; i < NS; ++i) {InitMyFloat(orbit[i]);}
      for (i = 0; i < NJ; ++i) {InitJet(jet_orbit[i]);}

      init=1;
    }
  // the final integration time is the 0th order of the jet jperiod
  InitMyFloat(te);
  AssignJetToFloat(te, jperiod);

  // set the initial jet_orbit by the input jx
  MakeMyFloatC(ti, "0.0", 0.0);

  AssignFloatToJet(jet_orbit[0], ti);
  AssignJetToJet(jet_orbit[1], jx[0]);
  AssignFloatToJet(jet_orbit[2], ti);
  AssignJetToJet(jet_orbit[3], jx[1]);

  // set the orbit (with no jets) from the jet_orbit
  for (i = 0; i < NS; ++i) {AssignJetToFloat(orbit[i], jet_orbit[i]);}

  // evaluate pxval with jets
  pxval(ti, orbit, NULL, jet_orbit, &jet_px);

  MultiplyMyFloatSI(orbit[2], ham, 2);
  AddJetFloatA(jet_px[0][0], jet_px[0][0], orbit[2]);

  if (MyFloatA_LT_B(MY_JET_DATA(jet_px[0][0], 0), ti)) {
      fprintf(stderr, "%s:%d outside domain\n", __FILE__, __LINE__);
      return 1;
    }
  sqrtJetA(jet_orbit[2], jet_px[0][0]);  // final px value
  AssignJetToFloat(orbit[2], jet_orbit[2]);

  while (taylor_step_ode(&ti, orbit, +1, 2,
                         log10abs, log10rel,
                         &te, NULL, NULL,
                         jet_orbit) != 1) {}

  // project the jet_orbit to the section {x=0}
  projection(NJ, GetJetVarDegree(), jperiod, orbit, jet_orbit, jet_px[0][0]);

  // set the output
  AssignJetToJet(jx[0], jet_orbit[1]);
  AssignJetToJet(jx[1], jet_orbit[3]);

  ClearMyFloat(te);
  return 0;
}
