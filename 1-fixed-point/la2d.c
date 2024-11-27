#include "la2d.h"

/**
 * @brief eigen2d
 *    Computes the eigenpairs of a 2-by-2 matrix
 *    (possibly complex)
 *
 * @param[in] n: It must be =2
 * @param[in] a: Matrix of n-by-n
 * @param[in] lda: Leading of a
 * @param[out] vaps: Eigenvalues
 * @param[out] veps: Matrix of unit eigenvectors
 * @param[in] ldv: Leading of veps
 * @return 0 if success
 *
 * @note a and veps can be the same
 */
int eigen2d(int n, MY_FLOAT a[], int lda, MY_FLOAT vaps[], MY_FLOAT veps[], int ldv)
{
  if (n != 2) {
      fprintf(stderr, "%s:%d dimension %d more than 2\n", __FILE__, __LINE__, n); fflush(stderr);
      exit(20);
    }
  int i, j;
  MY_FLOAT tmp, det, tr;
  InitMyFloat(tr);
  InitMyFloat(det);
  InitMyFloat(tmp);

  /* determinant */
  MultiplyMyFloatA(det, aij(a,lda, 0,0), aij(a,lda, 1,1));
  MultiplyMyFloatA(vaps[1], aij(a,lda, 1,0), aij(a,lda, 0,1));

  /* trace/2 */
  AddMyFloatA(tr, aij(a,lda, 0,0), aij(a,lda, 1,1));
  DivideMyFloatSI(tr, tr, 2);

  /* eigenvalues */
  MultiplyMyFloatA(vaps[0], tr, tr);
  SubtractMyFloatA(vaps[0], vaps[0], det);
  AddMyFloatA(vaps[0], vaps[0], vaps[1]);
  sqrtMyFloatA(vaps[1], vaps[0]); // possible complex
  AddMyFloatA(vaps[0], tr, vaps[1]);
  SubtractMyFloatA(vaps[1], tr, vaps[1]);

  /* eigenvectors */
  SubtractMyFloatA(aij(veps,ldv, 1,0), vaps[0], aij(a,lda, 0,0));
  SubtractMyFloatA(aij(veps,ldv, 1,1), vaps[1], aij(a,lda, 0,0));
  AssignMyFloat(aij(veps,ldv, 0,0), aij(a,lda, 0,1));
  AssignMyFloat(aij(veps,ldv, 0,1), aij(a,lda, 0,1));

  /* unitary eigenvectors in norm2 */
  for (j = 0; j < n; ++j)
    {
      MakeMyFloatC(det, "0.0", 0.0);
      for (i = 0; i < n; ++i)
        {
          fabsMyFloatA(tmp, aij(veps,ldv, i,j));
          MultiplyMyFloatA(tmp, tmp, tmp);
          AddMyFloatA(det, det, tmp);
        }
      sqrtMyFloatA(det, det); // norm of a vep
      for (i = 0; i < n; ++i)
        {DivideMyFloatA(aij(veps,ldv, i,j), aij(veps,ldv, i,j), tmp);}
    }

  ClearMyFloat(tmp);
  ClearMyFloat(tr);
  ClearMyFloat(det);
  return 0;
}



/**
 * @brief dgels2d
 *  Solve a 2-by-2 matrix of the form ax=b with
 *  possibly multiple columns in b
 *
 *  It assumes col-major order (it uses the
 *  internal macro aij)
 *
 *  It uses Cramer formula
 *
 * @param[in] a
 * @param[in] lda  leading coefficient (>= 2)
 * @param[inout] b
 * @param[in] ldb  leading coefficient (>= 2)
 * @param[out info 0 means successful
 */
void solve2d(int n, int ns, MY_FLOAT a[], int lda, MY_FLOAT b[], int ldb, double log10tol, int info[1])
{
  if (n != 2) {
      fprintf(stderr, "%s:%d dimension %d more than 2\n", __FILE__, __LINE__, n); fflush(stderr);
      exit(20);
    }
  int j;
  double dtmp;
  MY_FLOAT tmp, tmp2, det;

  InitMyFloat(tmp);
  InitMyFloat(tmp2);
  InitMyFloat(det);

  MultiplyMyFloatA(det, aij(a,lda, 0,0),aij(a,lda, 1,1));
  MultiplyMyFloatA(tmp, aij(a,lda, 0,1),aij(a,lda, 1,0));
  SubtractMyFloatA(det, det,tmp);
  fabsMyFloatA(tmp, det);
  MyFloatToDouble(dtmp, tmp);
  if (log10(dtmp) < log10tol) {*info=1; return;}

  for (j = 0; j < ns; ++j) {
      AssignMyFloat(tmp, aij(b,ldb, 0,j));

      MultiplyMyFloatA(aij(b,ldb, 0,j), aij(b,ldb, 0,j),aij(a,lda, 1,1));
      MultiplyMyFloatA(tmp2, aij(b,ldb, 1,j),aij(a,lda, 0,1));
      SubtractMyFloatA(aij(b,ldb, 0,j), aij(b,ldb, 0,j),tmp2);
      DivideMyFloatA(aij(b,ldb, 0,j), aij(b,ldb, 0,j),det);

      MultiplyMyFloatA(aij(b,ldb, 1,j), aij(b,ldb, 1,j),aij(a,lda, 0,0));
      MultiplyMyFloatA(tmp2, tmp,aij(a,lda, 1,0));
      SubtractMyFloatA(aij(b,ldb, 1,j), aij(b,ldb, 1,j),tmp2);
      DivideMyFloatA(aij(b,ldb, 1,j), aij(b,ldb, 1,j),det);
    }

  *info=0;
  ClearMyFloat(det);
  ClearMyFloat(tmp2);
  ClearMyFloat(tmp);
}
