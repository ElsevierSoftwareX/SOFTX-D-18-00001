/* Ergo, version 3.6, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2017 Elias Rudberg, Emanuel H. Rubensson, Pawel Salek,
 * and Anastasia Kruchinina.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

/** @file lapack_test.cc Tests some LAPACK functions
    such as template_lapack_???() etc to
    see that they are working properly and that they deliver
    the expected accuracy. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <vector>

#include "realtype.h"
#include "template_lapack_common.h"
#include "matInclude.h"


#if 0
static void print_matrix(int n, const ergo_real* A, const char* name)
{
  printf("printing matrix '%s', n = %3i\n", name, n);
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
	printf(" %8.3f", (double)A[j*n+i]);
      printf("\n");
    }      
}
#endif


static ergo_real get_maxabsdiff(int n, const ergo_real* x, const ergo_real* y)
{
  ergo_real maxabsdiff = 0;
  for(int i = 0; i < n; i++)
    {
      ergo_real diff = x[i] - y[i];
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      if(absdiff > maxabsdiff)
	maxabsdiff = absdiff;
    }
  return maxabsdiff;
}


static int test_gesv(ergo_real tolerance, bool verbose)
{
  if(verbose)
    printf("Testing solution of a linear system of equations using routine template_lapack_gesv()..\n");
  
  // Generate random matrix A and random vector x
  const int n = 77;

  std::vector<ergo_real> A(n*n);
  for(int i = 0; i < n*n; i++)
    A[i] = (ergo_real)rand() / (ergo_real)RAND_MAX;
  
  std::vector<ergo_real> x(n);
  for(int i = 0; i < n; i++)
    x[i] = (ergo_real)rand() / (ergo_real)RAND_MAX;

  // Now get right-hand-side b as A*x = b
  std::vector<ergo_real> b(n);
  for(int i = 0; i < n; i++)
    {
      ergo_real sum = 0;
      for(int j = 0; j < n; j++)
	sum += A[j*n+i] * x[j];
      b[i] = sum;
    }

  // Now use A and b to solve for x.
  
  int NRHS = 1;
  int n2 = n;
  int info = -1;
  std::vector<int> IPIV(n);
  std::vector<ergo_real> Atmp(n*n);
  std::vector<ergo_real> resultVector(n);
  memcpy(&Atmp[0], &A[0], n*n*sizeof(ergo_real));
  memcpy(&resultVector[0], &b[0], n*sizeof(ergo_real));
  
  template_lapack_gesv(&n2, &NRHS, &A[0], &n2, &IPIV[0], &resultVector[0], &n2, &info);
  if(info != 0)
    {
      printf("ERROR in template_lapack_gesv\n");
      return -1;
    }

  // Now compare resultVector with known x.
  ergo_real maxabsdiff = get_maxabsdiff(n, &resultVector[0], &x[0]);

  if(verbose)
    printf("maxabsdiff for template_lapack_gesv: %g\n", (double)maxabsdiff);
  int failed = 0;
  if(maxabsdiff > tolerance)
    {
      printf("template_lapack_gesv test FAILED.\n");
      failed = 1;
    }
  else {
    if(verbose)
      printf("template_lapack_gesv test OK.\n");
  }
  
  return failed;
}


static int test_potf2_trtri(ergo_real tolerance, bool verbose)
{
  int failed = 0;
  
  const int n = 13;

  // Create random upper triangular matrix U.

  std::vector<ergo_real> U(n*n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	ergo_real matrixElement = 0;
	if(i >= j)
	  matrixElement = (ergo_real)rand() / (ergo_real)RAND_MAX;
	if (i == j)
	  matrixElement = matrixElement + 2;
	U[i*n+j] = matrixElement;
      }

  
  // Compute matrix A = U' * U

  std::vector<ergo_real> A(n*n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < n; k++)
	  sum += U[i*n+k] * U[j*n+k];
	A[i*n+j] = sum;
      }


  // Compute Cholesky factorization of A using template_lapack_potf2()
  
  std::vector<ergo_real> Atmp(n*n);
  memcpy(&Atmp[0], &A[0], n*n*sizeof(ergo_real));
  int info = -1;

  template_lapack_potf2("U", &n, &Atmp[0], &n, &info);
  if(info != 0)
    {
      printf("error in template_lapack_potf2\n");
      return -1;
    }

  // Set rest to zero. This is needed because the potf2 operation leaves some 
  // garbage in the other part of the matrix space.
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	if(i < j)
	  Atmp[i*n+j] = 0;
      }

  // Now compare Atmp with U.
  ergo_real maxabsdiff_potf2 = get_maxabsdiff(n*n, &Atmp[0], &U[0]);
  if(verbose)
    printf("maxabsdiff for Cholesky decomposition template_lapack_potf2(): %g\n", (double)maxabsdiff_potf2);
  if(maxabsdiff_potf2 > tolerance)
    {
      printf("ERROR: template_lapack_potf2 not accurate enough.\n");
      printf("Error is %g, tolerance is %g\n", (double)maxabsdiff_potf2, (double)tolerance);
      failed++;
    }
  

  // Compute inverse of U using template_lapack_trtri().
  std::vector<ergo_real> Z(n*n);
  memcpy(&Z[0], &U[0], n*n*sizeof(ergo_real));
  info = -1;
  // use non-const strings uplo and diag needed to avoid compiler warnings. 
  char uplo[8];
  char diag[8];
  uplo[0] = 'U';
  uplo[1] = '\0';
  diag[0] = 'N';
  diag[1] = '\0';
  template_lapack_trtri(uplo, diag, &n, &Z[0], &n, &info);
  if(info != 0)
    {
      printf("error in template_lapack_trtri\n");
      return -1;
    }
  
  // Compute B = Z * U  ( B should be almost an identity matrix. )
  std::vector<ergo_real> B(n*n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < n; k++)
	  sum += Z[k*n+i] * U[j*n+k];
	B[j*n+i] = sum;
      }
  
  // Compute C = U * Z  ( C should be almost an identity matrix. )
  std::vector<ergo_real> C(n*n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < n; k++)
	  sum += U[k*n+i] * Z[j*n+k];
	C[j*n+i] = sum;
      }
  
  // Construct an identity matrix for comparison.
  std::vector<ergo_real> I(n*n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	I[i*n+j] = 0;
	if(i == j)
	  I[i*n+j] = 1;
      }
  
  ergo_real maxabsdiff_trtri_1 = get_maxabsdiff(n*n, &I[0], &B[0]);
  if(verbose)
    printf("maxabsdiff 1 for inverse template_lapack_trtri(): %g\n", (double)maxabsdiff_trtri_1);
  if(maxabsdiff_trtri_1 > tolerance)
    {
      printf("ERROR: template_lapack_trtri not accurate enough.\n");
      printf("Error is %g, tolerance is %g\n", (double)maxabsdiff_trtri_1, (double)tolerance);
      failed++;
    }

  ergo_real maxabsdiff_trtri_2 = get_maxabsdiff(n*n, &I[0], &C[0]);
  if(verbose)
    printf("maxabsdiff 2 for inverse template_lapack_trtri(): %g\n", (double)maxabsdiff_trtri_2);
  if(maxabsdiff_trtri_2 > tolerance)
    {
      printf("ERROR: template_lapack_trtri not accurate enough.\n");
      printf("Error is %g, tolerance is %g\n", (double)maxabsdiff_trtri_2, (double)tolerance);
      failed++;
    }

  if(!failed && verbose)
    printf("Tests of potf2 and trtri finished OK.\n");

  return failed;
}


static int test_stevr_small(ergo_real tolerance, bool verbose)
{
  int failed = 0;

  int size = 2;
  std::vector<ergo_real> alpha(size);
  std::vector<ergo_real> beta(size);
  alpha[0] = 1.419534;
  alpha[1] = 1.117599;
  beta[0] = 0.587297;
  beta[1] = 0.383218;
  int lowInd = 2;
  int uppInd = 2;
  int nEigsWanted = 1;
  int nEigsFound = 0;
  std::vector<ergo_real> eigArray(size);
  std::vector<ergo_real> eigVectors(size*nEigsWanted);
  std::vector<int> isuppz(2*nEigsWanted);
  int lwork = -1;
  int liwork = -1;
  ergo_real dummy = -1.0;
  ergo_real abstol = 0;
  int info;
  // First do a workspace query:
  ergo_real work_query;
  int iwork_query;
  template_lapack_stevr("V", "I", &size, &alpha[0], &beta[0],
			&dummy, &dummy, &lowInd, &uppInd, &abstol,
			&nEigsFound, &eigArray[0], &eigVectors[0], &size, &isuppz[0],
			&work_query, &lwork, &iwork_query, &liwork, &info);
  lwork = int(work_query);
  liwork = iwork_query;
  std::vector<ergo_real> work(lwork);
  std::vector<int> iwork(liwork);
  template_lapack_stevr("V", "I", &size, &alpha[0], &beta[0],
			&dummy, &dummy, &lowInd, &uppInd, &abstol,
			&nEigsFound, &eigArray[0], &eigVectors[0], &size, &isuppz[0],
			&work[0], &lwork, &iwork[0], &liwork, &info);
  if (info)
    std::cout << "info = " << info <<std::endl;
  if(info != 0 || nEigsFound != nEigsWanted) {
    printf("Error: (info != 0 || nEigsFound != nEigsWanted) after stevr call.\n");
    printf("info = %d\n", info);
    printf("nEigsWanted = %d\n", nEigsWanted);
    printf("nEigsFound  = %d\n", nEigsFound);
    failed++;
  }

  if(!failed && verbose)
    printf("Small test of stevr finished OK.\n");

  return failed;
}


int main(int argc, char *argv[])
{
  int failed = 0;
  int verbose = getenv("VERBOSE") != NULL;
  ergo_real machine_epsilon = mat::getMachineEpsilon<ergo_real>();
  
  printf("machine_epsilon = %g Run with env VERBOSE for more info.\n",
         (double)machine_epsilon);
  
  if(test_gesv(machine_epsilon*500, verbose) != 0)
    failed++;

  
  if(test_potf2_trtri(machine_epsilon*10000, verbose) != 0)
    failed++;


  if(test_stevr_small(machine_epsilon*500, verbose) != 0)
    failed++;
  

  if (!failed)
    puts("LAPACK test succeeded.");
  else
    puts("LAPACK test FAILED.");

  return failed;
}
