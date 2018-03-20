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

/** @file mmul_rectangular_test.cc Tests and measures timings for
    matrix-matrix multiplication using BLAS and compares to a naive
    implementation, for rectangular matrices. The idea is to run this
    linking to different BLAS variants with and without threading
    inside the BLAS gemm routine, to see how much speedup can be
    achieved from threading. */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "realtype.h"
#include "utilities.h"
#include "mat_gblas.h"


static void fill_matrix_with_random_numbers(int n, int m, std::vector<ergo_real> & A) {
  for(int i = 0; i < n; i++)
    for(int j = 0; j < m; j++) {
      ergo_real randomNumber = rand();
      A[i*m+j] = randomNumber / RAND_MAX;
    }
}

/* This routine computes C=A*B where A has dim n*m, B has dim m*p, and
   C has dim p*n. */
static void do_naive_mmul(std::vector<ergo_real> & C, 
			  const std::vector<ergo_real> & A, 
			  const std::vector<ergo_real> & B, 
			  int n, int m, int p) {
  for(int i = 0; i < n; i++)
    for(int j = 0; j < p; j++) {
      ergo_real sum = 0;
      for(int k = 0; k < m; k++)
	sum += A[i*m+k] * B[k*p+j];
      C[j*n+i] = sum;
    }
}

/* This routine verifies C=A*B where A has dim n*m, B has dim m*p, and
   C has dim p*n. */
static void verify_mmul_result(const std::vector<ergo_real> & A, 
			       const std::vector<ergo_real> & B, 
			       const std::vector<ergo_real> & C, 
			       int n, int m, int p) {
  double maxabsdiff = 0;
  int count = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < p; j++) {
      // Verify only certain elements to save time.
      if(n > 100 && (rand() % 1000) != 0)
	continue;
      ergo_real sum = 0;
      for(int k = 0; k < m; k++)
	sum += A[i*m+k] * B[k*p+j];
      double absdiff = template_blas_fabs(C[j*n+i] - sum);
      if(absdiff > maxabsdiff)
	maxabsdiff = absdiff;
      count++;
    }
  printf("verify_mmul_result, verified %d elements, maxabsdiff = %g\n", count, maxabsdiff);
}

ergo_real compare_matrices(const std::vector<ergo_real> & A, 
			   const std::vector<ergo_real> & B, 
			   int n, int m) {
  ergo_real maxabsdiff = 0;
  for(int i = 0; i < n*m; i++) {
    ergo_real absdiff = template_blas_fabs(A[i] - B[i]);
    if(absdiff > maxabsdiff)
      maxabsdiff = absdiff;
  }
  return maxabsdiff;
}

#if 0
static void print_matrix(int n, int m, const std::vector<ergo_real> & A, const char* s) {
  printf("Matrix '%s':\n", s);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      printf(" %9.4f", (double)A[j*n+i]);
    }
    printf("\n");
  }
}
#endif

int main(int argc, char *argv[])
{
  int n = 500;
  int m = 500;
  int p = 500;
  if(argc == 4) {
    n = atof(argv[1]);
    m = atof(argv[2]);
    p = atof(argv[3]);
  }
  printf("mmul_simple_test start, matrix sizes n = %6d, m = %6d, p = %6d.\n", n, m, p);
  // Generate matrices A and B filled with random numbers.
  std::vector<ergo_real> A(n*m);
  std::vector<ergo_real> B(m*p);
  fill_matrix_with_random_numbers(n, m, A);
  fill_matrix_with_random_numbers(m, p, B);
  printf("random matrices A and B generated OK.\n");

  // Compute matrix C = A*B using naive implementation.
  std::vector<ergo_real> C(n*p);
  Util::TimeMeter tm_naive;
  do_naive_mmul(C, A, B, n, m, p);
  double secondsTaken_naive_mmul = tm_naive.get_wall_seconds() - tm_naive.get_start_time_wall_seconds();
  printf("do_naive_mmul took   %6.3f wall seconds.\n", secondsTaken_naive_mmul);
  verify_mmul_result(A, B, C, n, m, p);

  // Now do the same computation by calling the BLAS gemm routine.
  ergo_real alpha = 1;
  ergo_real beta = 0;
  std::vector<ergo_real> C2(n*p);
  Util::TimeMeter tm_BLAS_gemm_1;
  int M = n;
  int N = p;
  int K = m;
  mat::gemm("T", "T", &M, &N, &K, &alpha, 
	    &A[0], &K, &B[0], &N,
	    &beta, &C2[0], &M);
  double secondsTaken_BLAS_gemm_1 = tm_BLAS_gemm_1.get_wall_seconds() - tm_BLAS_gemm_1.get_start_time_wall_seconds();
  printf("BLAS gemm call took   %6.3f wall seconds.\n", secondsTaken_BLAS_gemm_1);
  // Check that results are equal.
  ergo_real diff1 = compare_matrices(C, C2, n, p);
  printf("Max abs diff (elementwise) between naive and BLAS gemm results: %6.3g\n", (double)diff1);
  verify_mmul_result(A, B, C2, n, m, p);

  // Now do the same computation by again calling the BLAS gemm routine.
  std::vector<ergo_real> C3(n*p);
  Util::TimeMeter tm_BLAS_gemm_2;
  mat::gemm("T", "T", &M, &N, &K, &alpha, 
	    &A[0], &K, &B[0], &N,
	    &beta, &C3[0], &M);
  double secondsTaken_BLAS_gemm_2 = tm_BLAS_gemm_2.get_wall_seconds() - tm_BLAS_gemm_2.get_start_time_wall_seconds();
  printf("BLAS gemm call took   %6.3f wall seconds.\n", secondsTaken_BLAS_gemm_2);
  // Check that results are equal.
  ergo_real diff2 = compare_matrices(C, C3, n, p);
  printf("Max abs diff (elementwise) between naive and BLAS gemm results: %6.3g\n", (double)diff2);
  verify_mmul_result(A, B, C3, n, m, p);

  ergo_real tol = 2e-4;
  if(diff1 > tol || diff2 > tol) {
    printf("Error: too large diff between naive mmul and BLAS gemm results.\n");
    return -1;
  }

  printf("mmul_simple_test finished OK.\n");
  return 0;
}
