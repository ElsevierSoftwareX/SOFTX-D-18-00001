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

/** @file mmul_simple_many_test.cc Tests and measures timings for
    matrix-matrix multiplication using BLAS and compares to a naive
    implementation. The idea is to run this linking to different BLAS
    variants with and without threading inside the BLAS gemm routine,
    to see how much speedup can be achieved from threading. */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "realtype.h"
#include "utilities.h"
#include "mat_gblas.h"
#include "config.h" // Needed to get the PRECISION_SINGLE macro

static void fill_matrix_with_random_numbers(int n, std::vector<ergo_real> & A) {
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      ergo_real randomNumber = rand();
      A[i*n+j] = randomNumber / RAND_MAX;
    }
}

static void do_naive_mmul(std::vector<ergo_real> & C, 
			  const std::vector<ergo_real> & A, 
			  const std::vector<ergo_real> & B, 
			  int n) {
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      ergo_real sum = 0;
      for(int k = 0; k < n; k++)
	sum += A[i*n+k] * B[k*n+j];
      C[j*n+i] = sum;
    }
}

ergo_real compare_matrices(const std::vector<ergo_real> & A, 
			   const std::vector<ergo_real> & B, 
			   int n) {
  ergo_real maxabsdiff = 0;
  for(int i = 0; i < n*n; i++) {
    ergo_real absdiff = template_blas_fabs(A[i] - B[i]);
    if(absdiff > maxabsdiff)
      maxabsdiff = absdiff;
  }
  return maxabsdiff;
}

int main(int argc, char *argv[])
{
  int m = 10;
  int n = 500;
  if(argc == 3) {
    m = atof(argv[1]);
    n = atof(argv[2]);
  }
  printf("mmul_simple_many_test start, m = %5d, matrix size n = %6d.\n", m, n);
  // Generate matrix lists A and B of matrices filled with random numbers.
  std::vector< std::vector<ergo_real> > A(m);
  std::vector< std::vector<ergo_real> > B(m);
  for(int i = 0; i < m; i++) {
    A[i].resize(n*n);
    B[i].resize(n*n);
    fill_matrix_with_random_numbers(n, A[i]);
    fill_matrix_with_random_numbers(n, B[i]);
  }
  printf("random matrix lists A and B generated OK.\n");

  // Compute matrix C = A*B using naive implementation.
  std::vector< std::vector<ergo_real> > C(m);
  for(int i = 0; i < m; i++)
    C[i].resize(n*n);
  Util::TimeMeter tm_naive;
  for(int i = 0; i < m; i++)
    do_naive_mmul(C[i], A[i], B[i], n);
  double secondsTaken_naive_mmul = tm_naive.get_wall_seconds() - tm_naive.get_start_time_wall_seconds();
  printf("do_naive_mmul took   %6.3f wall seconds.\n", secondsTaken_naive_mmul);

  // Now do the same computation by calling the BLAS gemm routine.
  ergo_real alpha = 1;
  ergo_real beta = 0;
  std::vector< std::vector<ergo_real> > C2(m);
  for(int i = 0; i < m; i++)
    C2[i].resize(n*n);
  Util::TimeMeter tm_BLAS_gemm_1;
  for(int i = 0; i < m; i++)
    mat::gemm("T", "T", &n, &n, &n, &alpha, 
	      &A[i][0], &n, &B[i][0], &n,
	      &beta, &C2[i][0], &n);
  double secondsTaken_BLAS_gemm_1 = tm_BLAS_gemm_1.get_wall_seconds() - tm_BLAS_gemm_1.get_start_time_wall_seconds();
  printf("BLAS gemm call took   %6.3f wall seconds.\n", secondsTaken_BLAS_gemm_1);
  // Check that results are equal.
  ergo_real maxdiff1 = 0;
  for(int i = 0; i < m; i++) {
    ergo_real diff = compare_matrices(C[i], C2[i], n);
    if(diff > maxdiff1)
      maxdiff1 = diff;
  }
  printf("Max abs diff (elementwise) between naive and BLAS gemm results: %6.3g\n", (double)maxdiff1);

  // Now do the same computation again by calling the BLAS gemm routine.
  std::vector< std::vector<ergo_real> > C3(m);
  for(int i = 0; i < m; i++)
    C3[i].resize(n*n);
  Util::TimeMeter tm_BLAS_gemm_2;
  for(int i = 0; i < m; i++)
    mat::gemm("T", "T", &n, &n, &n, &alpha, 
	      &A[i][0], &n, &B[i][0], &n,
	      &beta, &C3[i][0], &n);
  double secondsTaken_BLAS_gemm_2 = tm_BLAS_gemm_2.get_wall_seconds() - tm_BLAS_gemm_2.get_start_time_wall_seconds();
  printf("BLAS gemm call took   %6.3f wall seconds.\n", secondsTaken_BLAS_gemm_2);
  // Check that results are equal.
  ergo_real maxdiff2 = 0;
  for(int i = 0; i < m; i++) {
    ergo_real diff = compare_matrices(C[i], C3[i], n);
    if(diff > maxdiff1)
      maxdiff2 = diff;
  }
  printf("Max abs diff (elementwise) between naive and BLAS gemm results: %6.3g\n", (double)maxdiff2);

  ergo_real tol = 1e-6;
#ifdef PRECISION_SINGLE
  tol = 5e-4;
#endif
  if(maxdiff1 > tol || maxdiff2 > tol) {
    printf("Error: too large diff between naive mmul and BLAS gemm results.\n");
    return -1;
  }
  printf("mmul_simple_many_test finished OK.\n");
  return 0;
}
