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

/** @file matrix_algebra.cc

    @brief A few matrix algebra routines for dense matrices.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>

#include "matrix_algebra.h"
#include "memorymanag.h"
#include "output.h"

#include "../matrix/mat_gblas.h"


#define USE_BLAS_MM


#if 0
#ifdef USE_BLAS_MM
#ifdef __cplusplus
extern "C" 
#endif
void dgemm_(const char *ta,const char *tb,
	    const int *n, const int *k, const int *l,
	    const double *alpha,const double *A,const int *lda,
	    const double *B, const int *ldb,
	    const double *beta,const double *C, const int *ldc);
#endif
#endif


/*
  Standard matrix multiplication.
*/
void 
multiply_matrices_general(int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB)
{
  int i, j;
  if(An2 != Bn1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, 
		"error in multiply_matrices_general: (An2 != Bn1)");
      exit(0);
    }
#ifdef USE_BLAS_MM
  if(An1 == 0 || An2 == 0 || Bn1 == 0 || Bn2 == 0)
    return;
  /*  call gemm */
  ergo_real alpha = 1;
  ergo_real beta = 0;
  ergo_real* ABtemp = (ergo_real*)ergo_malloc(An1*Bn2*sizeof(ergo_real));
  memset(ABtemp, 0, An1*Bn2*sizeof(ergo_real));
  mat::gemm("T", "T", &An1, &Bn2, &An2, &alpha, 
	    A, &An2, 
	    B, &Bn2, 
	    &beta, 
	    ABtemp, &An1);
  /*  do transpose of result */
  for(i = 0; i < An1; i++)
    for(j = 0; j < Bn2; j++)
      {
	AB[i*Bn2+j] = ABtemp[j*An1+i];
      }  
  ergo_free(ABtemp);
#else
  for(i = 0; i < An1; i++)
    for(j = 0; j < Bn2; j++)
      {
	ergo_real sum = 0;
	for(int k = 0; k < An2; k++)
	  sum += A[i*An2+k] * B[k*Bn2+j];
	AB[i*Bn2+j] = sum;
      }
#endif
}


/*
  Matrix multiplication when the first matrix is transposed.
*/
void 
multiply_matrices_general_T_1(int An1, int An2, int Bn1, int Bn2, const ergo_real* A, const ergo_real* B, ergo_real* AB)
{
  int i, j;
  if(An1 != Bn1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, "error in multiply_matrices_general_T_1: (An1 != Bn1)");
      exit(0);
    }
  ergo_real* At = (ergo_real*)ergo_malloc(An1*An2*sizeof(ergo_real));
  for(i = 0; i < An1; i++)
    for(j = 0; j < An2; j++)
      At[j*An1+i] = A[i*An2+j];
  multiply_matrices_general(An2, An1, Bn1, Bn2, At, B, AB);
  ergo_free(At);
}
