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

/** @file matrix_norm.cc

    @brief Code for computing Euclidean norm of a dense matrix.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <cmath>

#include "matrix_norm.h"
#include "memorymanag.h"
#include "output.h"

#include "../matrix/mat_gblas.h"

static ergo_real
get_largest_eigenvalue(int n, const ergo_real* M)
{
  int lwork = 3*n*n;
  ergo_real* work = (ergo_real*)ergo_malloc(lwork*sizeof(ergo_real));
  ergo_real* eigvalList = (ergo_real*)ergo_malloc(n*sizeof(ergo_real));
  ergo_real* A = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
  memcpy(A, M, n*n*sizeof(ergo_real));

  int info = -1;
  mat::syev("N", "U", &n, A,
	    &n, eigvalList, work, &lwork, 
	    &info);
  if(info != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_largest_eigenvalue, in syev, info = %i", info);
      exit(EXIT_FAILURE);
    }
  
  ergo_real largestEigenvalue = eigvalList[n-1];

  ergo_free(work);
  ergo_free(eigvalList);
  ergo_free(A);

  return largestEigenvalue;
}


ergo_real 
get_euclidean_norm(int m, int n, const ergo_real* A)
{
  if(n > m)
    {
      // Create transpose
      ergo_real* AT = (ergo_real*)ergo_malloc(n*m*sizeof(ergo_real));
      int i, j;
      for(i = 0; i < n; i++)
	for(j = 0; j < m; j++)
	  AT[j*n+i] = A[i*m+j];
      // Compute norm of AT, which is the same as norm of A
      ergo_real normOfTranspose = get_euclidean_norm(n, m, AT);
      ergo_free(AT);
      return normOfTranspose;
    }

  // Create matrix AT*A  ( n * n matrix )
  ergo_real* ATA = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
  int i, j, k;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
	ergo_real sum = 0;
	for(k = 0; k < m; k++)
	  sum += A[i*m+k]*A[j*m+k];
	ATA[i*n+j] = sum;
      }

  // Get largest abs eigenvalue of ATA
  ergo_real largestEigenvalue = get_largest_eigenvalue(n, ATA);

  // The Euclidean norm is given by 
  // the square root of the largest abs eigenvalue of ATA
  ergo_real euclideanNorm = template_blas_sqrt(largestEigenvalue);

  ergo_free(ATA);

  return euclideanNorm;
}
