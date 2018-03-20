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

/** @file csr_matrix_test.cc Tests the CSR matrix functionality in
    utilities_basic/csr_matrix by e.g. creating a matrix and checking
    that its values can be fetched correctly. */

#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <cstring>

#include "realtype.h"
#include "template_blas_common.h"
#include "matInclude.h"
#include "csr_matrix.h"
#include "output.h"

static void get_identical_rand_numbers(int n, int nNumbers, std::vector<int> & vec) {
  vec.reserve(nNumbers);
  for(int i = 0; i < nNumbers; i++) {
    int attemptCount = 0;
    int currNumber = -1;
    while(1) {
      currNumber = rand() % n;
      // Check that tmp is not identical to any of the previous numbers
      bool found = false;
      for(int j = 0; j < i; j++) {
	if(currNumber == vec[j])
	  found = true;
      }
      if(found == false)
	break;
      attemptCount++;
      if(attemptCount > 888)
	throw std::runtime_error("Error: too many attempts in get_identical_rand_numbers().");
    } // end while
    vec.push_back(currNumber);
  }
}

int main(int argc, char *argv[])
{
  int n = 222;
  int nnz_per_row = 7;
  if(argc == 3) {
    n = atoi(argv[1]);
    nnz_per_row = atoi(argv[2]);
  }
  printf("csr_matrix_test: n = %d, nnz_per_row = %d\n", n, nnz_per_row);

  //  enable_output(); // Do this if you want the ergoscf.out file, to see timings etc.

  long nnz = (long)n * (long)nnz_per_row;
  std::vector<int> rowind(nnz);
  std::vector<int> colind(nnz);
  long count = 0;
  for(int row = 0; row < n; row++) {
    std::vector<int> tmpVec;
    get_identical_rand_numbers(n, nnz_per_row, tmpVec);
    assert(tmpVec.size() == (size_t)nnz_per_row);
    for(int i = 0; i < nnz_per_row; i++) {
      rowind[count] = row;
      colind[count] = tmpVec[i];
      count++;
    }
  }
  assert(count == nnz);

  csr_matrix_struct A_CSR;
  memset(&A_CSR, 0, sizeof(csr_matrix_struct));
  int symmetryFlagForCSR = 0;

  if(ergo_CSR_create(&A_CSR, symmetryFlagForCSR, n, nnz, &rowind[0], &colind[0]) != 0) {
    printf("Error in ergo_CSR_create.\n");
    return EXIT_FAILURE;
  }

  for(long i = 0; i < nnz; i++) {
    int row = rowind[i];
    int col = colind[i];
    ergo_real element = ergo_CSR_get_element(&A_CSR, row, col);
    if(element != 0)
      throw std::runtime_error("Error: (element != 0).");
    ergo_real x = 0.5;
    if(ergo_CSR_add_to_element(&A_CSR, row, col, x) != 0)
      throw std::runtime_error("Error in ergo_CSR_add_to_element.");
    ergo_real y = ergo_CSR_get_element(&A_CSR, row, col);
    if(y != x)
      throw std::runtime_error("Error: (y != x).");
    if(ergo_CSR_add_to_element(&A_CSR, row, col, x) != 0)
      throw std::runtime_error("Error in ergo_CSR_add_to_element.");
    y = ergo_CSR_get_element(&A_CSR, row, col);
    if(y != (x*2))
      throw std::runtime_error("Error: (y != (x*2)).");
  }

  printf("csr_matrix_test finished OK.\n");

  return EXIT_SUCCESS;
}
