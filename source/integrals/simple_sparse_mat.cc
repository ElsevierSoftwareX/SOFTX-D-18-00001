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

/** @file simple_sparse_mat.cc

    @brief Simple sparse matrix implementation.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "simple_sparse_mat.h"
#include <cmath>
#include <stdexcept>
#include <cstdlib>

static int
compare_i_j_val_structs(const void* p1, const void* p2)
{
  i_j_val_struct* struct_1 = (i_j_val_struct*)p1;
  i_j_val_struct* struct_2 = (i_j_val_struct*)p2;
  const int returnValue1 = -1;
  const int returnValue2 =  1;
  if(struct_1->i < struct_2->i)
    return returnValue1;
  if(struct_1->i > struct_2->i)
    return returnValue2;
  if(struct_1->j < struct_2->j)
    return returnValue1;
  if(struct_1->j > struct_2->j)
    return returnValue2;
  return 0;
}

int spmat_sort_elements(i_j_val_struct* A, int nnzA) {
  qsort(A, nnzA, sizeof(i_j_val_struct), compare_i_j_val_structs);
  // Verify sort
  for(int k = 0; k < nnzA-1; k++) {
    if(A[k+1].i < A[k].i)
      throw std::runtime_error("ERROR in spmat_sort_elements: list not properly sorted.");
  }
  for(int k = 0; k < nnzA; k++) {
    // Check how many elements have the same i-value
    int kk = k;
    while(kk < nnzA) {
      if(A[kk].i != A[k].i)
	break;
      kk++;
    }
    // Now kk is the index of the first element found that does not have the same i-value (or kk==nnzA)
    A[k].same_i_count = kk - k;
  }
  return 0;
}

int spmat_multiply_matrices(const i_j_val_struct* A, int nnzA, const i_j_val_struct* B, int nnzB, i_j_val_struct* C, int M, int N) {
  i_j_val_struct* Cbuf = new i_j_val_struct[M*N];
  int counters[M];
  for(int i = 0; i < M; i++)
    counters[i] = 0;
  for(int idxA = 0; idxA < nnzA; idxA++)
    for(int idxB = 0; idxB < nnzB; idxB++) {
      if(A[idxA].j != B[idxB].i)
	continue;
      int i = A[idxA].i;
      int j = B[idxB].j;
      // OK, we have a contribution to C_i_j
      ergo_real contribValue = A[idxA].value * B[idxB].value;
      if(i < 0 || i >= M || j < 0 || j >= N)
	throw std::runtime_error("ERROR in multiply_matrices_sp: i, j out of bounds.");
      // Check if there is already an entry for this index pair.
      bool found = false;
      for(int k = 0; k < counters[i]; k++) {
	i_j_val_struct & p = Cbuf[i*N+k];
	if(p.i == i && p.j == j) {
	  p.value += contribValue;
	  found = true;
	  break;
	}
      }
      if(found == false) {
	if(counters[i] >= N)
	  throw std::runtime_error("ERROR in multiply_matrices_sp: (counters[i] >= N)");
	i_j_val_struct & p = Cbuf[i*N+counters[i]];
	counters[i]++;
	p.i = i;
	p.j = j;
	p.same_i_count = 1;
	p.value = contribValue;
      }
    }
  // Now create final result in C
  int count = 0;
  for(int i = 0; i < M; i++)
    for(int k = 0; k < counters[i]; k++) {
      C[count] = Cbuf[i*N+k];
      count++;
    }
  delete [] Cbuf;
  return count;
}

