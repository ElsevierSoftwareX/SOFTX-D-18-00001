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

/** @file simple_sparse_mat.h

    @brief Simple sparse matrix implementation.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef SIMPLE_SPARSE_MAT_HEADER
#define SIMPLE_SPARSE_MAT_HEADER

#include "realtype.h"

typedef struct {
  int i;
  int j;
  int same_i_count;
  ergo_real value;
} i_j_val_struct;

int spmat_multiply_matrices(const i_j_val_struct* A, int nnzA, 
			    const i_j_val_struct* B, int nnzB, 
			    i_j_val_struct* C, int M, int N);

int spmat_sort_elements(i_j_val_struct* A, int nnzA);

#endif
