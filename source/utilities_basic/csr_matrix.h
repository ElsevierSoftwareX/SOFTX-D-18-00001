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

/** @file csr_matrix.h

    @brief Compressed sparse row (CSR) matrix functionality.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef CSR_MATRIX_HEADER
#define CSR_MATRIX_HEADER

#include "realtype.h"


typedef struct
{
  long noOfElementsInRow;
  long firstElementIndex;
} csr_matrix_row_struct;

struct csr_matrix_struct
{
  int n;
  long nnz;
  int symmetryFlag;
  csr_matrix_row_struct* rowList;
  ergo_real* elementList;
  int* columnIndexList;
};



int ergo_CSR_create(csr_matrix_struct* csr, 
		    int symmetryFlag,
		    int n,
		    long nnz,
		    int* rowind,
		    int* colind);

int ergo_CSR_destroy(csr_matrix_struct* csr);

int ergo_CSR_copy(csr_matrix_struct* csrDest, const csr_matrix_struct* csrSource);

int ergo_CSR_add_equal_structure(csr_matrix_struct* csrDest, const csr_matrix_struct* csrSource);

int ergo_CSR_add_to_element(csr_matrix_struct* csr, 
			    int row,
			    int col,
			    ergo_real value);

ergo_real ergo_CSR_get_element(const csr_matrix_struct* csr, 
			       int row,
			       int col);

ergo_real ergo_CSR_get_max_abs_element(const csr_matrix_struct* csr);

long ergo_CSR_get_nvalues(const csr_matrix_struct* csr);

int ergo_CSR_get_values(const csr_matrix_struct* csr,
			int* rowind,
			int* colind,
			ergo_real* values,
			long nvalues);

int ergo_CSR_get_nvalues_singlerow(const csr_matrix_struct* csr,
				   int row);

int ergo_CSR_get_values_singlerow(const csr_matrix_struct* csr,
				  int row,
				  int* colind,
				  ergo_real* values,
				  int nvalues);



#endif
