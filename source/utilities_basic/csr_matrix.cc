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

/** @file csr_matrix.cc

    @brief Compressed sparse row (CSR) matrix functionality.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>

#include "csr_matrix.h"
#include "output.h"
#include "memorymanag.h"
#include "mat_gblas.h"


typedef struct
{
  int row;
  int col;
} csr_index_pair_struct;



static int
csr_compare_index_pairs_for_qsort(const void* p1, const void* p2)
{
  csr_index_pair_struct* pair_1 = (csr_index_pair_struct*)p1;
  csr_index_pair_struct* pair_2 = (csr_index_pair_struct*)p2;
  if(pair_1->row > pair_2->row)
    return 1;
  if(pair_1->row < pair_2->row)
    return -1;
  if(pair_1->col > pair_2->col)
    return 1;
  if(pair_1->col < pair_2->col)
    return -1;
  return 0;
}


int 
ergo_CSR_create(csr_matrix_struct* csr, 
		int symmetryFlag,
		int n,
		long nnz,
		int* rowind,
		int* colind)
{
  csr_index_pair_struct* list;

  csr->symmetryFlag = symmetryFlag;
  csr->n = n;
  csr->nnz = nnz;
  csr->rowList = new csr_matrix_row_struct[n];
  memset(csr->rowList, 0, n*sizeof(csr_matrix_row_struct));
  csr->columnIndexList = new int[nnz];
  memset(csr->columnIndexList, 0, nnz*sizeof(int));
  /* Note: we do not need to allocate elementList yet, better to wait
     until after the temporary index pairs list has been deleted.  */
  
  /* Create list of index pairs, and sort it. */
  list = new csr_index_pair_struct[nnz];

  for(long i = 0; i < nnz; i++)
    {
      int row = rowind[i];
      int col = colind[i];
      if(symmetryFlag == 1)
	{
	  if(row > col)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, "error in CSR_create: symmetry not satisfied.");
	      return -1;
	    }
	}
      list[i].row = row;
      list[i].col = col;
    }
  qsort(list, nnz, sizeof(csr_index_pair_struct), csr_compare_index_pairs_for_qsort);
  
  /* Now the list of index pairs is sorted by rows and within each row by columns. */
  /* Create CSR structure one row at a time. */
  long currIndex = 0;
  for(int row = 0; row < n; row++)
    {
      csr->rowList[row].firstElementIndex = currIndex;
      if(currIndex == nnz)
	{
	  /* no elements left; this row must be empty. */
	  csr->rowList[row].noOfElementsInRow = 0;
	  continue;
	}
      int count = 0;
      while(list[currIndex].row == row)
	{
	  csr->columnIndexList[currIndex] = list[currIndex].col;
	  count++;
	  currIndex++;
	  if(currIndex == nnz)
	    break;
	}
      csr->rowList[row].noOfElementsInRow = count;
      if(currIndex < nnz)
	{
	  if(list[currIndex].row < row)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, "error in CSR_create: list not sorted by rows.");
	      return -1;
	    }
	}
    } /* END FOR row */

  delete [] list;
  list = NULL;

  // Now allocate elementList. Better to do it now, after we have deleted the temporary list.
  csr->elementList = new ergo_real[nnz];
  memset(csr->elementList, 0, nnz*sizeof(ergo_real));

  /* Check that column indexes are properly sorted within each row. */
  for(int row = 0; row < n; row++)
    {
      for(int col = 1; col < csr->rowList[row].noOfElementsInRow; col++)
	{
	  long ii = csr->rowList[row].firstElementIndex;
	  if(csr->columnIndexList[ii+col] <= csr->columnIndexList[ii+col-1])
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_LOWLEVEL, "error in CSR_create: list not sorted within row.");
	      return -1;	      
	    }
	}
    }

  return 0;
}


int ergo_CSR_destroy(csr_matrix_struct* csr)
{
  delete [] csr->rowList;
  delete [] csr->elementList;
  delete [] csr->columnIndexList;
  memset(csr, 0, sizeof(csr_matrix_struct));
  return 0;
}


int ergo_CSR_copy(csr_matrix_struct* csrDest, const csr_matrix_struct* csrSource)
{
  int n = csrSource->n;
  long nnz = csrSource->nnz;

  csrDest->n = n;
  csrDest->nnz = nnz;
  csrDest->symmetryFlag = csrSource->symmetryFlag;

  csrDest->rowList = new csr_matrix_row_struct[n];
  csrDest->elementList = new ergo_real[nnz];
  csrDest->columnIndexList = new int[nnz];
  
  memcpy(csrDest->rowList, csrSource->rowList, n*sizeof(csr_matrix_row_struct));
  memcpy(csrDest->elementList, csrSource->elementList, nnz*sizeof(ergo_real));
  memcpy(csrDest->columnIndexList, csrSource->columnIndexList, nnz*sizeof(int));

  return 0;
}


int ergo_CSR_add_equal_structure(csr_matrix_struct* csrDest, const csr_matrix_struct* csrSource)
{
  int n = csrSource->n;
  long nnz = csrSource->nnz;
  
  // Check that matrices have identical structure
  if(csrDest->symmetryFlag != csrSource->symmetryFlag)
    return -1;
  if(memcmp(csrDest->rowList, csrSource->rowList, n*sizeof(csr_matrix_row_struct)) != 0)
    return -1;
  if(memcmp(csrDest->columnIndexList, csrSource->columnIndexList, nnz*sizeof(int)) != 0)
    return -1;

  // OK, add elements of source to elements of dest, store in dest.
  for(long i = 0; i < nnz; i++)
    csrDest->elementList[i] += csrSource->elementList[i];

  return 0;
}


static long
ergo_csr_find_index(const csr_matrix_struct* csr, int row, int col)
{
  if(row < 0 || row >= csr->n)
    throw std::runtime_error("Error: ergo_csr_find_index called with (row < 0 || row >= csr->n).");
  int n = csr->rowList[row].noOfElementsInRow;
  if(n <= 0) /* If row is empty we can return already here. */
    return -1;
  long baseIndex = csr->rowList[row].firstElementIndex;
  int* colList = &csr->columnIndexList[baseIndex];
  int lo = 0;
  int hi = n-1;
  while(lo < hi - 1) {
    int mid = (lo + hi) / 2;
    if(colList[mid] < col)
      lo = mid;
    else
      hi = mid;
  }
  if(colList[lo] == col)
    return baseIndex + lo;
  if(colList[hi] == col)
    return baseIndex + hi;
  /* Not found */
  return -1;
}


int 
ergo_CSR_add_to_element(csr_matrix_struct* csr, 
			int row,
			int col,
			ergo_real value)
{
  if(csr == NULL)
    throw std::runtime_error("Error: ergo_CSR_add_to_element called for (csr == NULL).");
  int row2 = row;
  int col2 = col;
  if(csr->symmetryFlag)
    {
      if(row > col)
	{
	  row2 = col;
	  col2 = row;
	}
    }
  long i = ergo_csr_find_index(csr, row2, col2);
  if(i < 0)
    return 0;
  csr->elementList[i] += value;
  return 0;
}


ergo_real 
ergo_CSR_get_element(const csr_matrix_struct* csr, 
		     int row,
		     int col)
{
  int row2 = row;
  int col2 = col;
  if(csr->symmetryFlag)
    {
      if(row > col)
	{
	  row2 = col;
	  col2 = row;
	}
    }
  long i = ergo_csr_find_index(csr, row2, col2);
  if(i < 0)
    return 0;
  return csr->elementList[i];
}


ergo_real 
ergo_CSR_get_max_abs_element(const csr_matrix_struct* csr)
{
  ergo_real maxabs = 0;
  for(long i = 0; i < csr->nnz; i++)
    {
      ergo_real absval = template_blas_fabs(csr->elementList[i]);
      if(absval > maxabs)
	maxabs = absval;
    }
  return maxabs;
}


long
ergo_CSR_get_nvalues(const csr_matrix_struct* csr)
{
  return csr->nnz;
}


int 
ergo_CSR_get_values(const csr_matrix_struct* csr,
		    int* rowind,
		    int* colind,
		    ergo_real* values,
		    long nvalues)
{
  int count = 0;
  
  if(nvalues != csr->nnz)
    return -1;
  for(int i = 0; i < csr->n; i++)
    {
      long baseIndex = csr->rowList[i].firstElementIndex;
      int* colList = &csr->columnIndexList[baseIndex];
      ergo_real* valueList = &csr->elementList[baseIndex];
      for(int j = 0; j < csr->rowList[i].noOfElementsInRow; j++)
	{
	  rowind[count] = i;
	  colind[count] = colList[j];
	  values[count] = valueList[j];
	  count++;
	}
    }
  if(count != nvalues)
    return -1;
  return 0;
}


int 
ergo_CSR_get_nvalues_singlerow(const csr_matrix_struct* csr,
			       int row)
{
  return csr->rowList[row].noOfElementsInRow;
}

int 
ergo_CSR_get_values_singlerow(const csr_matrix_struct* csr,
			      int row,
			      int* colind,
			      ergo_real* values,
			      int nvalues)
{
  if(nvalues != csr->rowList[row].noOfElementsInRow)
    return -1;
  long baseIndex = csr->rowList[row].firstElementIndex;
  int* colList = &csr->columnIndexList[baseIndex];
  ergo_real* valueList = &csr->elementList[baseIndex];
  for(int i = 0; i < nvalues; i++) {
    colind[i] = colList[i];
    values[i] = valueList[i];
  }
  return 0;
}

