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

/** @file hermite_conversion_prep.cc

    @brief Code for conversion between integrals computed for Hermite
    Gaussians and Cartesian Gaussians, using saved information from a
    preparatory step that can be done once and re-used many times.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "hermite_conversion_prep.h"
#include <cmath>
#include <memory.h>
#include <assert.h>
#include <vector>
#include "hermite_conversion_symb.h"
#include "monomial_info.h"
#include "utilities.h"
#include "mat_gblas.h"

#if BASIS_FUNC_POLY_MAX_DEGREE<6
const int MAX_NO_OF_CONTRIBS = 1000000;
#else
const int MAX_NO_OF_CONTRIBS = 10000000;
#endif


void hermite_conversion_info_struct::init(const monomial_info_struct & monomial_info) {
  Util::TimeMeter timeMeter;
  const int nmax = HERMITE_CONVERSION_MAX_N;

  // Allocate work memory buffers.
  std::vector<hermite_conversion_contrib_struct> currlist(MAX_NO_OF_CONTRIBS);
  std::vector<hermite_conversion_element_struct> currlist_elements(MAX_NO_OF_CONTRIBS);
  const int nMonMax = monomial_info.no_of_monomials_list[nmax];
  std::vector<symb_matrix_element> list(nMonMax*nMonMax);
  const int invFlag = 1;

  // Get list_right
  for(int n1 = 0; n1 <= nmax; n1++) {
    get_hermite_conversion_matrix_symb(&monomial_info, n1, invFlag, &list[0]);
    for(int n2 = 0; n2 <= nmax; n2++) {
      int count = 0;
      int nMon1 = monomial_info.no_of_monomials_list[n1];
      int nMon2 = monomial_info.no_of_monomials_list[n2];
      for(int j = 0; j < nMon1; j++) {
	for(int k = 0; k < nMon1; k++) {
	  int idx = j*nMon1+k;
	  if(template_blas_fabs(list[idx].coeff) > 1e-5) {
	    for(int i = 0; i < nMon2; i++) {
	      assert(count < MAX_NO_OF_CONTRIBS);
	      currlist[count].destIndex = j*nMon2+i;
	      currlist[count].sourceIndex = k*nMon2+i;
	      currlist[count].a_power = list[idx].ia;
	      currlist[count].coeff = list[idx].coeff;
	      count++;
	    }
	  }
	}
      }
      list_right[n1][n2] = new hermite_conversion_contrib_struct[count];
      memcpy(list_right[n1][n2], &currlist[0], 
	     count*sizeof(hermite_conversion_contrib_struct));
      counters_right[n1][n2] = count;
    }
  }

  // Get list_left
  for(int n2 = 0; n2 <= nmax; n2++) {
    get_hermite_conversion_matrix_symb(&monomial_info, n2, invFlag, &list[0]);
    for(int n1 = 0; n1 <= nmax; n1++) {
      int count = 0;
      int nMon1 = monomial_info.no_of_monomials_list[n1];
      int nMon2 = monomial_info.no_of_monomials_list[n2];
      for(int i = 0; i < nMon2; i++) {
	for(int k = 0; k < nMon2; k++) {
	  int idx = i*nMon2+k;
	  if(template_blas_fabs(list[idx].coeff) > 1e-5) {
	    for(int j = 0; j < nMon1; j++) {
	      assert(count < MAX_NO_OF_CONTRIBS);
	      currlist[count].destIndex = j*nMon2+i;
	      currlist[count].sourceIndex = j*nMon2+k;
	      currlist[count].a_power = list[idx].ia;
	      currlist[count].coeff = list[idx].coeff;
	      count++;
	    }
	  }
	}
      }
      list_left[n1][n2] = new hermite_conversion_contrib_struct[count];
      memcpy(list_left[n1][n2], &currlist[0], 
	     count*sizeof(hermite_conversion_contrib_struct));
      counters_left[n1][n2] = count;
    }
  }

  // Get list_right_simple
  for(int n1 = 0; n1 <= nmax; n1++) {
    get_hermite_conversion_matrix_symb(&monomial_info, n1, invFlag, &list[0]);
    int nMon1 = monomial_info.no_of_monomials_list[n1];
    int count = 0;
    for(int j = 0; j < nMon1; j++) {
      for(int k = 0; k < nMon1; k++) {
	int idx = j*nMon1+k;
	if(template_blas_fabs(list[idx].coeff) > 1e-5) {
	  assert(count < MAX_NO_OF_CONTRIBS);
	  currlist_elements[count].idx_j = j;
	  currlist_elements[count].idx_k = k;
	  currlist_elements[count].a_power = list[idx].ia;
	  currlist_elements[count].dummy = 0;
	  currlist_elements[count].coeff = list[idx].coeff;
	  count++;
	}
      }
    }
    list_right_simple[n1] = new hermite_conversion_element_struct[count];
    memcpy(list_right_simple[n1], &currlist_elements[0], count*sizeof(hermite_conversion_element_struct));
    counters_right_simple[n1] = count;
  }

  // Get list_left_simple
  for(int n1 = 0; n1 <= nmax; n1++) {
    get_hermite_conversion_matrix_symb(&monomial_info, n1, invFlag, &list[0]);
    int nMon1 = monomial_info.no_of_monomials_list[n1];
    int count = 0;
    for(int j = 0; j < nMon1; j++) {
      for(int k = 0; k < nMon1; k++) {
	int idx = k*nMon1+j; // reverse meaning of j k here compared to "right" case
	if(template_blas_fabs(list[idx].coeff) > 1e-5) {
	  assert(count < MAX_NO_OF_CONTRIBS);
	  currlist_elements[count].idx_j = j;
	  currlist_elements[count].idx_k = k;
	  currlist_elements[count].a_power = list[idx].ia;
	  currlist_elements[count].dummy = 0;
	  currlist_elements[count].coeff = list[idx].coeff;
	  count++;
	}
      }
    }
    list_left_simple[n1] = new hermite_conversion_element_struct[count];
    memcpy(list_left_simple[n1], &currlist_elements[0], count*sizeof(hermite_conversion_element_struct));
    counters_left_simple[n1] = count;
  }

  timeMeter.print(LOG_AREA_INTEGRALS, 
		  "hermite_conversion_info_struct constructor");
} // end hermite_conversion_info_struct constructor

void hermite_conversion_info_struct::clear_lists() {
  const int nmax = HERMITE_CONVERSION_MAX_N;
  for(int n1 = 0; n1 <= nmax; n1++)
    for(int n2 = 0; n2 <= nmax; n2++) {
      list_right[n1][n2] = NULL;
      list_left [n1][n2] = NULL;
      counters_right[n1][n2] = 0;
      counters_left[n1][n2] = 0;
    }
  for(int n1 = 0; n1 <= nmax; n1++) {
    list_right_simple[n1] = NULL;
    list_left_simple [n1] = NULL;
    counters_right_simple[n1] = 0;
    counters_left_simple[n1] = 0;
  }
}

hermite_conversion_info_struct::hermite_conversion_info_struct() {
  clear_lists();
}

hermite_conversion_info_struct::~hermite_conversion_info_struct() {
  const int nmax = HERMITE_CONVERSION_MAX_N;
  for(int n1 = 0; n1 <= nmax; n1++)
    for(int n2 = 0; n2 <= nmax; n2++) {
      delete []list_right[n1][n2];
      delete []list_left [n1][n2];
    }
  for(int n1 = 0; n1 <= nmax; n1++) {
    delete [] list_right_simple[n1];
    delete [] list_left_simple [n1];
  }
}


int hermite_conversion_info_struct::multiply_by_hermite_conversion_matrix_from_right(const monomial_info_struct & monomial_info,
										     int n1max,        
										     int n2max,        
										     ergo_real a,      
										     ergo_real* A,     
										     ergo_real* result) const
{
  int noOfContribs = counters_right[n1max][n2max];
  hermite_conversion_contrib_struct* list = list_right[n1max][n2max];
  assert(list != NULL);
  
  int nMon1 = monomial_info.no_of_monomials_list[n1max];
  int nMon2 = monomial_info.no_of_monomials_list[n2max];

  int Ntot = n1max + n2max;
  ergo_real apowlist[Ntot+1];
  apowlist[0] = 1;
  for(int i = 1; i <= Ntot; i++)
    apowlist[i] = apowlist[i-1] * a;
  
  for(int i = 0; i < nMon1*nMon2; i++)
    result[i] = 0;
  
  for(int i = 0; i < noOfContribs; i++)
    result[list[i].destIndex] += A[list[i].sourceIndex] * list[i].coeff * apowlist[-list[i].a_power];
  
  return 0;
}


int hermite_conversion_info_struct::multiply_by_hermite_conversion_matrix_from_left(const monomial_info_struct & monomial_info,
										    int n1max,        
										    int n2max,        
										    ergo_real a,      
										    ergo_real* A,     
										    ergo_real* result) const
{
  int noOfContribs = counters_left[n1max][n2max];
  hermite_conversion_contrib_struct* list = list_left[n1max][n2max];
  assert(list != NULL);
  
  int nMon1 = monomial_info.no_of_monomials_list[n1max];
  int nMon2 = monomial_info.no_of_monomials_list[n2max];
  
  int Ntot = n1max + n2max;
  ergo_real apowlist[Ntot+1];
  apowlist[0] = 1;
  for(int i = 1; i <= Ntot; i++)
    apowlist[i] = apowlist[i-1] * a;
  
  for(int i = 0; i < nMon1*nMon2; i++)
    result[i] = 0;
  
  for(int i = 0; i < noOfContribs; i++)
    result[list[i].destIndex] += A[list[i].sourceIndex] * list[i].coeff * apowlist[-list[i].a_power];
  
  return 0;
}


int hermite_conversion_info_struct::get_hermite_conversion_matrix_right(const monomial_info_struct & monomial_info,
									int nmax,
									ergo_real a,
									ergo_real* result) const {
  int noOfContribs = counters_right_simple[nmax];
  hermite_conversion_element_struct* list = list_right_simple[nmax];  
  assert(list != NULL);
  int nMon1 = monomial_info.no_of_monomials_list[nmax];
  int Ntot = 2 * nmax;
  ergo_real apowlist[Ntot+1];
  apowlist[0] = 1;
  for(int i = 1; i <= Ntot; i++)
    apowlist[i] = apowlist[i-1] * a;  
  for(int i = 0; i < nMon1*nMon1; i++)
    result[i] = 0;
  for(int i = 0; i < noOfContribs; i++) {
    int j = list[i].idx_j;
    int k = list[i].idx_k;
    result[j*nMon1+k] = list[i].coeff * apowlist[-list[i].a_power];
  }  
  return 0;
}

int hermite_conversion_info_struct::get_hermite_conversion_matrix_right_sparse(const monomial_info_struct & monomial_info,
									       int nmax,
									       ergo_real a,
									       i_j_val_struct* result) const {
  int noOfContribs = counters_right_simple[nmax];
  hermite_conversion_element_struct* list = list_right_simple[nmax];  
  assert(list != NULL);
  int Ntot = 2 * nmax;
  ergo_real apowlist[Ntot+1];
  apowlist[0] = 1;
  for(int i = 1; i <= Ntot; i++)
    apowlist[i] = apowlist[i-1] * a;
  int count = 0;
  for(int i = 0; i < noOfContribs; i++) {
    int j = list[i].idx_j;
    int k = list[i].idx_k;
    result[count].i = j;
    result[count].j = k;
    result[count].same_i_count = 1;
    result[count].value = list[i].coeff * apowlist[-list[i].a_power];
    count++;
  }
  return count;
}

int hermite_conversion_info_struct::get_hermite_conversion_matrix_left(const monomial_info_struct & monomial_info,
								       int nmax,
								       ergo_real a,
								       ergo_real* result) const {
  int noOfContribs = counters_left_simple[nmax];
  hermite_conversion_element_struct* list = list_left_simple[nmax];  
  assert(list != NULL);
  int nMon1 = monomial_info.no_of_monomials_list[nmax];
  int Ntot = 2 * nmax;
  ergo_real apowlist[Ntot+1];
  apowlist[0] = 1;
  for(int i = 1; i <= Ntot; i++)
    apowlist[i] = apowlist[i-1] * a;  
  for(int i = 0; i < nMon1*nMon1; i++)
    result[i] = 0;
  for(int i = 0; i < noOfContribs; i++) {
    int j = list[i].idx_j;
    int k = list[i].idx_k;
    result[j*nMon1+k] = list[i].coeff * apowlist[-list[i].a_power];
  }  
  return 0;
}


/** Function needed for Chunks&Tasks usage. */
hermite_conversion_info_struct::hermite_conversion_info_struct(const hermite_conversion_info_struct & other) {
  // ELIAS NOTE 2015-06-17: we start by calling clear_lists to make sure there are no pointers to garbage etc.
  // There was a segfault problem when object was deleted because the list_left_simple pointer was not allocated and not null,
  // calling clear_lists() here should avoid that.
  clear_lists();
  memcpy(counters_right, other.counters_right, sizeof(counters_right));
  memcpy(counters_left, other.counters_left, sizeof(counters_left));
  const int nmax = HERMITE_CONVERSION_MAX_N;
  for(int n1 = 0; n1 <= nmax; n1++)
    for(int n2 = 0; n2 <= nmax; n2++) {
      list_right[n1][n2] = new hermite_conversion_contrib_struct[counters_right[n1][n2]];
      memcpy(list_right[n1][n2], other.list_right[n1][n2], counters_right[n1][n2]*sizeof(hermite_conversion_contrib_struct));
      list_left[n1][n2] = new hermite_conversion_contrib_struct[counters_left[n1][n2]];
      memcpy(list_left[n1][n2], other.list_left[n1][n2], counters_left[n1][n2]*sizeof(hermite_conversion_contrib_struct));
    }
  // Take care of list_right_simple and list_left_simple
  memcpy(counters_right_simple, other.counters_right_simple, sizeof(counters_right_simple));
  memcpy(counters_left_simple, other.counters_left_simple, sizeof(counters_left_simple));
  for(int n1 = 0; n1 <= nmax; n1++) {
    list_right_simple[n1] = new hermite_conversion_element_struct[counters_right_simple[n1]];
    memcpy(list_right_simple[n1], other.list_right_simple[n1], counters_right_simple[n1]*sizeof(hermite_conversion_element_struct));
    list_left_simple[n1] = new hermite_conversion_element_struct[counters_left_simple[n1]];
    memcpy(list_left_simple[n1], other.list_left_simple[n1], counters_left_simple[n1]*sizeof(hermite_conversion_element_struct));
  }
}

/** Function needed for Chunks&Tasks usage. */
void hermite_conversion_info_struct::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error in hermite_conversion_info_struct::write_to_buffer: bufferSize too small.");
  // counters_right
  memcpy(p, counters_right, sizeof(counters_right));
  p += sizeof(counters_right);
  // counters_left
  memcpy(p, counters_left, sizeof(counters_left));
  p += sizeof(counters_left);
  // list_right and list_left
  const int nmax = HERMITE_CONVERSION_MAX_N;
  for(int n1 = 0; n1 <= nmax; n1++)
    for(int n2 = 0; n2 <= nmax; n2++) {
      memcpy(p, list_right[n1][n2], counters_right[n1][n2]*sizeof(hermite_conversion_contrib_struct));
      p += counters_right [n1][n2]*sizeof(hermite_conversion_contrib_struct);
      memcpy(p, list_left [n1][n2], counters_left [n1][n2]*sizeof(hermite_conversion_contrib_struct));
      p += counters_left  [n1][n2]*sizeof(hermite_conversion_contrib_struct);
    }
  // counters_right_simple
  memcpy(p, counters_right_simple, sizeof(counters_right_simple));
  p += sizeof(counters_right_simple);
  // counters_left_simple
  memcpy(p, counters_left_simple, sizeof(counters_left_simple));
  p += sizeof(counters_left_simple);
  // list_right_simple and list_left_simple
  for(int n1 = 0; n1 <= nmax; n1++) {
    memcpy(p, list_right_simple[n1], counters_right_simple[n1]*sizeof(hermite_conversion_element_struct));
    p += counters_right_simple[n1]*sizeof(hermite_conversion_element_struct);
    memcpy(p, list_left_simple[n1], counters_left_simple[n1]*sizeof(hermite_conversion_element_struct));
    p += counters_left_simple[n1]*sizeof(hermite_conversion_element_struct);
  }
  assert((size_t)(p - dataBuffer) == get_size());
}

/** Function needed for Chunks&Tasks usage. */
size_t hermite_conversion_info_struct::get_size() const {
  size_t size = 0;
  const int nmax = HERMITE_CONVERSION_MAX_N;
  for(int n1 = 0; n1 <= nmax; n1++)
    for(int n2 = 0; n2 <= nmax; n2++) {
      size += counters_right[n1][n2]*sizeof(hermite_conversion_contrib_struct);
      size += counters_left [n1][n2]*sizeof(hermite_conversion_contrib_struct);
    }
  for(int n1 = 0; n1 <= nmax; n1++) {
    size += counters_right_simple[n1]*sizeof(hermite_conversion_contrib_struct);
    size += counters_left_simple[n1]*sizeof(hermite_conversion_contrib_struct);
  }
  size += sizeof(counters_right);
  size += sizeof(counters_left);
  size += sizeof(counters_right_simple);
  size += sizeof(counters_left_simple);
  return size;
}

/** Function needed for Chunks&Tasks usage. */
void hermite_conversion_info_struct::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  // counters_right
  memcpy(counters_right, p, sizeof(counters_right));
  p += sizeof(counters_right);
  // counters_left
  memcpy(counters_left, p, sizeof(counters_left));
  p += sizeof(counters_left);
  // list_right and list_left
  const int nmax = HERMITE_CONVERSION_MAX_N;
  for(int n1 = 0; n1 <= nmax; n1++)
    for(int n2 = 0; n2 <= nmax; n2++) {
      list_right[n1][n2] = new hermite_conversion_contrib_struct[counters_right[n1][n2]];
      memcpy(list_right[n1][n2], p, counters_right[n1][n2]*sizeof(hermite_conversion_contrib_struct));
      p += counters_right [n1][n2]*sizeof(hermite_conversion_contrib_struct);
      list_left[n1][n2] = new hermite_conversion_contrib_struct[counters_left[n1][n2]];
      memcpy(list_left[n1][n2], p, counters_left[n1][n2]*sizeof(hermite_conversion_contrib_struct));
      p += counters_left  [n1][n2]*sizeof(hermite_conversion_contrib_struct);
    }
  // counters_right_simple
  memcpy(counters_right_simple, p, sizeof(counters_right_simple));
  p += sizeof(counters_right_simple);
  // counters_left_simple
  memcpy(counters_left_simple, p, sizeof(counters_left_simple));
  p += sizeof(counters_left_simple);
  // list_right_simple and list_left_simple
  for(int n1 = 0; n1 <= nmax; n1++) {
    list_right_simple[n1] = new hermite_conversion_element_struct[counters_right_simple[n1]];
    memcpy(list_right_simple[n1], p, counters_right_simple[n1]*sizeof(hermite_conversion_element_struct));
    p += counters_right_simple[n1]*sizeof(hermite_conversion_element_struct);
    list_left_simple[n1] = new hermite_conversion_element_struct[counters_left_simple[n1]];
    memcpy(list_left_simple[n1], p, counters_left_simple[n1]*sizeof(hermite_conversion_element_struct));
    p += counters_left_simple[n1]*sizeof(hermite_conversion_element_struct);
  }
  assert((size_t)(p - dataBuffer) == get_size());
  // DONE!
  if(static_cast<size_t>(p-dataBuffer) > bufferSize)
    throw std::runtime_error("Error: (p > bufferSize).");  
}

