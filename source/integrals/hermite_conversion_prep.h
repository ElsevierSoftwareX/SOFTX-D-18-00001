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

/** @file hermite_conversion_prep.h

    @brief Code for conversion between integrals computed for Hermite
    Gaussians and Cartesian Gaussians, using saved information from a
    preparatory step that can be done once and re-used many times.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef HERMITE_CONVERSION_PREP_HEADER
#define HERMITE_CONVERSION_PREP_HEADER

#include <cstring>
#include "realtype.h"
#include "polydegree.h"
#include "monomial_info.h"
#include "simple_sparse_mat.h"

struct hermite_conversion_contrib_struct
{
  int destIndex;
  int sourceIndex;
  int a_power;
  int dummy;
  ergo_real coeff;
};

typedef hermite_conversion_contrib_struct * hermite_conversion_contrib_struct_ptr;

struct hermite_conversion_element_struct {
  int idx_j;
  int idx_k;
  int a_power;
  int dummy;
  ergo_real coeff;
};

typedef hermite_conversion_element_struct * hermite_conversion_element_struct_ptr;

const int HERMITE_CONVERSION_MAX_N = BASIS_FUNC_POLY_MAX_DEGREE*2;

class hermite_conversion_info_struct {
 private:
  hermite_conversion_contrib_struct_ptr list_right[HERMITE_CONVERSION_MAX_N+1][HERMITE_CONVERSION_MAX_N+1];
  hermite_conversion_contrib_struct_ptr list_left [HERMITE_CONVERSION_MAX_N+1][HERMITE_CONVERSION_MAX_N+1];
  int counters_right[HERMITE_CONVERSION_MAX_N+1][HERMITE_CONVERSION_MAX_N+1];
  int counters_left [HERMITE_CONVERSION_MAX_N+1][HERMITE_CONVERSION_MAX_N+1];
  // Simple lists used to represent conversion matrices
  hermite_conversion_element_struct_ptr list_right_simple[HERMITE_CONVERSION_MAX_N+1];
  hermite_conversion_element_struct_ptr list_left_simple [HERMITE_CONVERSION_MAX_N+1];
  int counters_right_simple[HERMITE_CONVERSION_MAX_N+1];
  int counters_left_simple [HERMITE_CONVERSION_MAX_N+1];
  void clear_lists();
 public:
  void init(const monomial_info_struct & monomial_info);
  hermite_conversion_info_struct();
  ~hermite_conversion_info_struct();
  int multiply_by_hermite_conversion_matrix_from_right(const monomial_info_struct & monomial_info,
						       int n1max,        
						       int n2max,        
						       ergo_real a,      
						       ergo_real* A,     
						       ergo_real* result) const;
  int multiply_by_hermite_conversion_matrix_from_left(const monomial_info_struct & monomial_info,
						      int n1max,        
						      int n2max,        
						      ergo_real a,      
						      ergo_real* A,     
						      ergo_real* result) const;
  int get_hermite_conversion_matrix_right(const monomial_info_struct & monomial_info,
					  int nmax,
					  ergo_real a,
					  ergo_real* result) const;
  int get_hermite_conversion_matrix_left(const monomial_info_struct & monomial_info,
					 int nmax,
					 ergo_real a,
					 ergo_real* result) const;

  int get_hermite_conversion_matrix_right_sparse(const monomial_info_struct & monomial_info,
						 int nmax,
						 ergo_real a,
						 i_j_val_struct* result) const;

  // Stuff needed for Chunks&Tasks usage
  hermite_conversion_info_struct(const hermite_conversion_info_struct & other);
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};


#endif
