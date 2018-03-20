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

/** @file operator_matrix.h

    @brief Functions for computing the matrix of a
    dipole/quadrupole/etc operator. Full and sparse versions.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef OPERATOR_MATRIX_HEADER
#define OPERATOR_MATRIX_HEADER

#include "integral_info.h"
#include "basisinfo.h"


int compute_overlap_matrix(const BasisInfoStruct & basisInfoA,
			   const BasisInfoStruct & basisInfoB,
			   ergo_real* result);

int compute_operator_matrix_full(const BasisInfoStruct & basisInfoA,
				 const BasisInfoStruct & basisInfoB,
				 int pow_x,
				 int pow_y,
				 int pow_z,
				 ergo_real* result);

int compute_operator_matrix_sparse(const BasisInfoStruct & basisInfoA,
				   const BasisInfoStruct & basisInfoB,
				   int pow_x,
				   int pow_y,
				   int pow_z,
				   int n_A,
				   int n_B,
				   std::vector<int> & nvaluesList,
				   std::vector< std::vector<int> > & colindList,
				   std::vector< std::vector<ergo_real> > & valuesList);


#endif
