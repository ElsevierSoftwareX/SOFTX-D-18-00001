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

/** @file basis_func_pair_list_1el.h

    @brief Functions for setting up lists of non-negligible basis
    function pairs, for 1-electron integrals.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef BASIS_FUNC_PAIR_LIST_1EL_HEADER
#define BASIS_FUNC_PAIR_LIST_1EL_HEADER


#include "realtype.h"
#include "basisinfo.h"



typedef struct
{
  int index_1;
  int index_2;
} basis_func_index_pair_struct_1el;

int
get_basis_func_pair_list_1el_for_V(const BasisInfoStruct & basisInfo,
				   ergo_real threshold,
				   ergo_real boxSize,
				   ergo_real maxCharge,
				   basis_func_index_pair_struct_1el* result_basisFuncPairList,
				   int resultMaxCount);

int
get_basis_func_pair_list_simple(const BasisInfoStruct & basisInfo,
				ergo_real threshold,
				ergo_real boxSize,
				basis_func_index_pair_struct_1el* result_basisFuncPairList,
				int resultMaxCount);


#endif
