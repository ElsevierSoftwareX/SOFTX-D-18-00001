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

/** @file basis_func_extent.h

    @brief Code for determining extent of basis functions, for
    2-electron integral evaluation.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef BASIS_FUNC_EXTENT_HEADER
#define BASIS_FUNC_EXTENT_HEADER


#include "realtype.h"
#include "integral_info.h"
#include "basisinfo.h"


int
compute_extent_for_all_basis_funcs_2el(const IntegralInfo & integralInfo,
				       const BasisInfoStruct & basisInfo, 
				       ergo_real* basisFuncExtentList, 
				       ergo_real threshold,
				       ergo_real maxAbsDensMatElement);

#endif
