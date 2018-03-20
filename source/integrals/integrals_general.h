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

/** @file integrals_general.h

    @brief General functionality related to computation of integrals
    involving Gaussian basis functions.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef INTEGRALS_GENERAL_HEADER
#define INTEGRALS_GENERAL_HEADER

#include "integral_info.h"
#include "basisinfo.h"

#ifndef BASIS_FUNC_POLY_MAX_DEGREE
#error The constant BASIS_FUNC_POLY_MAX_DEGREE must be defined.
#endif
#if BASIS_FUNC_POLY_MAX_DEGREE<6
const int POLY_PRODUCT_MAX_DISTRS = 10000;
#else
const int POLY_PRODUCT_MAX_DISTRS = 20000;
#endif

typedef struct{
  ergo_real a0;
  ergo_real a1;
} polydeg1struct;

int get_product_simple_prims(const DistributionSpecStruct& primA,
			     const DistributionSpecStruct& primB,
			     DistributionSpecStruct resultList[],
			     int maxCount,
			     ergo_real threshold);

int get_product_simple_primitives(const BasisInfoStruct & basisInfoA, int iA,
				  const BasisInfoStruct & basisInfoB, int iB,
				  DistributionSpecStruct resultList[],
				  int maxCount,
				  ergo_real threshold);

ergo_real compute_integral_of_simple_prim(const DistributionSpecStruct & distr);

int multiply_polynomials(ergo_real result[], 
			 polydeg1struct* polydeg1, 
			 int dim, 
			 ergo_real a[]);

ergo_real get_largest_simple_integral(const BasisInfoStruct & basisInfo);

ergo_real get_max_basis_func_abs_value(const BasisInfoStruct & basisInfo);

int get_basis_func_extent_list(const BasisInfoStruct & basisInfo, 
			       ergo_real* basisFuncExtentList, 
			       ergo_real maxAbsValue);

#endif
