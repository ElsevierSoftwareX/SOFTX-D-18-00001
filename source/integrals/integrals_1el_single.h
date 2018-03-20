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

/** @file integrals_1el_single.h

    @brief Functionality for computing a single 1-electron integral,
    for a given primitive Gaussian distribution and a given point
    charge.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef INTEGRALS_1EL_SINGLE_HEADER
#define INTEGRALS_1EL_SINGLE_HEADER


#include "basisinfo.h"

ergo_real do_1e_repulsion_integral_using_symb_info(const DistributionSpecStruct & psi,
						   ergo_real pointCharge,
						   const ergo_real* pointChargeCoords,
						   const IntegralInfo & integralInfo);

std::vector<ergo_real> do_1e_repulsion_integral_derivatives_using_symb_info(const DistributionSpecStruct* psi,
									    ergo_real pointCharge,
									    const ergo_real* pointChargeCoords,
									    const IntegralInfo & integralInfo);

ergo_real do_1e_repulsion_integral_list_using_symb_info(DistributionSpecStruct* psi,
							const ergo_real* pointChargeList,
							const ergo_real* pointChargeCoordsList,
							int noOfPointCharges,
							const IntegralInfo & integralInfo);


#endif
