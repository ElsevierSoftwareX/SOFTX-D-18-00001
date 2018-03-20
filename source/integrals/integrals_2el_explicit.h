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

/** @file integrals_2el_explicit.h

    @brief Code for explicit computation of 4-index 2-electron
    integrals.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef INTEGRALS_2EL_EXPLICIT_HEADER
#define INTEGRALS_2EL_EXPLICIT_HEADER


#include "basisinfo.h"

ergo_real do_2e_integral(int mu, 
			 int nu, 
			 int la, 
			 int si, 
			 const BasisInfoStruct & basisInfo, 
			 const IntegralInfo & integralInfo);

ergo_real do_2e_integral_general(int mu, 
				 int nu, 
				 int la, 
				 int si, 
				 const BasisInfoStruct & basisInfo_mu, 
				 const BasisInfoStruct & basisInfo_nu, 
				 const BasisInfoStruct & basisInfo_la, 
				 const BasisInfoStruct & basisInfo_si, 
				 const IntegralInfo & integralInfo);

int compute_2e_matrix_list_explicit(const BasisInfoStruct & basisInfo,
				    const IntegralInfo & integralInfo,
				    ergo_real** resultList,
				    ergo_real** densList,
				    int noOfMatrices,
				    ergo_real threshold);

int compute_2e_matrix_simple(const BasisInfoStruct & basisInfo,
			     const IntegralInfo & integralInfo,
			     ergo_real hf_weight,
			     ergo_real* result,
			     const ergo_real* dens);


#endif /* INTEGRALS_2EL_EXPLICIT_HEADER */
