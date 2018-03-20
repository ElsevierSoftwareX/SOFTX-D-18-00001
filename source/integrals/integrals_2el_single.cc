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

/** @file integrals_2el_single.cc

    @brief Functionality for computing a single 2-electron integral,
    for two given primitive Gaussian distributions.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include "integrals_2el_single.h"
#include "pi.h"
#include "boysfunction.h"

#include "integrals_hermite.h"


static ergo_real
do_2e_integral_using_symb_info_h(const JK::ExchWeights & CAM_params,
				 const DistributionSpecStruct* psi1,
				 const DistributionSpecStruct* psi2,
				 const IntegralInfo & integralInfo)
{
  const ergo_real twoTimesPiToPow5half = 2 * std::pow(pi, 2.5);
  ergo_real alpha1 = psi1->exponent;
  ergo_real alpha2 = psi2->exponent;
  ergo_real alphasum = alpha1 + alpha2;
  ergo_real alphaproduct = alpha1 * alpha2;
  ergo_real alpha0 = alphaproduct / alphasum;
  
  int n1 = 0;
  int n2 = 0;
  for(int i = 0; i < 3; i++)
    {
      n1 += psi1->monomialInts[i];
      n2 += psi2->monomialInts[i];
    }
  int n1x = psi1->monomialInts[0];
  int n1y = psi1->monomialInts[1];
  int n1z = psi1->monomialInts[2];
  int n2x = psi2->monomialInts[0];
  int n2y = psi2->monomialInts[1];
  int n2z = psi2->monomialInts[2];

  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1];
  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2];
  
  ergo_real dx0 = psi2->centerCoords[0] - psi1->centerCoords[0];
  ergo_real dx1 = psi2->centerCoords[1] - psi1->centerCoords[1];
  ergo_real dx2 = psi2->centerCoords[2] - psi1->centerCoords[2];

  ergo_real resultPreFactor = twoTimesPiToPow5half / (alphaproduct*template_blas_sqrt(alphasum));

  ergo_real primitiveIntegralList_h[noOfMonomials_1*noOfMonomials_2];
  ergo_real primitiveIntegralList_tmp[noOfMonomials_1*noOfMonomials_2];
  ergo_real primitiveIntegralList[noOfMonomials_1*noOfMonomials_2];
  
  get_related_integrals_hermite(integralInfo,
				CAM_params,
				n1, noOfMonomials_1,
				n2, noOfMonomials_2,
				dx0, 
				dx1, 
				dx2, 
				alpha0,
				resultPreFactor,
				primitiveIntegralList_h);

  integralInfo.multiply_by_hermite_conversion_matrix_from_right(n1,
								n2,
								1.0/alpha1,
								primitiveIntegralList_h,
								primitiveIntegralList_tmp);

  integralInfo.multiply_by_hermite_conversion_matrix_from_left(n1,
							       n2,
							       1.0/alpha2,
							       primitiveIntegralList_tmp,
							       primitiveIntegralList);

  int monomialIndex1 = integralInfo.monomial_info.monomial_index_list[n1x][n1y][n1z];  
  int monomialIndex2 = integralInfo.monomial_info.monomial_index_list[n2x][n2y][n2z];  

  ergo_real result = psi1->coeff * psi2->coeff * primitiveIntegralList[monomialIndex1*noOfMonomials_2+monomialIndex2];
  
  return result;
}


ergo_real
do_2e_integral_using_symb_info(const JK::ExchWeights & CAM_params,
			       const DistributionSpecStruct* psi1,
			       const DistributionSpecStruct* psi2,
			       const IntegralInfo & integralInfo)
{
  return do_2e_integral_using_symb_info_h(CAM_params, psi1, psi2, integralInfo);
}

