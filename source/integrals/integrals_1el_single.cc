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

/** @file integrals_1el_single.cc

    @brief Functionality for computing a single 1-electron integral,
    for a given primitive Gaussian distribution and a given point
    charge.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include "integrals_1el_single.h"
#include "pi.h"
#include "boysfunction.h"

#include "integrals_hermite.h"


static ergo_real
do_1e_repulsion_integral_using_symb_info_h(const DistributionSpecStruct & psi,
					   ergo_real pointCharge,
					   const ergo_real* pointChargeCoords,
					   const IntegralInfo & integralInfo)
{
  // Let the distr be 1 and the pointcharge 2
  ergo_real alpha1 = psi.exponent;
  // alpha2 is ~ infinity
  ergo_real alpha0 = alpha1;
  // alpham is ~ infinity
  // g is 1 (not needed)
  
  int n1 = 0;
  int n2 = 0;
  for(int i = 0; i < 3; i++)
    n1 += psi.monomialInts[i];
  int n1x = psi.monomialInts[0];
  int n1y = psi.monomialInts[1];
  int n1z = psi.monomialInts[2];

  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1];
  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2];

  ergo_real dx0 = pointChargeCoords[0] - psi.centerCoords[0];
  ergo_real dx1 = pointChargeCoords[1] - psi.centerCoords[1];
  ergo_real dx2 = pointChargeCoords[2] - psi.centerCoords[2];

  ergo_real resultPreFactor = 2 * pi / alpha1;

  ergo_real primitiveIntegralList_h[noOfMonomials_1*noOfMonomials_2];
  ergo_real primitiveIntegralList_2[noOfMonomials_1*noOfMonomials_2];

  const JK::ExchWeights CAM_params_not_used;
  
  get_related_integrals_hermite(integralInfo,
				CAM_params_not_used,
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
								primitiveIntegralList_2);
  int monomialIndex = integralInfo.monomial_info.monomial_index_list[n1x][n1y][n1z];

  ergo_real result = psi.coeff * pointCharge * primitiveIntegralList_2[monomialIndex];

  return result;
}


/* This routine is supposed to compute derivatives of integrals w.r.t. changes in the pointCharge coordinates.  */
std::vector<ergo_real> do_1e_repulsion_integral_derivatives_using_symb_info(const DistributionSpecStruct* psi,
									    ergo_real pointCharge,
									    const ergo_real* pointChargeCoords,
									    const IntegralInfo & integralInfo) {
  // Let the distr be 1 and the pointcharge 2
  ergo_real alpha1 = psi->exponent;
  // alpha2 is ~ infinity
  ergo_real alpha0 = alpha1;
  // alpham is ~ infinity
  // g is 1 (not needed)
  
  int n1 = 0;
  int n2 = 1;
  for(int i = 0; i < 3; i++)
    n1 += psi->monomialInts[i];
  int n1x = psi->monomialInts[0];
  int n1y = psi->monomialInts[1];
  int n1z = psi->monomialInts[2];

  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1];
  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2];

  ergo_real dx0 = pointChargeCoords[0] - psi->centerCoords[0];
  ergo_real dx1 = pointChargeCoords[1] - psi->centerCoords[1];
  ergo_real dx2 = pointChargeCoords[2] - psi->centerCoords[2];

  ergo_real resultPreFactor = 2 * pi / alpha1;

  ergo_real primitiveIntegralList_h[noOfMonomials_1*noOfMonomials_2];

  const JK::ExchWeights CAM_params_not_used;
  
  get_related_integrals_hermite(integralInfo,
				CAM_params_not_used,
				n1, noOfMonomials_1,
				n2, noOfMonomials_2,
				dx0, 
				dx1, 
				dx2, 
				alpha0,
				resultPreFactor,
				primitiveIntegralList_h);

  int n1b = n1;
  int n2b = 0;
  ergo_real primitiveIntegralList_h_components[3][noOfMonomials_1];
  int monomialIndex_x = integralInfo.monomial_info.monomial_index_list[1][0][0];
  int monomialIndex_y = integralInfo.monomial_info.monomial_index_list[0][1][0];
  int monomialIndex_z = integralInfo.monomial_info.monomial_index_list[0][0][1];
  for(int i = 0; i < noOfMonomials_1; i++) {
    primitiveIntegralList_h_components[0][i] = primitiveIntegralList_h[i*noOfMonomials_2+monomialIndex_x];
    primitiveIntegralList_h_components[1][i] = primitiveIntegralList_h[i*noOfMonomials_2+monomialIndex_y];
    primitiveIntegralList_h_components[2][i] = primitiveIntegralList_h[i*noOfMonomials_2+monomialIndex_z];
  }
  ergo_real primitiveIntegralList_2_components[3][noOfMonomials_1];
  for(int i = 0; i < 3; i++)
    integralInfo.multiply_by_hermite_conversion_matrix_from_right(n1b, n2b, 1.0/alpha1, primitiveIntegralList_h_components[i], primitiveIntegralList_2_components[i]);

  int monomialIndex = integralInfo.monomial_info.monomial_index_list[n1x][n1y][n1z];

  ergo_real result_x = psi->coeff * pointCharge * primitiveIntegralList_2_components[0][monomialIndex];
  ergo_real result_y = psi->coeff * pointCharge * primitiveIntegralList_2_components[1][monomialIndex];
  ergo_real result_z = psi->coeff * pointCharge * primitiveIntegralList_2_components[2][monomialIndex];

  std::vector<ergo_real> v(3);
  v[0] = result_x;
  v[1] = result_y;
  v[2] = result_z;
  return v;
}





ergo_real 
do_1e_repulsion_integral_using_symb_info(const DistributionSpecStruct & psi,
					 ergo_real pointCharge,
					 const ergo_real* pointChargeCoords,
					 const IntegralInfo & integralInfo)
{
  return do_1e_repulsion_integral_using_symb_info_h(psi,
						    pointCharge,
						    pointChargeCoords,
						    integralInfo);
}



