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

/** @file integrals_hermite.h

    @brief Code for computation of Coulomb integrals of Hermite
    Gaussians, using the the McMurchie-Davidson scheme as described in
    the book "Molecular electronic-structure theory" by Trygve
    Helgaker, Poul Jorgensen, and Jeppe Olsen.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef INTEGRALS_HERMITE_HEADER
#define INTEGRALS_HERMITE_HEADER

#include "realtype.h"
#include "integral_info.h"

int get_related_integrals_hermite(const IntegralInfo & integralInfo,
				  const JK::ExchWeights & paramsCAM,
				  int n1max, int noOfMonomials_1,
				  int n2max, int noOfMonomials_2,
				  ergo_real dx0, 
				  ergo_real dx1, 
				  ergo_real dx2, 
				  ergo_real alpha0,
				  ergo_real resultPreFactor,
				  ergo_real* primitiveIntegralList);


#endif
