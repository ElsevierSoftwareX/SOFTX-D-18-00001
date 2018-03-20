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

/** @file hermite_conversion_explicit.cc

    @brief Code for explicitly computing the matrix for conversion
    between integrals computed for Hermite Gaussians and Cartesian
    Gaussians, for given nmax and exponent values.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <cmath>
#include "hermite_conversion_symb.h"
#include "mat_gblas.h"

int 
get_hermite_conversion_matrix(const monomial_info_struct* monomial_info,
			      int nmax,
			      int inverseFlag,
			      ergo_real exponent,
			      ergo_real* result)
{
  int noOfMonomials = monomial_info->no_of_monomials_list[nmax];
  std::vector<symb_matrix_element> result_symb(noOfMonomials*noOfMonomials);
  get_hermite_conversion_matrix_symb(monomial_info,
				     nmax,
				     inverseFlag,
				     &result_symb[0]);
  for(int i = 0; i < noOfMonomials*noOfMonomials; i++)
    {
      ergo_real factor = template_blas_pow(exponent, (ergo_real)(result_symb[i].ia));
      result[i] = result_symb[i].coeff * factor;
    }
  return 0;
}
