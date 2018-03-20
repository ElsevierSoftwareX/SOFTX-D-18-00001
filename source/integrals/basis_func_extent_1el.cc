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

/** @file basis_func_extent_1el.cc

    @brief Code for determining extent of basis functions, for
    1-electron integral evaluation.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basis_func_extent_1el.h"
#include "output.h"
#include "pi.h"
#include "exponent_list.h"


int
compute_extent_for_all_basis_funcs_1el(const BasisInfoStruct & basisInfo, 
				       ergo_real* basisFuncExtentList, 
				       ergo_real maxCharge,
				       ergo_real threshold)
{
  // Create discretized list of exponents with maxAbsCoeff for each unique exponent.
  ExponentList exponentList;
  if(exponentList.get_list_of_available_exponents(basisInfo) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_available_exponents");
      return -1;
    }
  int n = basisInfo.noOfBasisFuncs;
  for(int i = 0; i < n; i++)
    {
      BasisFuncStruct* basisFunc = &basisInfo.basisFuncList[i];
      ergo_real largestExtentSoFar = 0;
      // go through all primitives for this basis function.
      int nPrims = basisFunc->noOfSimplePrimitives;
      int start  = basisFunc->simplePrimitiveIndex;
      for(int j = 0; j < nPrims; j++)
	{
	  DistributionSpecStruct* prim = &basisInfo.simplePrimitiveList[start + j];
	  ergo_real currExponent = prim->exponent;
	  ergo_real currAbsCoeff = template_blas_fabs(prim->coeff);
	  ergo_real a = currExponent;
	  ergo_real c_a = currAbsCoeff;
	  // now go through all available exponents
	  for(int ii = 0; ii < exponentList.noOfExponents; ii++)
	    {
	      ergo_real b = exponentList.list[ii].exponent;
	      ergo_real c_b = exponentList.list[ii].maxAbsCoeff;
	      if(c_b > 0)
		{
		  // This extent definition was derived from the estimate of the largest possible contribution to V given by (2 * pi * coeff / exponent).
		  ergo_real R2 = -1 * ((a+b)/(a*b)) * template_blas_log(threshold*(a+b)/(2*pi*c_a*c_b*maxCharge));
		  if(R2 < 0)
		    continue;
		  ergo_real R = template_blas_sqrt(R2);
		  if(R > largestExtentSoFar)
		    largestExtentSoFar = R;
		}
	    } // END FOR ii
	} // END FOR j
      if( largestExtentSoFar <= 0 )
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_extent_for_all_basis_funcs_1el: (largestExtentSoFar <= 0)");
	  return -1;
	}
      basisFuncExtentList[i] = largestExtentSoFar;
    } // END FOR i
  return 0;
}

