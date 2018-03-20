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

/** @file exponent_list.cc

    @brief ExponentList class for keeping track of a set of unique
    exponents in Gaussian functions, together with a max coefficient
    value for each exponent.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "exponent_list.h"
#include "output.h"


int ExponentList::get_list_of_available_exponents(const BasisInfoStruct & basisInfo)
{
  int n = basisInfo.noOfBasisFuncs;
  int count = 0;
  for(int i = 0; i < n; i++)
    {
      BasisFuncStruct* basisFunc = &basisInfo.basisFuncList[i];
      // go through all primitives for this basis function.
      int nPrims = basisFunc->noOfSimplePrimitives;
      int start  = basisFunc->simplePrimitiveIndex;
      for(int j = 0; j < nPrims; j++)
	{
	  DistributionSpecStruct* prim = &basisInfo.simplePrimitiveList[start + j];
	  ergo_real currExponent = prim->exponent;
	  ergo_real currAbsCoeff = template_blas_fabs(prim->coeff);
	  // now go through list to check if we already have this exponent.
	  int foundIndex = -1;
	  for(int k = 0; k < count; k++)
	    {
	      ergo_real absDiff = template_blas_fabs(list[k].exponent - currExponent);
	      if(absDiff < CONST_EXPONENT_DIFF_TOLERANCE)
		{
		  foundIndex = k;
		  break;
		} // END IF found
	    } // END FOR k
	  if(foundIndex >= 0)
	    {
	      // OK, we already have this exponent in list.
	      // Update maxAbsCoeff if needed.
	      if(currAbsCoeff > list[foundIndex].maxAbsCoeff)
		list[foundIndex].maxAbsCoeff = currAbsCoeff;
	    }
	  else
	    {
	      // Add new exponent to list.
	      if(count >= MAX_NO_OF_UNIQUE_EXPONENTS)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_available_exponents: (count >= MAX_NO_OF_UNIQUE_EXPONENTS)");
		  return -1;
		}
	      list[count].exponent = currExponent;
	      list[count].maxAbsCoeff = currAbsCoeff;
	      count++;
	    }
	} // END FOR j
    } // END FOR i
  noOfExponents = count;
  return 0;
}
