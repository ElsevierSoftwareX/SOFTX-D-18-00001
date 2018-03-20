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

/** @file integrals_2el_repeating.cc

    @brief Functionality for keeping track of certain kinds of
    integrals that are computed repeatedly, saving the computed values
    instead of recomputing them.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include <assert.h>
#include "pi.h"
#include "integrals_hermite.h"
#include "integrals_2el_repeating.h"
#include "realtype.h"


IntegratorCase::IntegratorCase(int Nmax_in, 
			       int noOfMonomials, 
			       ergo_real exponent_in, 
			       const ergo_real* newList)
{
  list = new ergo_real[noOfMonomials];
  Nmax = Nmax_in;
  exponent = exponent_in;
  for(int i = 0; i < noOfMonomials; i++)
    list[i] = newList[i];
}

IntegratorCase::~IntegratorCase()
{
  delete []list;
}


typedef IntegratorCase* IntegratorCasePtr;

const int MAX_NO_OF_CASES = 44444;

IntegratorWithMemory::IntegratorWithMemory(const IntegralInfo* b)
{
  integralInfo = b;
  noOfCases = 0;
  caseList = new IntegratorCasePtr[MAX_NO_OF_CASES];
}

IntegratorWithMemory::~IntegratorWithMemory()
{
  for(int i = 0; i < noOfCases; i++)
    delete caseList[i];
  delete []caseList;
}

ergo_real
IntegratorWithMemory::do_2e_integral(const DistributionSpecStruct* psi)
{
  const ergo_real twoTimesPiToPow5half = 2 * pitopow52;
  int Ntot;
  ergo_real alpha1 = psi->exponent;
  ergo_real alpha2 = psi->exponent;
  ergo_real alphasum = alpha1 + alpha2;
  ergo_real alphaproduct = alpha1 * alpha2;
  ergo_real alpha0 = alphaproduct / alphasum;
  
  int n1 = 0;
  int n2 = 0;
  for(int i = 0; i < 3; i++)
    {
      n1 += psi->monomialInts[i];
      n2 += psi->monomialInts[i];
    }
  int nx = psi->monomialInts[0];
  int ny = psi->monomialInts[1];
  int nz = psi->monomialInts[2];

  Ntot = n1 + n2;
  int Nmax = Ntot;

  int noOfMonomials = integralInfo->monomial_info.no_of_monomials_list[n1];

  ergo_real dx0 = 0;
  ergo_real dx1 = 0;
  ergo_real dx2 = 0;

  ergo_real resultPreFactor = twoTimesPiToPow5half / (alphaproduct*template_blas_sqrt(alphasum));

  int monomialIndex = integralInfo->monomial_info.monomial_index_list[nx][ny][nz];  

  // Check if this type of integral is already known
  // That is, have we already calculated integrals for this Nmax (or higher Nmax) and for this exponent?
  int foundIndex = -1;
  if(noOfCases > 0)
    {
      int lo = 0;
      int hi = noOfCases - 1;
      while(lo < hi-1)
	{
	  int mid = (lo + hi) / 2;
	  if(caseList[mid]->exponent < alpha0)
	    lo = mid;
	  else
	    hi = mid;
	} // END WHILE
      const ergo_real exponent_diff_limit = 1e-11; // FIXME: SHOULD NOT USE HARD-CODED VALUE HERE!
      ergo_real exponentDiff1 = template_blas_fabs(caseList[lo]->exponent - alpha0);
      if(exponentDiff1 < exponent_diff_limit)
	foundIndex = lo;
      ergo_real exponentDiff2 = template_blas_fabs(caseList[hi]->exponent - alpha0);
      if(exponentDiff2 < exponent_diff_limit)
	foundIndex = hi;
    }
  if(foundIndex >= 0)
    {
      if(caseList[foundIndex]->Nmax >= Nmax)
	{
	  // OK, found it!
	  return resultPreFactor * psi->coeff * psi->coeff * caseList[foundIndex]->list[monomialIndex];
	}
    }


  // No, not found. Create new case.

  ergo_real primitiveIntegralList_h[noOfMonomials*noOfMonomials];
  ergo_real primitiveIntegralList_tmp[noOfMonomials*noOfMonomials];
  ergo_real primitiveIntegralList[noOfMonomials*noOfMonomials];

  const JK::ExchWeights CAM_params_not_used;
  
  get_related_integrals_hermite(*integralInfo,
				CAM_params_not_used,
				n1, noOfMonomials,
				n2, noOfMonomials,
				dx0, 
				dx1, 
				dx2, 
				alpha0,
				1.0,
				primitiveIntegralList_h);

  integralInfo->multiply_by_hermite_conversion_matrix_from_right(n1,
								 n2,
								 1.0/alpha1,
								 primitiveIntegralList_h,
								 primitiveIntegralList_tmp);

  integralInfo->multiply_by_hermite_conversion_matrix_from_left(n1,
								n2,
								1.0/alpha2,
								primitiveIntegralList_tmp,
								primitiveIntegralList);

  ergo_real newList[noOfMonomials];
  for(int i = 0; i < noOfMonomials; i++)
    newList[i] = primitiveIntegralList[i*noOfMonomials+i];
  
  assert(noOfCases < MAX_NO_OF_CASES);


  IntegratorCase* newCase = new IntegratorCase(Nmax, noOfMonomials, alpha0, newList);

  // Check if this exponent is already present in list. If so, replace it with this new higher Nmax.
  if(foundIndex >= 0)
    {
      delete caseList[foundIndex];
      caseList[foundIndex] = newCase;
    }
  else
    {
      // Insert new entry at proper place in list.
      // First skip all entries with too small exponent.
      int nSkipped = 0;
      for(int i = 0; i < noOfCases; i++)
	{
	  if(caseList[i]->exponent < alpha0)
	    nSkipped++;
	  else
	    break;
	} // END FOR i 
      int newIndex = nSkipped;
      // Now move all remaining entries one step down in list.
      for(int i = noOfCases-1; i >= nSkipped; i--)
	caseList[i+1] = caseList[i];
      // Now add new entry.
      caseList[newIndex] = newCase;
      noOfCases++;
    }

  // We have now modified the list so that the needed case is present. Call this function again.
  return do_2e_integral(psi);
}


