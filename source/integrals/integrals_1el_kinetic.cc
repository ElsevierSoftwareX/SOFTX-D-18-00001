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

/** @file integrals_1el_kinetic.cc

    @brief Code for 1-electron integrals, computation of
    kinetic-energy matrix T.

    @author: Elias Rudberg <em>responsible</em>
*/

/* Written by Elias Rudberg, KTH, Stockholm */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>

#include "integrals_1el_kinetic.h"
#include "memorymanag.h"
#include "pi.h"
#include "output.h"
#include "utilities.h"
#include "boysfunction.h"
#include "integral_info.h"
#include "integrals_general.h"
#include "box_system.h"
#include "multipole.h"
#include "integrals_2el_single.h"
#include "integrals_1el_single.h"
#include "basis_func_pair_list_1el.h"


/* FIXME do not use this hard-coded value! */
static const ergo_real MATRIX_ELEMENT_THRESHOLD_VALUE = 1e-12;



static void
do_derivative_of_simple_prim(const DistributionSpecStruct& prim,
			     DistributionSpecStruct* resultList,
			     int coord)
{
  /* first term */
  if(prim.monomialInts[coord] > 0)
    {
      memcpy(&resultList[0], &prim, sizeof(DistributionSpecStruct));
      resultList[0].coeff *= prim.monomialInts[coord];
      resultList[0].monomialInts[coord] -= 1;
    }
  else
    {
      /* first term is zero */
      resultList[0].coeff = 0;
    }
  /* second term */
  memcpy(&resultList[1], &prim, sizeof(DistributionSpecStruct));
  resultList[1].coeff *= -2*prim.exponent;
  resultList[1].monomialInts[coord] += 1;
}

/** Computes the contribution to kinetic energy integral along the
   cartesian coordinate coord between two distributions prim1 and
   prim2. Note that this function is *not* strict wrt the
   effectiveThreshold parameter, the approximation is only
   proportional to its value but it can exceed it. */
ergo_real 
simplePrimTintegral(const DistributionSpecStruct& prim1,
		    const DistributionSpecStruct& prim2,
		    int coord,
		    ergo_real threshold)
{
  const int maxDistrsInTempList = 888;
  DistributionSpecStruct tempList[maxDistrsInTempList];
  int i, k, nNewPrims;
  ergo_real sum;
  DistributionSpecStruct list1[2];
  DistributionSpecStruct list2[4];
  do_derivative_of_simple_prim(prim2, list1, coord);
  if(list1[0].coeff != 0)
    {
      do_derivative_of_simple_prim(list1[0], &list2[0], coord);
    }
  else
    {
      list2[0].coeff = 0;
      list2[1].coeff = 0;
    }
  if(list1[1].coeff != 0)
    {
      do_derivative_of_simple_prim(list1[1], &list2[2], coord);
    }
  else
    {
      list2[2].coeff = 0;
      list2[3].coeff = 0;
    }
  /* now the resulting 4 terms are stored in list2 */

  sum = 0;
  for(i = 0; i < 4; i++)
    {
      if(list2[i].coeff == 0)
	continue;

      nNewPrims = get_product_simple_prims(prim1, 
					   list2[i], 
					   tempList,
					   maxDistrsInTempList,
					   threshold);
      if(nNewPrims < 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_prims");
	  return -1;
	}

      for(k = 0; k < nNewPrims; k++)
	{
	  const DistributionSpecStruct & currDistr = tempList[k];
	  sum += compute_integral_of_simple_prim(currDistr);
	} /* END FOR k */
    }
  
  return sum;
}


int
compute_T_matrix_sparse_linear(const BasisInfoStruct& basisInfo,
			       ergo_real threshold,
			       ergo_real boxSize,
			       int* nvaluesList,
			       int** colindList,
			       ergo_real** valuesList)
{
  int internal_error = 0;
  int n = basisInfo.noOfBasisFuncs;

  int noOfBasisFuncIndexPairs = get_basis_func_pair_list_simple(basisInfo, threshold, boxSize, NULL, 2000000000);
  if(noOfBasisFuncIndexPairs <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in get_basis_func_pair_list_simple, noOfBasisFuncIndexPairs = %i", noOfBasisFuncIndexPairs);
    return -1;
  }
  std::vector<basis_func_index_pair_struct_1el> basisFuncIndexPairList(noOfBasisFuncIndexPairs);
  noOfBasisFuncIndexPairs = get_basis_func_pair_list_simple(basisInfo, threshold, boxSize, &basisFuncIndexPairList[0], noOfBasisFuncIndexPairs);
  if(noOfBasisFuncIndexPairs <= 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in get_basis_func_pair_list_simple, noOfBasisFuncIndexPairs = %i", noOfBasisFuncIndexPairs);
    return -1;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "compute_T_matrix_sparse_linear: n = %d, threshold = %g, boxSize = %f",
	    n, (double)threshold, (double)boxSize);
  do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "compute_T_matrix_sparse_linear: noOfBasisFuncIndexPairs = %i ==> storing %6.2f %% of a full matrix", 
	    noOfBasisFuncIndexPairs, (double)100*noOfBasisFuncIndexPairs/((double)n*n));

  // To reduce scaling we want some kind of "extent" for each basis function.
  // Start by getting largest simple integral for each of the two basis sets.
  ergo_real A = get_largest_simple_integral(basisInfo);
  std::vector<ergo_real> basisFuncExtentList(n);
  get_basis_func_extent_list(basisInfo, &basisFuncExtentList[0], MATRIX_ELEMENT_THRESHOLD_VALUE / A);

  std::vector<int> offsetVec(n);
  std::vector<int> countVec(n);
  int currOffset = 0;
  int countSumToVerify = 0;
  for(int i = 0; i < n; i++) {
    int savedOffset = currOffset;
    while(currOffset < noOfBasisFuncIndexPairs && basisFuncIndexPairList[currOffset].index_1 == i)
      currOffset++;
    int count = currOffset - savedOffset;
    offsetVec[i] = savedOffset;
    countVec[i] = count;
    countSumToVerify += count;
  }
  assert(currOffset == noOfBasisFuncIndexPairs);
  assert(countSumToVerify == noOfBasisFuncIndexPairs);
  
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Allocate vector for results for one row.
    std::vector<ergo_real> rowValueList(n);
    std::vector<int> row_nu_list(n);
#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
    for(int mu = 0; mu < n; mu++) {
      int no_of_nu_values = countVec[mu];
      int startOffset = offsetVec[mu];
      int count = 0;
      BasisFuncStruct* basisFunc_mu = &basisInfo.basisFuncList[mu];
      int n_mu = basisFunc_mu->noOfSimplePrimitives;
      int start_prim_mu = basisFunc_mu->simplePrimitiveIndex;
      DistributionSpecStruct* list_mu = &basisInfo.simplePrimitiveList[start_prim_mu];
      for(int nuCounter = 0; nuCounter < no_of_nu_values; nuCounter++) {
	int nu = basisFuncIndexPairList[startOffset+nuCounter].index_2;
	assert(mu == basisFuncIndexPairList[startOffset+nuCounter].index_1);
	assert(nu <= mu);
	// Compute distance between basis function centers
	ergo_real dx = basisInfo.basisFuncList[mu].centerCoords[0] - basisInfo.basisFuncList[nu].centerCoords[0];
	ergo_real dy = basisInfo.basisFuncList[mu].centerCoords[1] - basisInfo.basisFuncList[nu].centerCoords[1];
	ergo_real dz = basisInfo.basisFuncList[mu].centerCoords[2] - basisInfo.basisFuncList[nu].centerCoords[2];
	ergo_real distance = template_blas_sqrt(dx*dx + dy*dy + dz*dz);
	// We can skip if distance is greater than sum of extents.
	if(distance > basisFuncExtentList[mu] + basisFuncExtentList[nu])
	  continue;
	BasisFuncStruct* basisFunc_nu = &basisInfo.basisFuncList[nu];
	int n_nu = basisFunc_nu->noOfSimplePrimitives;
	int start_prim_nu = basisFunc_nu->simplePrimitiveIndex;
	DistributionSpecStruct* list_nu = &basisInfo.simplePrimitiveList[start_prim_nu];
	/* compute matrix element [mu,nu] */
	ergo_real sum = 0;
	int i, j, k;
	for(j = 0; j < n_mu; j++) {
	  const DistributionSpecStruct& prim_mu_j = list_mu[j];
	  for(k = 0; k < n_nu; k++) {
	    const DistributionSpecStruct& prim_nu_k = list_nu[k];
	    ergo_real effectiveThreshold = 2.0*threshold/(n_mu*n_nu*3);
	    /* now loop over coordinates */
	    for(i = 0; i < 3; i++) {
	      /* Note that this function is not strict wrt the
		 effectiveThreshold parameter, the
		 approximation is only proportional to its
		 value but it can exceed it. */
	      sum += simplePrimTintegral(prim_mu_j,
					 prim_nu_k,
					 i,
					 effectiveThreshold);
	    } /* END FOR i */
	  } /* END FOR k */
	} /* END FOR j */
	rowValueList[count] = -0.5 * sum;
	row_nu_list[count] = nu;
	if(template_blas_fabs(rowValueList[count]) > MATRIX_ELEMENT_THRESHOLD_VALUE)
	  count++;
      } /* END FOR nuCounter */
      // OK, this row done.
      // Now go through results to check which elements need to be saved.
      nvaluesList[mu] = count;
      // Now allocate result vectors for this row.
      colindList[mu] = ergo_new(count, int);
      valuesList[mu] = ergo_new(count, ergo_real);
      for(int j = 0; j < count; j++) {
	colindList[mu][j] = row_nu_list[j];
	valuesList[mu][j] = rowValueList[j];
      }
    } /* END FOR mu */
  }

  return internal_error ? -1 : 0;
}



int
compute_T_matrix_full(const BasisInfoStruct& basisInfo,
		      ergo_real threshold,
		      ergo_real* result)
{
  int n = basisInfo.noOfBasisFuncs;
  int* nvaluesList = ergo_new(n, int);
  int** colindList = ergo_new(n, int*);
  ergo_real** valuesList = ergo_new(n, ergo_real*);

  ergo_real boxSize = 6.3;
  if(compute_T_matrix_sparse_linear(basisInfo,
				    threshold,
				    boxSize,
				    nvaluesList,
				    colindList,
				    valuesList) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_T_matrix_sparse");
      return -1;
    }
  
  // Now populate full result matrix
  memset(result, 0, n*n*sizeof(ergo_real));
  int i;
  for(i = 0; i < n; i++)
    {
      int count = nvaluesList[i];
      int* colind = colindList[i];
      ergo_real* values = valuesList[i];
      int j;
      for(j = 0; j < count; j++)
	{
	  int row = i;
	  int col = colind[j];
	  ergo_real value = values[j];
	  result[row*n+col] = value;
	  result[col*n+row] = value;
	}
    } // END FOR i
  
  // Remember to free memory allocated inside compute_T_matrix_sparse.
  for(i = 0; i < n; i++)
    {
      ergo_free(colindList[i]);
      ergo_free(valuesList[i]);
    }

  ergo_free(nvaluesList);
  ergo_free(colindList);
  ergo_free(valuesList);

  return 0;
}



