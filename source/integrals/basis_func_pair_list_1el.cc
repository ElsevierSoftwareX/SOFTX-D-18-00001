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

/** @file basis_func_pair_list_1el.cc

    @brief Functions for setting up lists of non-negligible basis
    function pairs, for 1-electron integrals.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "basis_func_pair_list_1el.h"
#include "basis_func_extent_1el.h"
#include "output.h"
#include "pi.h"
#include "box_system.h"
#include "utilities.h"
#include "integrals_general.h"

int
get_basis_func_pair_list_1el_for_V(const BasisInfoStruct & basisInfo,
				   ergo_real threshold,
				   ergo_real boxSize,
				   ergo_real maxCharge,
				   basis_func_index_pair_struct_1el* result_basisFuncPairList,
				   int resultMaxCount)
{
  int n = basisInfo.noOfBasisFuncs;
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "get_basis_func_pair_list_1el_for_V, n = %6i", n);
  
  Util::TimeMeter timeMeter;
  
  // compute extent for all basis functions
  ergo_real* basisFuncExtentList = new ergo_real[n];
  if(compute_extent_for_all_basis_funcs_1el(basisInfo, 
					    basisFuncExtentList, 
					    maxCharge,
					    threshold) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_extent_for_all_basis_funcs_1el");
      return -1;
    }
  ergo_real maxExtent = 0;
  for(int i = 0; i < n; i++)
    {
      ergo_real currExtent = basisFuncExtentList[i];
      if(currExtent > maxExtent)
	maxExtent = currExtent;
    }

  // Create box system for basisInfo.
  box_item_struct* itemList = new box_item_struct[n];
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < 3; j++)
	itemList[i].centerCoords[j] = basisInfo.basisFuncList[i].centerCoords[j];
      itemList[i].originalIndex = i;
    }
  ergo_real toplevelBoxSize = boxSize;
  BoxSystem boxSystem;
  if(boxSystem.create_box_system(itemList,
				 n,
				 toplevelBoxSize) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system.");
      return -1;
    }

  int* orgIndexList = new int [n];
  
  int count = 0;
  for(int i = 0; i < n; i++)
    {
      // Now, instead of looping again over all n basis functions, we use box system to find relevant ones.
      ergo_real maxDistance = basisFuncExtentList[i] + maxExtent;
      ergo_real coords[3];
      for(int coordNo = 0; coordNo < 3; coordNo++)
	coords[coordNo] = basisInfo.basisFuncList[i].centerCoords[coordNo];
      int nRelevant = boxSystem.get_items_near_point(itemList, coords, maxDistance, orgIndexList);
      for(int jRelevant = 0; jRelevant < nRelevant; jRelevant++)
	{
	  int j = orgIndexList[jRelevant];
	  if(j < i)
	    continue;
	  // Now we are concerned with basis functions i and j.
	  // If they are far enough apart, we can skip this pair.
	  ergo_real dx = basisInfo.basisFuncList[i].centerCoords[0] - basisInfo.basisFuncList[j].centerCoords[0];
	  ergo_real dy = basisInfo.basisFuncList[i].centerCoords[1] - basisInfo.basisFuncList[j].centerCoords[1];
	  ergo_real dz = basisInfo.basisFuncList[i].centerCoords[2] - basisInfo.basisFuncList[j].centerCoords[2];
	  ergo_real distance = template_blas_sqrt(dx*dx + dy*dy + dz*dz);
	  if(distance > basisFuncExtentList[i] + basisFuncExtentList[j])
	    continue;
	  // There may be some overlap between these two basis functions.
	  // However, the extent check is rather rough.
	  // To check more carefully, compute product explicitly.
	  int currProductLargeEnough = 0;
	  const int maxCountProduct = POLY_PRODUCT_MAX_DISTRS;
	  DistributionSpecStruct psi_list[maxCountProduct];
	  /* form product of basisfuncs i and j, store product in psi_list */
	  int n_psi = get_product_simple_primitives(basisInfo, i,
						    basisInfo, j,
						    psi_list,
						    maxCountProduct,
						    0);
	  if(n_psi < 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives");
	      return -1;
	    }
	  for(int k = 0; k < n_psi; k++)
	    {
	      // Check if this distr can give contribution larger than threshold.
	      ergo_real estimatedMaxContrib = template_blas_fabs(2 * pi * maxCharge * psi_list[k].coeff / psi_list[k].exponent);
	      if(estimatedMaxContrib > threshold)
		{
		  // This product distr is large enough.
		  currProductLargeEnough = 1;
		  break;
		} // END IF above threshold
	    } // END FOR k
	  if(currProductLargeEnough == 1)
	    {
	      // Include this pair in the list
	      if(result_basisFuncPairList != NULL)
		{
		  if(count >= resultMaxCount)
		    {
		      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_basis_func_pair_list: (count >= resultMaxCount)");
		      return -1;
		    }
		  result_basisFuncPairList[count].index_1 = i;
		  result_basisFuncPairList[count].index_2 = j;
		}
	      count++;
	    } // END IF include this pair in the list
	} // END FOR jRelevant
    } // END FOR i
  
  delete [] basisFuncExtentList;
  delete [] itemList;
  delete [] orgIndexList;

  timeMeter.print(LOG_AREA_INTEGRALS, "get_basis_func_pair_list_1el_for_V");
  
  return count;
}


/* get_basis_func_pair_list_simple: Get list of non-negligible basis
   function index pairs (i,j) such that j<=i. */
int
get_basis_func_pair_list_simple(const BasisInfoStruct & basisInfo,
				ergo_real threshold,
				ergo_real boxSize,
				basis_func_index_pair_struct_1el* result_basisFuncPairList,
				int resultMaxCount)
{
  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "get_basis_func_pair_list_simple, n = %6i", n);

  Util::TimeMeter timeMeter;

  // compute extent for all basis functions
  ergo_real A = get_largest_simple_integral(basisInfo);
  std::vector<ergo_real> basisFuncExtentList(n);
  get_basis_func_extent_list(basisInfo, &basisFuncExtentList[0], threshold / A);

  ergo_real maxExtent = 0;
  for(int i = 0; i < n; i++) {
    ergo_real currExtent = basisFuncExtentList[i];
    if(currExtent > maxExtent)
      maxExtent = currExtent;
  }

  // Create box system for basisInfo.
  std::vector<box_item_struct> itemList(n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < 3; j++)
      itemList[i].centerCoords[j] = basisInfo.basisFuncList[i].centerCoords[j];
    itemList[i].originalIndex = i;
  }
  ergo_real toplevelBoxSize = boxSize;
  BoxSystem boxSystem;
  if(boxSystem.create_box_system(&itemList[0], n, toplevelBoxSize) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system.");
    return -1;
  }

  std::vector<int> orgIndexList(n);
  
  int count = 0;
  for(int i = 0; i < n; i++) {
    // Now, instead of looping again over all n basis functions, we use box system to find relevant ones.
    ergo_real maxDistance = basisFuncExtentList[i] + maxExtent;
    ergo_real coords[3];
    for(int coordNo = 0; coordNo < 3; coordNo++)
      coords[coordNo] = basisInfo.basisFuncList[i].centerCoords[coordNo];
    int nRelevant = boxSystem.get_items_near_point(&itemList[0], coords, maxDistance, &orgIndexList[0]);
    for(int jRelevant = 0; jRelevant < nRelevant; jRelevant++) {
      int j = orgIndexList[jRelevant];
      if(j > i) // Here we decide to include only j<=i cases, skip if j>i
	continue;
      // Now we are concerned with basis functions i and j.
      // If they are far enough apart, we can skip this pair.
      ergo_real dx = basisInfo.basisFuncList[i].centerCoords[0] - basisInfo.basisFuncList[j].centerCoords[0];
      ergo_real dy = basisInfo.basisFuncList[i].centerCoords[1] - basisInfo.basisFuncList[j].centerCoords[1];
      ergo_real dz = basisInfo.basisFuncList[i].centerCoords[2] - basisInfo.basisFuncList[j].centerCoords[2];
      ergo_real distance = template_blas_sqrt(dx*dx + dy*dy + dz*dz);
      if(distance > basisFuncExtentList[i] + basisFuncExtentList[j])
	continue;
      // Include this pair in the list
      if(result_basisFuncPairList != NULL) {
	if(count >= resultMaxCount) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_basis_func_pair_list_simple: (count >= resultMaxCount)");
	  return -1;
	}
	result_basisFuncPairList[count].index_1 = i;
	result_basisFuncPairList[count].index_2 = j;
      }
      count++;
    } // END FOR jRelevant
  } // END FOR i

  timeMeter.print(LOG_AREA_INTEGRALS, "get_basis_func_pair_list_simple");

  return count;
}


