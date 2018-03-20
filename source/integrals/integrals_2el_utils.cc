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

/** @file integrals_2el_utils.cc

    \brief Code for various utilities used by 2-electron integral
    computation (i.e. computation of J and K matrices).

    @author: Elias Rudberg <em>responsible</em>.
*/

#include <string.h>
#include <stdexcept>

#include "integrals_2el_utils.h"
#include "integrals_hermite.h"
#include "template_blas_common.h"
#include "basis_func_extent.h"
#include "integrals_2el_repeating.h"
#include "integrals_general.h"


/* ----- ResultMatContrib implementation ----- */

ResultMatContrib::ResultMatContrib() {
  currVecIndex = 0;
  currContribCount = 0;
  indexInCurrVec = 0;
  currVecReservedSize = 20;
  for(int i = 0; i < nVectorsMax; i++)
    vList[i] = 0;
  vList[currVecIndex] = new std::vector<RowColVal>;
  vList[currVecIndex]->reserve(currVecReservedSize);
}

ResultMatContrib::~ResultMatContrib() {
  for(int i = 0; i <= currVecIndex; i++)
    delete vList[i];
}

void ResultMatContrib::addContrib(int row, int col, ergo_real value) {
  // Check if there is room in current vector
  assert(indexInCurrVec <= currVecReservedSize);
  if(indexInCurrVec == currVecReservedSize) {
    // Current vector is full. Reserve space in next vector and use it.
    currVecIndex++;
    if(currVecIndex >= nVectorsMax)
      throw std::runtime_error("Error in ResultMatContrib::addContrib: (currVecIndex >= nVectorsMax).");
    indexInCurrVec = 0;
    currVecReservedSize *= 2;
    vList[currVecIndex] = new std::vector<RowColVal>;
    vList[currVecIndex]->reserve(currVecReservedSize);
  }
  RowColVal tmp = {row, col, value};
  vList[currVecIndex]->push_back(tmp);
  indexInCurrVec++;
  currContribCount++;
}

const ResultMatContrib::RowColVal & ResultMatContrib::fetchNextContrib(int & currVecIndexForFetch, int & indexInCurrVecForFetch) const {
  assert(currVecIndexForFetch >= 0 && currVecIndexForFetch < nVectorsMax);
  int currVecSize = vList[currVecIndexForFetch]->size();
  assert(indexInCurrVecForFetch >= 0 && indexInCurrVecForFetch < currVecSize);
  const RowColVal & result = (*vList[currVecIndexForFetch])[indexInCurrVecForFetch];
  indexInCurrVecForFetch++;
  if(indexInCurrVecForFetch == currVecSize) {
    currVecIndexForFetch++;
    indexInCurrVecForFetch = 0;
  }
  return result;
}

/* ----- End of ResultMatContrib implementation ----- */


ergo_real
get_max_abs_vector_element(int n, const ergo_real* vector)
{
  ergo_real maxabs = 0;
  for(int i = 0; i < n; i++)
    {
      ergo_real absval = template_blas_fabs(vector[i]);
      if(absval > maxabs)
	maxabs = absval;
    }
  return maxabs;
}




box_struct::box_struct()
{
  memset(branchIndexListForJ, 0, MAX_NO_OF_BRANCHES*sizeof(int));
  memset(branchCountListForJ, 0, MAX_NO_OF_BRANCHES*sizeof(int));
}


void
allocate_buffers_needed_by_integral_code(const IntegralInfo & integralInfo, 
					 int maxNoOfMonomials,
					 int basisFuncListCount_max,
					 JK_contribs_buffer_struct* bufferStruct)
{
  bufferStruct->summedIntegralList = new ergo_real[MAX_NO_OF_BASIS_FUNC_PAIRS_PER_BATCH*MAX_NO_OF_BASIS_FUNC_PAIRS_PER_BATCH];
  bufferStruct->primitiveIntegralList = new ergo_real[maxNoOfMonomials*maxNoOfMonomials];
  bufferStruct->primitiveIntegralList_work = new ergo_real[maxNoOfMonomials*maxNoOfMonomials];
  if(basisFuncListCount_max > 0)
    {
      bufferStruct->partial_dmat_1 = new ergo_real[basisFuncListCount_max*basisFuncListCount_max];
      bufferStruct->partial_K_1    = new ergo_real[basisFuncListCount_max*basisFuncListCount_max];
      // FIXME: only allocate _2 buffers if nonsymm case.
      bufferStruct->partial_dmat_2 = new ergo_real[basisFuncListCount_max*basisFuncListCount_max];
      bufferStruct->partial_K_2    = new ergo_real[basisFuncListCount_max*basisFuncListCount_max];
    }
  else
    {
      bufferStruct->partial_dmat_1 = NULL;
      bufferStruct->partial_K_1    = NULL;
      bufferStruct->partial_dmat_2 = NULL;
      bufferStruct->partial_K_2    = NULL;
    }
}

void
free_buffers_needed_by_integral_code(JK_contribs_buffer_struct* bufferStruct)
{
  delete [] bufferStruct->summedIntegralList;
  delete [] bufferStruct->primitiveIntegralList;
  delete [] bufferStruct->primitiveIntegralList_work;
  if(bufferStruct->partial_dmat_1)
    delete [] bufferStruct->partial_dmat_1;
  if(bufferStruct->partial_K_1)
    delete [] bufferStruct->partial_K_1;
  if(bufferStruct->partial_dmat_2)
    delete [] bufferStruct->partial_dmat_2;
  if(bufferStruct->partial_K_2)
    delete [] bufferStruct->partial_K_2;
  bufferStruct->summedIntegralList = NULL;
  bufferStruct->primitiveIntegralList = NULL;
  bufferStruct->primitiveIntegralList_work = NULL;
  bufferStruct->partial_dmat_1 = NULL;
  bufferStruct->partial_K_1 = NULL;
  bufferStruct->partial_dmat_2 = NULL;
  bufferStruct->partial_K_2 = NULL;
}


int
get_related_integrals_h(
			const IntegralInfo & integralInfo,
			const JK::ExchWeights & CAM_params,
			int n1max, int noOfMonomials_1,
			int n2max, int noOfMonomials_2,
			ergo_real dx0, 
			ergo_real dx1, 
			ergo_real dx2, 
			ergo_real alpha1, 
			ergo_real alpha2,
			ergo_real alpha0,
			ergo_real* primitiveIntegralList,
			ergo_real* primitiveIntegralList_work,
			ergo_real resultPreFactor
			)
{
  get_related_integrals_hermite(integralInfo,
				CAM_params,
				n1max, noOfMonomials_1,
				n2max, noOfMonomials_2,
				dx0, 
				dx1, 
				dx2, 
				alpha0,
				resultPreFactor,
				primitiveIntegralList);

  integralInfo.multiply_by_hermite_conversion_matrix_from_right(n1max,
								n2max,
								1.0/alpha1,
								primitiveIntegralList,
								primitiveIntegralList_work);
  integralInfo.multiply_by_hermite_conversion_matrix_from_left(n1max,
							       n2max,
							       1.0/alpha2,
							       primitiveIntegralList_work,
							       primitiveIntegralList);
  
  return 0;
}



static ergo_real
erfc_inverse(ergo_real x, ergo_real requested_accuracy)
{
  ergo_real y_min = 0.0;
  ergo_real y_max = 10.0;
  if(template_blas_erfc(y_max) > x)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in erfc_inverse: (erfc(y_max) > x)");
      exit(EXIT_FAILURE);
    }
  int count = 0;
  ergo_real y = 0;
  while(y_max - y_min > requested_accuracy)
    {
      y = (y_min + y_max) / 2;
      if(template_blas_erfc(y) > x)
	y_min = y;
      else
	y_max = y;
      count++;
      if(count > 222)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in erfc_inverse: too many iterations.");
	  exit(EXIT_FAILURE);
	}
    } // END WHILE requested accuracy not reached
  return y;
}




void
compute_extent_for_list_of_distributions(int n,
					 DistributionSpecStructLabeled* distrList,
					 ergo_real threshold,
					 ergo_real maxLimitingFactor,
					 ergo_real maxabsDmatelement)
{
  ergo_real requested_erfcinv_accuracy = 0.001;
  for(int i = 0; i < n; i++)
    {
      ergo_real erfc_inverse_value = erfc_inverse(threshold / (distrList[i].limitingFactor * maxLimitingFactor * maxabsDmatelement), requested_erfcinv_accuracy);
      distrList[i].distr.extent = erfc_inverse_value / template_blas_sqrt(distrList[i].distr.exponent);
    }
}



int
get_list_of_labeled_distrs_maxLimitingFactor(const BasisInfoStruct & basisInfo,
					     const IntegralInfo & integralInfo,
					     ergo_real threshold,
					     ergo_real* resultMaxLimitingFactor,
					     ergo_real maxDensityMatrixElement)
{
  int n = basisInfo.noOfBasisFuncs;

  std::vector<ergo_real> basisFuncExtentList(n);
  if(compute_extent_for_all_basis_funcs_2el(integralInfo,
					    basisInfo, 
					    &basisFuncExtentList[0], 
					    threshold,
					    maxDensityMatrixElement) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_extent_for_all_basis_funcs");
      return -1;
    }

  ergo_real maxExtent = 0;
  for(int i = 0; i < n; i++)
    {
      ergo_real currExtent = basisFuncExtentList[i];
      if(currExtent > maxExtent)
	maxExtent = currExtent;
    }
  std::vector<int> orgIndexList(n);

  // Create box system for basisInfo.
  std::vector<box_item_struct> itemList(n);
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < 3; j++)
	itemList[i].centerCoords[j] = basisInfo.basisFuncList[i].centerCoords[j];
      itemList[i].originalIndex = i;
    }
  ergo_real toplevelBoxSize = 7.0;
  BoxSystem boxSystem;
  if(boxSystem.create_box_system(&itemList[0],
				 n,
				 toplevelBoxSize) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system.");
      return -1;
    }

  IntegratorWithMemory integrator(&integralInfo);

  ergo_real maxLimitingFactor = 0;
  for(int i = 0; i < n; i++)
    {
      // Now, instead of looping again over all n basis functions, we use box system to find relevant ones.
      ergo_real maxDistance = basisFuncExtentList[i] + maxExtent;
      ergo_real coords[3];
      for(int coordNo = 0; coordNo < 3; coordNo++)
	coords[coordNo] = basisInfo.basisFuncList[i].centerCoords[coordNo];
      int nRelevant = boxSystem.get_items_near_point(&itemList[0], coords, maxDistance, &orgIndexList[0]);
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
	  
	  const int maxCountProduct = 10000;
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
	      ergo_real limitingFactor = template_blas_sqrt(integrator.do_2e_integral(&psi_list[k]));
	      if(limitingFactor > maxLimitingFactor)
		maxLimitingFactor = limitingFactor;
	    } // END FOR k
	} // END FOR j
    } // END FOR i
  *resultMaxLimitingFactor = maxLimitingFactor;

  return 0;
}



int
get_list_of_labeled_distrs(const BasisInfoStruct & basisInfo,
			   const IntegralInfo & integralInfo,
			   ergo_real threshold,
			   DistributionSpecStructLabeled* resultList,
			   int maxCountDistrs,
			   ergo_real maxLimitingFactor,
			   const ergo_real* dens,
			   ergo_real maxDensityMatrixElement)
{
  int n = basisInfo.noOfBasisFuncs;
  
  // compute extent for all basis functions
  std::vector<ergo_real> basisFuncExtentList(n);
  if(compute_extent_for_all_basis_funcs_2el(integralInfo,
					    basisInfo, 
					    &basisFuncExtentList[0], 
					    threshold,
					    maxDensityMatrixElement) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_extent_for_all_basis_funcs");
      return -1;
    }

  ergo_real maxExtent = 0;
  for(int i = 0; i < n; i++)
    {
      ergo_real currExtent = basisFuncExtentList[i];
      if(currExtent > maxExtent)
	maxExtent = currExtent;
    }
  std::vector<int> orgIndexList(n);

  // Create box system for basisInfo.
  std::vector<box_item_struct> itemList(n);
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < 3; j++)
	itemList[i].centerCoords[j] = basisInfo.basisFuncList[i].centerCoords[j];
      itemList[i].originalIndex = i;
    }
  ergo_real toplevelBoxSize = 7.0;
  BoxSystem boxSystem;
  if(boxSystem.create_box_system(&itemList[0],
				 n,
				 toplevelBoxSize) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system.");
      return -1;
    }

  IntegratorWithMemory integrator(&integralInfo);

  // create list of product primitives, with labels
  int distrCount = 0;
  for(int i = 0; i < n; i++)
    {
      // Now, instead of looping again over all n basis functions, we use box system to find relevant ones.
      ergo_real maxDistance = basisFuncExtentList[i] + maxExtent;
      ergo_real coords[3];
      for(int coordNo = 0; coordNo < 3; coordNo++)
	coords[coordNo] = basisInfo.basisFuncList[i].centerCoords[coordNo];
      int nRelevant = boxSystem.get_items_near_point(&itemList[0], coords, maxDistance, &orgIndexList[0]);
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
	  
	  // Set dmatElement if dens given, otherwise just set it to zero.
	  ergo_real dmatElement = 0;
	  if(dens != NULL)
	    dmatElement = dens[i*n+j];

	  const int maxCountProduct = 10000;
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
	      ergo_real limitingFactor = template_blas_sqrt(integrator.do_2e_integral(&psi_list[k]));
	      if(limitingFactor*maxLimitingFactor*maxDensityMatrixElement > threshold)
		{
		  if(maxCountDistrs > 0 && distrCount >= maxCountDistrs)
		    {
		      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_list_of_labeled_distrs: (maxCountDistrs > 0 && distrCount >= maxCountDistrs)");
		      return -1;
		    }
		  if(resultList != NULL)
		    {
		      resultList[distrCount].distr = psi_list[k];
		      resultList[distrCount].basisFuncIndex_1 = i;
		      resultList[distrCount].basisFuncIndex_2 = j;
		      resultList[distrCount].pairIndex = -1; // not used
		      resultList[distrCount].limitingFactor = limitingFactor;
		      resultList[distrCount].dmatElement = dmatElement;
		    }
		  distrCount++;
		} // END IF above threshold
	    } // END FOR k
	} // END FOR j
    } // END FOR i

  return distrCount;
}



static void
create_item_list_from_list_of_distributions(int n, 
					    const DistributionSpecStructLabeled* distrList, 
					    box_item_struct* itemList)
{
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < 3; j++)
	itemList[i].centerCoords[j] = distrList[i].distr.centerCoords[j];      
      itemList[i].originalIndex = i;
    } // END FOR i
}


int
create_box_system_and_reorder_distrs(int distrCount,
				     DistributionSpecStructLabeled* distrList,
				     ergo_real toplevelBoxSize,
				     BoxSystem & boxSystem)
{
  std::vector<box_item_struct> itemList(distrCount);
  create_item_list_from_list_of_distributions(distrCount, &distrList[0], &itemList[0]);
  if(boxSystem.create_box_system(&itemList[0],
				 distrCount,
				 toplevelBoxSize) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in create_box_system");
    return -1;
  }
  // reorder list of labeled distrs, where they are ordered box by box
  // at the level of smallest boxes.  Since the distr struct is rather
  // big we want to do this reordering without using a whole new
  // vector of structs. Instead we use two vectors of int, one
  // containing the current location of each distr and the other
  // containing the contents of each entry.
  std::vector<int> locations(distrCount);
  for(int i = 0; i < distrCount; i++)
    locations[i] = i;
  std::vector<int> contents(distrCount);
  for(int i = 0; i < distrCount; i++)
    contents[i] = i;
  for(int i = 0; i < distrCount; i++) {
    // Now we want to copy the correct entry into position i. We also
    // need to copy the entry that was occupying this spot to
    // somewhere else, and update the locations and contents vectors
    // accordingly.
    int otherIdx = locations[itemList[i].originalIndex];
    // Switch place of distr structs.
    DistributionSpecStructLabeled temp = distrList[otherIdx];
    distrList[otherIdx] = distrList[i];
    distrList[i] = temp;
    // Update locations and contents vectors.
    int orgIdxToMove = contents[i];
    contents[i] = itemList[i].originalIndex;
    locations[itemList[i].originalIndex] = i;
    contents[otherIdx] = orgIdxToMove;
    locations[orgIdxToMove] = otherIdx;
  }
  return 0;
}
