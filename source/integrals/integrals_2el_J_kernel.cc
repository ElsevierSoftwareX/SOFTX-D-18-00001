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

/** @file integrals_2el_J_kernel.cc

    \brief Code for computational kernel for computing the Coulomb
    matrix J.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include "integrals_2el_J_kernel.h"
#include "pi.h"
#include "integrals_hermite.h"
#include "integrals_2el_util_funcs.h"

static void transfer_to_result_J_list(int nBatchs,
				      const batch_struct* batchList,
				      const basis_func_pair_struct* basisFuncPairList,
				      const ergo_real* result_J_list_local,
				      ergo_real* result_J_list) {
  for(int batch_i = 0; batch_i < nBatchs; batch_i++) {
    int noOfBasisFuncPairs = batchList[batch_i].noOfBasisFuncPairs;
    for(int i = 0; i < noOfBasisFuncPairs; i++) {
      int k = batchList[batch_i].basisFuncPairListIndex+i;
      int J_list_index = basisFuncPairList[k].pairIndex;
      result_J_list[J_list_index] += result_J_list_local[k];
    }
  }
}

static void transfer_to_resultMatContrib(int nBatchs,
					 const batch_struct* batchList,
					 const basis_func_pair_struct* basisFuncPairList,
					 const ergo_real* result_J_list_local,
					 ResultMatContrib* resultMatContrib) {
  for(int batch_i = 0; batch_i < nBatchs; batch_i++) {
    int noOfBasisFuncPairs = batchList[batch_i].noOfBasisFuncPairs;
    for(int i = 0; i < noOfBasisFuncPairs; i++) {
      int k = batchList[batch_i].basisFuncPairListIndex+i;
      int a = basisFuncPairList[k].index_1;
      int b = basisFuncPairList[k].index_2;
      resultMatContrib->addContrib(a, b, result_J_list_local[k]);
    }
  }
}

int 
get_J_contribs_from_2_interacting_boxes(const IntegralInfo & integralInfo,
					ergo_real* result_J_list, // NULL if not used
					ResultMatContrib* resultMatContrib, // NULL if not used
					const distr_org_struct & distr_org_struct_1,
					const distr_org_struct & distr_org_struct_2,
					int interactionWithSelf,
					ergo_real threshold,
					JK_contribs_buffer_struct* bufferStructPtr)
{
  const JK::ExchWeights CAM_params_not_used;

  const ergo_real twoTimesPiToPow5half = 2 * pitopow52;
  ergo_real* summedIntegralList = bufferStructPtr->summedIntegralList;
  ergo_real* primitiveIntegralList = bufferStructPtr->primitiveIntegralList;

  const distr_group_struct* groupList_1 = &distr_org_struct_1.groupList[0];
  const distr_group_struct* groupList_2 = &distr_org_struct_2.groupList[0];
  const cluster_struct* clusterList_1 = &distr_org_struct_1.clusterList[0];
  const cluster_struct* clusterList_2 = &distr_org_struct_2.clusterList[0];
  const batch_struct* batchList_1 = &distr_org_struct_1.batchList[0];
  int nBatchs_1 = distr_org_struct_1.batchList.size();
  const batch_struct* batchList_2 = &distr_org_struct_2.batchList[0];
  int nBatchs_2 = distr_org_struct_2.batchList.size();
  const basis_func_pair_struct* basisFuncPairList_1 = &distr_org_struct_1.basisFuncPairList[0];
  const basis_func_pair_struct* basisFuncPairList_2 = &distr_org_struct_2.basisFuncPairList[0];
  const i_j_val_struct* spMatElementList_1 = &distr_org_struct_1.spMatElementList[0];
  const int* spMatCountList_1 = &distr_org_struct_1.spMatCountList[0];
  const int* spMatIdxList_1 = &distr_org_struct_1.spMatIdxList[0];
  const i_j_val_struct* spMatElementList_2 = &distr_org_struct_2.spMatElementList[0];
  const int* spMatCountList_2 = &distr_org_struct_2.spMatCountList[0];
  const int* spMatIdxList_2 = &distr_org_struct_2.spMatIdxList[0];

  // Prepare local result lists
  int result_J_list_local_1_size = distr_org_struct_1.basisFuncPairList.size();
  std::vector<ergo_real> result_J_list_local_1(result_J_list_local_1_size);
  memset(&result_J_list_local_1[0], 0x00, result_J_list_local_1_size*sizeof(ergo_real));
  int result_J_list_local_2_size = distr_org_struct_2.basisFuncPairList.size();
  std::vector<ergo_real> result_J_list_local_2(result_J_list_local_2_size);
  memset(&result_J_list_local_2[0], 0x00, result_J_list_local_2_size*sizeof(ergo_real));

  for(int batch_i = 0; batch_i < nBatchs_1; batch_i++)
    {
      int batch_j_start = 0;
      if(interactionWithSelf == 1)
	batch_j_start = batch_i;
      for(int batch_j = batch_j_start; batch_j < nBatchs_2; batch_j++)
	{
	  int noOfBasisFuncPairs_1 = batchList_1[batch_i].noOfBasisFuncPairs;
	  int noOfBasisFuncPairs_2 = batchList_2[batch_j].noOfBasisFuncPairs;
	  // set integral list to zero
	  memset(summedIntegralList, 0, noOfBasisFuncPairs_1*noOfBasisFuncPairs_2*sizeof(ergo_real));

	  // get largest dmat element
	  ergo_real maxabsdmatelement = 0;
	  for(int i = 0; i < noOfBasisFuncPairs_1; i++)
	    for(int j = 0; j < noOfBasisFuncPairs_2; j++)
	      {
		ergo_real D_ab = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+i].dmatElement;
		ergo_real D_cd = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+j].dmatElement;
		ergo_real absval;
		absval = template_blas_fabs(D_ab);
		if(absval > maxabsdmatelement)
		  maxabsdmatelement = absval;
		absval = template_blas_fabs(D_cd);
		if(absval > maxabsdmatelement)
		  maxabsdmatelement = absval;
	      } // END FOR i j get largest dmat element
	  int cluster_i_start = batchList_1[batch_i].clusterStartIndex;
	  int clusterCount1 = batchList_1[batch_i].noOfClusters;
	  for(int cluster_i = cluster_i_start; cluster_i < cluster_i_start + clusterCount1; cluster_i++)
	    {
	      int cluster_j_start = batchList_2[batch_j].clusterStartIndex;
	      int clusterCount2 = batchList_2[batch_j].noOfClusters;
	      int cluterIndexEnd2 = cluster_j_start + clusterCount2;
	      if(interactionWithSelf == 1 && batch_i == batch_j)
		cluster_j_start = cluster_i;
	      for(int cluster_j = cluster_j_start; cluster_j < cluterIndexEnd2; cluster_j++)
		{
		  // check if we can skip this combination of clusters
		  if(clusterList_1[cluster_i].maxLimitingFactorForCluster * clusterList_2[cluster_j].maxLimitingFactorForCluster * maxabsdmatelement < threshold)
		    continue;

		  int group_i_start = clusterList_1[cluster_i].groupStartIndex;
		  int group_i_end = group_i_start + clusterList_1[cluster_i].noOfGroups;
		  int group_j_start = clusterList_2[cluster_j].groupStartIndex;
		  int group_j_end = group_j_start + clusterList_2[cluster_j].noOfGroups;

		  int n1max = clusterList_1[cluster_i].nmax;
		  int n2max = clusterList_2[cluster_j].nmax;

		  // Now we can precompute things that depend only on exponents
		  ergo_real alpha_1 = groupList_1[group_i_start].exponent;
		  ergo_real alpha_2 = groupList_2[group_j_start].exponent;
		  ergo_real alphasum = alpha_1 + alpha_2;
		  ergo_real alphaproduct = alpha_1 * alpha_2;
		  ergo_real alpha_0 = alphaproduct / alphasum;

		  ergo_real resultPreFactor = twoTimesPiToPow5half / (alphaproduct*template_blas_sqrt(alphasum));

		  for(int group_i = group_i_start; group_i < group_i_end; group_i++)
		    {
		      if(interactionWithSelf == 1 && batch_i == batch_j && cluster_i == cluster_j)
			group_j_start = group_i;
		      for(int group_j = group_j_start; group_j < group_j_end; group_j++)
			{
			  // Only J is considered; we can use maxAbsDmatElementGroup
			  ergo_real maxabs_1 = groupList_1[group_i].maxAbsDmatElementGroup;
			  ergo_real maxabs_2 = groupList_2[group_j].maxAbsDmatElementGroup;
			  if((groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabs_1 < threshold) && 
			     (groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabs_2 < threshold))
			    continue;

			  // now we can do all integrals needed for this pair of groups
			  ergo_real dx = groupList_2[group_j].centerCoords[0] - groupList_1[group_i].centerCoords[0];
			  ergo_real dy = groupList_2[group_j].centerCoords[1] - groupList_1[group_i].centerCoords[1];
			  ergo_real dz = groupList_2[group_j].centerCoords[2] - groupList_1[group_i].centerCoords[2];

			  // now we have dx dy dz alpha0 alpha1 n1max n2max. Get all integrals for this case.
			  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1max];
			  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2max];

			  get_related_integrals_hermite(integralInfo,
							CAM_params_not_used,
							n1max, noOfMonomials_1,
							n2max, noOfMonomials_2,
							dx, dy, dz, alpha_0,
							resultPreFactor,
							primitiveIntegralList);

			  if(interactionWithSelf == 1 && group_j == group_i && batch_i == batch_j && cluster_i == cluster_j) {
			    do_summedIntegralList_contribs_self(&spMatElementList_1[spMatIdxList_1[group_i]], spMatCountList_1[group_i],
								&spMatElementList_2[spMatIdxList_2[group_j]], spMatCountList_2[group_j],
								noOfMonomials_1, noOfMonomials_2,
								primitiveIntegralList,
								noOfBasisFuncPairs_1, noOfBasisFuncPairs_2,
								summedIntegralList);
			  }
			  else {
			    do_summedIntegralList_contribs_std(&spMatElementList_1[spMatIdxList_1[group_i]], spMatCountList_1[group_i],
							       &spMatElementList_2[spMatIdxList_2[group_j]], spMatCountList_2[group_j],
							       noOfMonomials_1, noOfMonomials_2,
							       primitiveIntegralList,
							       noOfBasisFuncPairs_1, noOfBasisFuncPairs_2,
							       summedIntegralList);
			  }

			} // END FOR group_j
		    } // END FOR group_i
		} // END FOR cluster_j
	    } // END FOR cluster_i

	  for(int idx_1 = 0; idx_1 < noOfBasisFuncPairs_1; idx_1++)
	    for(int idx_2 = 0; idx_2 < noOfBasisFuncPairs_2; idx_2++)
	      {
		int a = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_1;
		int b = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_2;
		int c = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_1;
		int d = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_2;
		ergo_real integralValueCurr = summedIntegralList[idx_1*noOfBasisFuncPairs_2 + idx_2];

		ergo_real D_ab = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].dmatElement;
		ergo_real D_cd = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].dmatElement;

		int i1 = batchList_1[batch_i].basisFuncPairListIndex+idx_1;
		int J_list_index_ab = basisFuncPairList_1[i1].pairIndex;
		int i2 = batchList_2[batch_j].basisFuncPairListIndex+idx_2;
		int J_list_index_cd = basisFuncPairList_2[i2].pairIndex;

		// Multiply integralValueCurr by 2 if ab and cd refer to the same index pair.
		// This is done differently depending on whether J_list_index_ab info is available (it may be -1)
		// FIXME ELIAS: put this factor of 2 in the if-statements below instead, that should be possible and avoid this complication wuth (J_list_index_ab >= 0) here.
		// (The places below where the factor 2 could be added is those places where there is only one contribution to result_J_list instead of two.)
		if(J_list_index_ab >= 0) {
		  if(J_list_index_ab == J_list_index_cd)
		    integralValueCurr *= 2;
		}
		else {
		  if((a == c && b == d) || (a == d && b == c))
		    integralValueCurr *= 2;
		}

		if(template_blas_fabs(integralValueCurr)*maxabsdmatelement < threshold)
		  continue;

		// Place results in result_J_list_local_1 and result_J_list_local_2
		if(a != b && c != d && (a != c || b != d)) {
		  result_J_list_local_1[i1] += 2 * D_cd * integralValueCurr;
		  result_J_list_local_2[i2] += 2 * D_ab * integralValueCurr;
		}
		else if(a != b && c != d && a == c && b == d) {
		  result_J_list_local_1[i1] += 2 * D_cd * integralValueCurr;
		}
		else if(a == b && c != d) {
		  result_J_list_local_1[i1] += 2 * D_cd * integralValueCurr;
		  result_J_list_local_2[i2] += 1 * D_ab * integralValueCurr;
		}
		else if(a != b && c == d) {
		  result_J_list_local_1[i1] += 1 * D_cd * integralValueCurr;
		  result_J_list_local_2[i2] += 2 * D_ab * integralValueCurr;
		}
		else if(a == b && c == d && a != c) {
		  result_J_list_local_1[i1] += D_cd * integralValueCurr;
		  result_J_list_local_2[i2] += D_ab * integralValueCurr;
		}
		else if(a == b && c == d && a == c) {
		  result_J_list_local_1[i1] += D_cd * integralValueCurr;
		}
		else {
		  return -1; // This should never happen
		}

	      } // END FOR idx_1 idx_2
	} // END FOR batch_j
    } // END FOR batch_i

  // Transfer results from local result lists to final result location.
  if(result_J_list) {
    // Place results in result_J_list
    transfer_to_result_J_list(nBatchs_1, batchList_1, basisFuncPairList_1, &result_J_list_local_1[0], result_J_list);
    transfer_to_result_J_list(nBatchs_2, batchList_2, basisFuncPairList_2, &result_J_list_local_2[0], result_J_list);
  }
  else {
    // Place results in resultMatContrib
    assert(resultMatContrib != NULL);
    transfer_to_resultMatContrib(nBatchs_1, batchList_1, basisFuncPairList_1, &result_J_list_local_1[0], resultMatContrib);
    transfer_to_resultMatContrib(nBatchs_2, batchList_2, basisFuncPairList_2, &result_J_list_local_2[0], resultMatContrib);
  }

  return 0;
}

