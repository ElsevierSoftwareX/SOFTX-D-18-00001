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

/** @file integrals_2el_J_mm_kernel.cc

    \brief Code for multipole method computational kernel for
    computing the Coulomb matrix J.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include "integrals_2el_J_mm_kernel.h"

int
do_multipole_interaction_between_2_boxes_branches(const distr_list_description_struct & distrDescription_1,
						  const multipole_struct_large & branchMultipole,
						  const multipole_struct_small* multipoleList_1,
						  ergo_real* result_J_list, // NULL if not used
						  ResultMatContrib* resultMatContrib, // NULL if not used
						  ergo_real threshold,
						  int* largest_L_used_so_far, // optional output, NULL if not used
						  MMInteractor & interactor,
						  const MMLimitTable & mmLimitTable)
{
  const batch_struct* batchList_1                 = &distrDescription_1.org.batchList[0];
  const cluster_struct* clusterList_1             = &distrDescription_1.org.clusterList[0];
  const distr_group_struct* groupList_1           = &distrDescription_1.org.groupList[0];
  const minimal_distr_struct* minimalDistrList_1  = &distrDescription_1.org.minimalDistrList[0];
  int batchCount_1                                = distrDescription_1.org.batchList.size();
  const basis_func_pair_struct* basisFuncPairList = &distrDescription_1.org.basisFuncPairList[0];
  // Prepare local result list, if needed
  std::vector<ergo_real> result_J_list_local;
  if(resultMatContrib) {
    int result_J_list_local_size = distrDescription_1.org.basisFuncPairList.size();
    result_J_list_local.resize(result_J_list_local_size);
    memset(&result_J_list_local[0], 0x00, result_J_list_local_size*sizeof(ergo_real));
  }
  int distrCountTot = 0;
  for(int batchIndex_1 = 0; batchIndex_1 < batchCount_1; batchIndex_1++) {
    int clusterCount_1 = batchList_1[batchIndex_1].noOfClusters;
    int cluster_start_1 = batchList_1[batchIndex_1].clusterStartIndex;
    for(int clusterIndex_1 = cluster_start_1; clusterIndex_1 < cluster_start_1 + clusterCount_1; clusterIndex_1++) {
      int group_start_1 = clusterList_1[clusterIndex_1].groupStartIndex;
      int group_end_1 = group_start_1 + clusterList_1[clusterIndex_1].noOfGroups;
      for(int groupIndex_1 = group_start_1; groupIndex_1 < group_end_1; groupIndex_1++) {
	const distr_group_struct* currGroup_1 = &groupList_1[groupIndex_1];
	const multipole_struct_small & multipoleCurrGroup = distrDescription_1.org_mm.multipoleListForGroups[groupIndex_1];
	ergo_real dx = branchMultipole.centerCoords[0] - currGroup_1->centerCoords[0];
	ergo_real dy = branchMultipole.centerCoords[1] - currGroup_1->centerCoords[1];
	ergo_real dz = branchMultipole.centerCoords[2] - currGroup_1->centerCoords[2];
	ergo_real r = template_blas_sqrt(dx*dx + dy*dy + dz*dz);

	// loop over distrs of 1 (and at the same time over multipoles for those distrs)
	// in order to find largest norm for each subvector.
	int distr_start = currGroup_1->startIndex;
	int distr_end = distr_start + currGroup_1->distrCount;
	ergo_real maxMomentVectorNormForDistrsListCurrGroup[MAX_MULTIPOLE_DEGREE_BASIC+1];
	for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
	  maxMomentVectorNormForDistrsListCurrGroup[l] = 0;
	int maxDegreeForDistrs = 0;
	for(int distrIndex = distr_start; distrIndex < distr_end; distrIndex++) {
	  const multipole_struct_small* distrMultipole = &multipoleList_1[distrCountTot + distrIndex - distr_start];
	  if(distrMultipole->degree > maxDegreeForDistrs)
	    maxDegreeForDistrs = distrMultipole->degree;
	  for(int l = 0; l <= distrMultipole->degree; l++) {
	    int startIndex = l*l;
	    int endIndex = (l+1)*(l+1);
	    ergo_real sum = 0;
	    for(int A = startIndex; A < endIndex; A++)
	      sum += distrMultipole->momentList[A]*distrMultipole->momentList[A];
	    ergo_real subNorm = template_blas_sqrt(sum);
	    if(subNorm > maxMomentVectorNormForDistrsListCurrGroup[l])
	      maxMomentVectorNormForDistrsListCurrGroup[l] = subNorm;
	  }
	}

	// check which degree is needed
	int degreeNeeded = mmLimitTable.get_minimum_multipole_degree_needed(r, &branchMultipole, maxDegreeForDistrs,
									    maxMomentVectorNormForDistrsListCurrGroup, threshold);
	if(degreeNeeded < 0)
	  return -1;
	if(largest_L_used_so_far != NULL) {
	  if(degreeNeeded > *largest_L_used_so_far)
	    *largest_L_used_so_far = degreeNeeded;
	}
	int branchNoOfMoments = (degreeNeeded+1)*(degreeNeeded+1);

	// create interaction matrix
	ergo_real T[multipoleCurrGroup.noOfMoments * branchNoOfMoments];
	interactor.getInteractionMatrix(dx, dy, dz, multipoleCurrGroup.degree, degreeNeeded, T);
	ergo_real tempVector[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
	for(int A = 0; A < multipoleCurrGroup.noOfMoments; A++) {
	  ergo_real sum = 0;
	  for(int B = 0; B < branchNoOfMoments; B++)
	    sum += branchMultipole.momentList[B] * T[A*branchNoOfMoments+B];
	  tempVector[A] = sum;
	}

	// loop over distrs of 1 (and at the same time over multipoles for those distrs)
	for(int distrIndex = distr_start; distrIndex < distr_end; distrIndex++) {
	  const multipole_struct_small* distrMultipole = &multipoleList_1[distrCountTot];
	  distrCountTot++;
	  ergo_real sum = 0;
	  for(int A = 0; A < distrMultipole->noOfMoments; A++)
	    sum += tempVector[A] * distrMultipole->momentList[A];
	  int basisFuncPairIndex = minimalDistrList_1[distrIndex].basisFuncPairIndex;
	  int i1 = batchList_1[batchIndex_1].basisFuncPairListIndex+basisFuncPairIndex;
	  int pairIndex = basisFuncPairList[i1].pairIndex;
	  if(result_J_list)
	    result_J_list[pairIndex] += sum;
	  else
	    result_J_list_local[i1] += sum;
	} // END FOR distrIndex
      }
    }
  }
  if(result_J_list == NULL) {
    // Transfer results from result_J_list_local to resultMatContrib
    assert(resultMatContrib != NULL);
    for(int batchIndex_1 = 0; batchIndex_1 < batchCount_1; batchIndex_1++) {
      int noOfBasisFuncPairs = batchList_1[batchIndex_1].noOfBasisFuncPairs;
      for(int i = 0; i < noOfBasisFuncPairs; i++) {
	int k = batchList_1[batchIndex_1].basisFuncPairListIndex+i;
	int a = basisFuncPairList[k].index_1;
	int b = basisFuncPairList[k].index_2;
	ergo_real currContrib = result_J_list_local[k];
	if(currContrib != 0)
	  resultMatContrib->addContrib(a, b, currContrib);
      }
    }
  }
  return 0;
}
