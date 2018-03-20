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

/** @file integrals_2el_K_kernel.cc

    \brief Code for computational kernel for computing the
    Hartree-Fock exchange matrix K.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include "integrals_2el_K_kernel.h"
#include "pi.h"
#include "integrals_hermite.h"
#include "integrals_2el_util_funcs.h"

pthread_mutex_t K_CSR_shared_access_mutex = PTHREAD_MUTEX_INITIALIZER;

static inline int
ergo_csr_find_index_inline(const csr_matrix_struct* csr, int row, int col) {
  int n = csr->rowList[row].noOfElementsInRow;
  int baseIndex = csr->rowList[row].firstElementIndex;
  int* colList = &csr->columnIndexList[baseIndex];
  int lo = 0;
  int hi = n-1;
  while(lo < hi - 1) {
    int mid = (lo + hi) / 2;
    if(colList[mid] < col)
      lo = mid;
    else
      hi = mid;
  }
  if(colList[lo] == col)
    return baseIndex + lo;
  if(colList[hi] == col)
    return baseIndex + hi;
  // Not found
  return -1;
}

static inline ergo_real 
ergo_CSR_get_element_inline(const csr_matrix_struct* csr, int row, int col) {
  int row2 = row;
  int col2 = col;
  if(csr->symmetryFlag) {
    if(row > col) {
      row2 = col;
      col2 = row;
    }
  }
  int i = ergo_csr_find_index_inline(csr, row2, col2);
  if(i < 0)
    return 0;
  return csr->elementList[i];
}

struct abcd_struct {
  int a, b, c, d;
  int poly_ab_index;
  int poly_cd_index;
  int idx1;
  int idx2;
  ergo_real densValue;
};

#define set_abcd_list_item_macro(i,A,B,C,D,v,i1,i2)			\
  list[i].a = A; list[i].b = B; list[i].c = C; list[i].d = D; list[i].densValue = v; list[i].idx1 = i1; list[i].idx2 = i2; 

int
get_K_contribs_from_2_interacting_boxes(const IntegralInfo & integralInfo,
					const JK::ExchWeights & CAM_params,
					int maxNoOfMonomials,
					csr_matrix_struct* result_K_CSR_shared, // NULL if not used
					ResultMatContrib* resultMatContrib, // NULL if not used
					const csr_matrix_struct* dens_CSR,
					int symmetryFlag,
					const distr_org_struct & distr_org_struct_1,
					const distr_org_struct & distr_org_struct_2,
					int interactionWithSelf,
					ergo_real threshold,
					JK_contribs_buffer_struct* bufferStructPtr,
					int use_multipole_screening_for_clusters,
					ergo_real boxDistance)
{
  const ergo_real twoTimesPiToPow5half = 2 * pitopow52;// = 2 * pow(pi, 2.5);
  ergo_real* summedIntegralList = bufferStructPtr->summedIntegralList;
  ergo_real* primitiveIntegralList = bufferStructPtr->primitiveIntegralList;

  int nBatchs_1 = distr_org_struct_1.batchList.size();
  int nBatchs_2 = distr_org_struct_2.batchList.size();
  const batch_struct* batchList_1 = &distr_org_struct_1.batchList[0];
  const batch_struct* batchList_2 = &distr_org_struct_2.batchList[0];
  const cluster_struct* clusterList_1 = &distr_org_struct_1.clusterList[0];
  const cluster_struct* clusterList_2 = &distr_org_struct_2.clusterList[0];
  const distr_group_struct* groupList_1 = &distr_org_struct_1.groupList[0];
  const distr_group_struct* groupList_2 = &distr_org_struct_2.groupList[0];
  const basis_func_pair_struct* basisFuncPairList_1 = &distr_org_struct_1.basisFuncPairList[0];
  const basis_func_pair_struct* basisFuncPairList_2 = &distr_org_struct_2.basisFuncPairList[0];
  const int* basisFuncListForBatchs_map_1 = &distr_org_struct_1.basisFuncListForBatchs_map[0];
  const int* basisFuncListForBatchs_map_2 = &distr_org_struct_2.basisFuncListForBatchs_map[0];
  const int* basisFuncList_1 = &distr_org_struct_1.basisFuncList[0];
  int basisFuncList_1_count = distr_org_struct_1.basisFuncList.size();
  const int* basisFuncList_2 = &distr_org_struct_2.basisFuncList[0];
  int basisFuncList_2_count = distr_org_struct_2.basisFuncList.size();

  const i_j_val_struct* spMatElementList_1 = &distr_org_struct_1.spMatElementList[0];
  const int* spMatCountList_1 = &distr_org_struct_1.spMatCountList[0];
  const int* spMatIdxList_1 = &distr_org_struct_1.spMatIdxList[0];
  const i_j_val_struct* spMatElementList_2 = &distr_org_struct_2.spMatElementList[0];
  const int* spMatCountList_2 = &distr_org_struct_2.spMatCountList[0];
  const int* spMatIdxList_2 = &distr_org_struct_2.spMatIdxList[0];

  // Set up "partial box-box density matrix"
  int nnn1 = basisFuncList_1_count;
  int nnn2 = basisFuncList_2_count;

  ergo_real* partial_dmat_1 = bufferStructPtr->partial_dmat_1;
  ergo_real* partial_dmat_2 = bufferStructPtr->partial_dmat_2;

  // If dens_CSR is NULL then we assume that partial_dmat_1 and
  // partial_dmat_2 are already prepared, otherwise we prepare them
  // now.
  if(dens_CSR != NULL) {
    for(int i1 = 0; i1 < nnn1; i1++)
      for(int i2 = 0; i2 < nnn2; i2++) {
	int a = basisFuncList_1[i1];
	int b = basisFuncList_2[i2];
	partial_dmat_1[i1*nnn2+i2] = ergo_CSR_get_element_inline(dens_CSR, a, b);
      }
    if(symmetryFlag == 0) {
      for(int i1 = 0; i1 < nnn1; i1++)
	for(int i2 = 0; i2 < nnn2; i2++) {
	  int a = basisFuncList_1[i1];
	  int b = basisFuncList_2[i2];
	  partial_dmat_2[i1*nnn2+i2] = ergo_CSR_get_element_inline(dens_CSR, b, a);
	}
    }
  }

  ergo_real* partial_K_1 = bufferStructPtr->partial_K_1;
  ergo_real* partial_K_2 = bufferStructPtr->partial_K_2;

  for(int i1 = 0; i1 < nnn1; i1++)
    for(int i2 = 0; i2 < nnn2; i2++) {
      partial_K_1[i1*nnn2+i2] = 0;
      if(symmetryFlag == 0)
	partial_K_2[i1*nnn2+i2] = 0;
    }

  for(int batch_i = 0; batch_i < nBatchs_1; batch_i++) {
    int batch_j_start = 0;
    if(interactionWithSelf == 1)
      batch_j_start = batch_i;
    for(int batch_j = batch_j_start; batch_j < nBatchs_2; batch_j++) {
      int noOfBasisFuncPairs_1 = batchList_1[batch_i].noOfBasisFuncPairs;
      int noOfBasisFuncPairs_2 = batchList_2[batch_j].noOfBasisFuncPairs;
      // set integral list to zero
      memset(summedIntegralList, 0, noOfBasisFuncPairs_1*noOfBasisFuncPairs_2*sizeof(ergo_real));

      // Set up "local density matrix" for this pair of batchs.
      int nn1 = batchList_1[batch_i].basisFuncForBatchCount;
      int nn2 = batchList_2[batch_j].basisFuncForBatchCount;
      ergo_real local_dmat_1[nn1][nn2];
      ergo_real local_dmat_2[nn1][nn2];
      ergo_real maxabsdmatelement = 0;
      for(int i1 = 0; i1 < nn1; i1++)
	for(int i2 = 0; i2 < nn2; i2++) {
	  int a2 =  basisFuncListForBatchs_map_1[batchList_1[batch_i].basisFuncForBatchsIndex+i1];
	  int b2 =  basisFuncListForBatchs_map_2[batchList_2[batch_j].basisFuncForBatchsIndex+i2];
	  local_dmat_1[i1][i2] = partial_dmat_1[a2*nnn2+b2];
	  if(symmetryFlag == 0)
	    local_dmat_2[i1][i2] = partial_dmat_2[a2*nnn2+b2];
	  ergo_real absval = template_blas_fabs(local_dmat_1[i1][i2]);
	  if(absval > maxabsdmatelement)
	    maxabsdmatelement = absval;
	  if(symmetryFlag == 0) {
	    ergo_real absval = template_blas_fabs(local_dmat_2[i1][i2]);
	    if(absval > maxabsdmatelement)
	      maxabsdmatelement = absval;
	  }
	}

      int cluster_i_start = batchList_1[batch_i].clusterStartIndex;
      int clusterCount1 = batchList_1[batch_i].noOfClusters;
      for(int cluster_i = cluster_i_start; cluster_i < cluster_i_start + clusterCount1; cluster_i++) {
	int cluster_j_start = batchList_2[batch_j].clusterStartIndex;
	int clusterCount2 = batchList_2[batch_j].noOfClusters;
	int cluterIndexEnd2 = cluster_j_start + clusterCount2;
	if(interactionWithSelf == 1 && batch_i == batch_j)
	  cluster_j_start = cluster_i;
	for(int cluster_j = cluster_j_start; cluster_j < cluterIndexEnd2; cluster_j++) {
	  // check if we can skip this combination of clusters
	  if(clusterList_1[cluster_i].maxLimitingFactorForCluster * clusterList_2[cluster_j].maxLimitingFactorForCluster * maxabsdmatelement < threshold)
	    continue;

	  if(use_multipole_screening_for_clusters == 1) {
	    // Try multipole screening
	    int maxDegree = 2;
	    ergo_real maxAbsContributionFromMultipole = integralInfo.GetMMLimitTable().get_max_abs_mm_contrib(maxDegree,
													      clusterList_1[cluster_i].multipoleEuclNormListForK,
													      maxDegree,
													      clusterList_2[cluster_j].multipoleEuclNormListForK,
													      boxDistance);
	    if(maxAbsContributionFromMultipole * maxabsdmatelement < threshold)
	      continue;
	  } // END IF try multipole screening
		  
	  int group_i_start = clusterList_1[cluster_i].groupStartIndex;
	  int group_i_end = group_i_start + clusterList_1[cluster_i].noOfGroups;
	  int group_j_start = clusterList_2[cluster_j].groupStartIndex;
	  int group_j_end = group_j_start + clusterList_2[cluster_j].noOfGroups;

	  int n1max = clusterList_1[cluster_i].nmax;
	  int n2max = clusterList_2[cluster_j].nmax;
	  int noOfMonomials_1 = integralInfo.monomial_info.no_of_monomials_list[n1max];
	  int noOfMonomials_2 = integralInfo.monomial_info.no_of_monomials_list[n2max];

	  // Now we can precompute things that depend only on exponents
	  ergo_real alpha_1 = groupList_1[group_i_start].exponent;
	  ergo_real alpha_2 = groupList_2[group_j_start].exponent;
	  ergo_real alphasum = alpha_1 + alpha_2;
	  ergo_real alphaproduct = alpha_1 * alpha_2;
	  ergo_real alpha_0 = alphaproduct / alphasum;
	  ergo_real resultPreFactor = twoTimesPiToPow5half / (alphaproduct*template_blas_sqrt(alphasum));

	  for(int group_i = group_i_start; group_i < group_i_end; group_i++) {
	    if(interactionWithSelf == 1 && batch_i == batch_j && cluster_i == cluster_j)
	      group_j_start = group_i;
	    for(int group_j = group_j_start; group_j < group_j_end; group_j++) {
	      // Try Cauchy-Schwartz screening
	      if(groupList_1[group_i].maxLimitingFactorGroup * groupList_2[group_j].maxLimitingFactorGroup * maxabsdmatelement < threshold)
		continue;

	      ergo_real dx = groupList_2[group_j].centerCoords[0] - groupList_1[group_i].centerCoords[0];
	      ergo_real dy = groupList_2[group_j].centerCoords[1] - groupList_1[group_i].centerCoords[1];
	      ergo_real dz = groupList_2[group_j].centerCoords[2] - groupList_1[group_i].centerCoords[2];

	      // Check if multipole screening can be used
	      ergo_real distance = template_blas_sqrt(dx*dx+dy*dy+dz*dz);
	      if(distance > groupList_1[group_i].maxExtentGroup + groupList_2[group_j].maxExtentGroup) {
		// Try multipole screening
		int maxDegree = 2;
		ergo_real maxAbsContributionFromMultipole = integralInfo.GetMMLimitTable().get_max_abs_mm_contrib(maxDegree,
														  groupList_1[group_i].multipoleEuclNormListForK,
														  maxDegree,
														  groupList_2[group_j].multipoleEuclNormListForK,
														  distance);
		if(maxAbsContributionFromMultipole * maxabsdmatelement < threshold)
		  continue;
	      } // END IF try multipole screening

	      // now we can do all integrals needed for this pair of groups
	      // now we have dx dy dz alpha0 alpha1 n1max n2max. Get all integrals for this case.
	      get_related_integrals_hermite(integralInfo,
					    CAM_params,
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
	for(int idx_2 = 0; idx_2 < noOfBasisFuncPairs_2; idx_2++) {

	  int a =     basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_1;
	  int b =     basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_2;
	  int c =     basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_1;
	  int d =     basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_2;
	  ergo_real integralValueCurr = summedIntegralList[idx_1*noOfBasisFuncPairs_2 + idx_2];

	  int a_mod = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_1_mod;
	  int b_mod = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_2_mod;
	  int c_mod = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_1_mod;
	  int d_mod = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_2_mod;

	  int a_mod2 = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_inbox_1;
	  int b_mod2 = basisFuncPairList_1[batchList_1[batch_i].basisFuncPairListIndex+idx_1].index_inbox_2;
	  int c_mod2 = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_inbox_1;
	  int d_mod2 = basisFuncPairList_2[batchList_2[batch_j].basisFuncPairListIndex+idx_2].index_inbox_2;

	  if(a == c && b == d)
	    integralValueCurr *= 2;

	  if(template_blas_fabs(integralValueCurr)*maxabsdmatelement < threshold)
	    continue;

	  ergo_real dens_ac = local_dmat_1[a_mod][c_mod];
	  ergo_real dens_ad = local_dmat_1[a_mod][d_mod];
	  ergo_real dens_bc = local_dmat_1[b_mod][c_mod];
	  ergo_real dens_bd = local_dmat_1[b_mod][d_mod];

	  ergo_real dens_ca = dens_ac;
	  ergo_real dens_da = dens_ad;
	  ergo_real dens_cb = dens_bc;
	  ergo_real dens_db = dens_bd;
		
	  if(symmetryFlag == 0) {
	    dens_ca = local_dmat_2[a_mod][c_mod];
	    dens_da = local_dmat_2[a_mod][d_mod];
	    dens_cb = local_dmat_2[b_mod][c_mod];
	    dens_db = local_dmat_2[b_mod][d_mod];
	  }

	  if(symmetryFlag) {

	    if(a != b && c != d && a != c && a != d && b != c && b != d) {
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
	    }
	    else if(a == b && c != d && a != c && a != d && b != c && b != d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
	    }
	    else if(a != b && c == d && a != c && a != d && b != c && b != d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
	    }
	    else if(a != b && c != d && a == c && a != d && b != c && b != d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr * 2.0;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
	    }
	    else if(a != b && c != d && a != c && a == d && b != c && b != d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr * 2.0;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
	    }
	    else if(a != b && c != d && a != c && a != d && b == c && b != d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr * 2.0;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
	    }
	    else if(a != b && c != d && a != c && a != d && b != c && b == d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr * 2.0;
	    }
	    else if(a != b && c != d && a == c && b == d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
	    }
	    else if(a == b && c == d && a != c && a != d && b != c && b != d) { // OK
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	    }
	    else if(a == b && c == d && a == c && a == d && b == c && b == d) { // OK
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	    }
	    else if(a == b && c != d && a == c && a != d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr*2.0;
	    }
	    else if(a == b && c != d && a != c && a == d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr*2.0;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	    }

	    else if(a != b && c == d && a == c && b != d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr*2.0;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
	    }
	    else if(a != b && c == d && b == c && a != d) { // OK
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr*2.0;
	    }
	    else {
	      return -1;
	    }

	  }
	  else if(a != b && c != d && a != c && a != d && b != c && b != d) {
	    if(symmetryFlag) {
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;
	    }
	    else {
	      partial_K_1[a_mod2*nnn2+d_mod2] += -0.5 * dens_bc * integralValueCurr;
	      partial_K_1[a_mod2*nnn2+c_mod2] += -0.5 * dens_bd * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+c_mod2] += -0.5 * dens_ad * integralValueCurr;
	      partial_K_1[b_mod2*nnn2+d_mod2] += -0.5 * dens_ac * integralValueCurr;

	      partial_K_2[a_mod2*nnn2+d_mod2] += -0.5 * dens_cb * integralValueCurr;
	      partial_K_2[a_mod2*nnn2+c_mod2] += -0.5 * dens_db * integralValueCurr;
	      partial_K_2[b_mod2*nnn2+c_mod2] += -0.5 * dens_da * integralValueCurr;
	      partial_K_2[b_mod2*nnn2+d_mod2] += -0.5 * dens_ca * integralValueCurr;
	    }
	  }
	  else {
	    /* general case, not covered by special cases above (this should be rare, probably not performance-critical). */
	    abcd_struct list[8];
		    
	    /* determine unique configurations */
	    set_abcd_list_item_macro(0, a, b, c, d, dens_bc, a_mod2, d_mod2);
	    set_abcd_list_item_macro(1, a, b, d, c, dens_bd, a_mod2, c_mod2);
	    set_abcd_list_item_macro(2, b, a, c, d, dens_ac, b_mod2, d_mod2);
	    set_abcd_list_item_macro(3, b, a, d, c, dens_ad, b_mod2, c_mod2);

	    set_abcd_list_item_macro(4, c, d, a, b, dens_da, b_mod2, c_mod2);
	    set_abcd_list_item_macro(5, d, c, a, b, dens_ca, b_mod2, d_mod2);
	    set_abcd_list_item_macro(6, c, d, b, a, dens_db, a_mod2, c_mod2);
	    set_abcd_list_item_macro(7, d, c, b, a, dens_cb, a_mod2, d_mod2);

	    int ccc = 0;
	  
	    for(int ii = 0; ii < 8; ii++) {
	      abcd_struct* abcd = &list[ii];
	      int aa, dd;

	      /* check if this is a new unique configuration */
	      int unique = 1;
	      for(int jj = 0; jj < ii; jj++) {
		if(abcd->a == list[jj].a && 
		   abcd->b == list[jj].b && 
		   abcd->c == list[jj].c && 
		   abcd->d == list[jj].d)
		  unique = 0;
	      }
	      if(unique == 0)
		continue;
	      /* now we know that this configuration is unique. */
	      aa = abcd->a;
	      dd = abcd->d;

	      ccc++;

	      if(symmetryFlag) {
		if(dd >= aa)
		  {
		    partial_K_1[abcd->idx1*nnn2+abcd->idx2] += -0.5 * abcd->densValue * integralValueCurr;
		  }
	      }
	      else {
		if(ii <= 3)
		  partial_K_1[abcd->idx1*nnn2+abcd->idx2] += -0.5 * abcd->densValue * integralValueCurr;
		else
		  partial_K_2[abcd->idx1*nnn2+abcd->idx2] += -0.5 * abcd->densValue * integralValueCurr;
	      }

	    } /* END FOR ii go through 8 configurations */
	  }

	} // END FOR idx_1 idx_2
    } // END FOR batch_j
  } // END FOR batch_i

  if(result_K_CSR_shared) {
    if(result_K_CSR_shared->n) {
      pthread_mutex_lock(&K_CSR_shared_access_mutex);
      // Now move results from partial_K to K.
      for(int i1 = 0; i1 < nnn1; i1++)
	for(int i2 = 0; i2 < nnn2; i2++) {
	  int a = basisFuncList_1[i1];
	  int b = basisFuncList_2[i2];
	  ergo_CSR_add_to_element(result_K_CSR_shared,
				  a,
				  b,
				  partial_K_1[i1*nnn2+i2]);
	  if(symmetryFlag == 0) {
	    ergo_CSR_add_to_element(result_K_CSR_shared,
				    b,
				    a,
				    partial_K_2[i1*nnn2+i2]);
	  }
	}
      pthread_mutex_unlock(&K_CSR_shared_access_mutex);
    }
  } // end if use result_K_CSR_shared
  else {
    // Place results in resultMatContrib
    if(resultMatContrib == NULL)
      return -1;
    for(int i1 = 0; i1 < nnn1; i1++)
      for(int i2 = 0; i2 < nnn2; i2++) {
	int a = basisFuncList_1[i1];
	int b = basisFuncList_2[i2];
	if(symmetryFlag == 1) {
	  if(a <= b)
	    resultMatContrib->addContrib(a, b, partial_K_1[i1*nnn2+i2]);
	  else
	    resultMatContrib->addContrib(b, a, partial_K_1[i1*nnn2+i2]);
	}
	else {
	  resultMatContrib->addContrib(a, b, partial_K_1[i1*nnn2+i2]);
	  resultMatContrib->addContrib(b, a, partial_K_2[i1*nnn2+i2]);
	}
      }
  }

  return 0;
}

