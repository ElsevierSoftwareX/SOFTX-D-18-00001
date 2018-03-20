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

/** @file integrals_2el_J_mm_utils.cc

    \brief Utility functions related to multipole method, used in
    construction of the Coulomb matrix J.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include "integrals_2el_J_mm_utils.h"
#include "box_system.h"


/* Returns 0 or 1 (or -1 on error) */
int
check_if_multipoles_can_be_used(const IntegralInfo & integralInfo,
				ergo_real threshold,
				const ergo_real* boxCenterCoords_1,
				const ergo_real* boxCenterCoords_2,
				ergo_real boxWidth,
				const distr_org_struct & org_1,
				const distr_org_mm_struct & org_mm_1,
				const distr_org_struct & org_2,
				const distr_org_mm_struct & org_mm_2) {
  // check if multipoles can be used.
  // start by computing the minimum distance between the boxes.
  // We assume that both boxes have the same width.
  bool sameBox = true;
  ergo_real dxList[3];
  for(int coordIndex = 0; coordIndex< 3; coordIndex++) {
    ergo_real x1 = boxCenterCoords_1[coordIndex];
    ergo_real x2 = boxCenterCoords_2[coordIndex];
    ergo_real dx = template_blas_fabs(x1 - x2);
    ergo_real width = boxWidth;
    if(dx > width)
      dxList[coordIndex] = dx - width;
    else
      dxList[coordIndex] = 0;
    // If boxes 1 and 2 are the same we get dxList[coordIndex]=0, so to detect that case we also check if dx is significant, then we know it is not the same box.
    if(dx > width/2)
      sameBox = false;
  }
  ergo_real sumOfSquares = 0;
  for(int coordIndex = 0; coordIndex< 3; coordIndex++)
    sumOfSquares += dxList[coordIndex] * dxList[coordIndex];
  ergo_real distance = template_blas_sqrt(sumOfSquares);
  ergo_real maxDistanceOutsideBox_1 = org_1.data.maxDistanceOutsideBox;
  ergo_real maxDistanceOutsideBox_2 = org_2.data.maxDistanceOutsideBox;
  int useMultipoleDescription = 0;

  if(sameBox == false && distance >= maxDistanceOutsideBox_1 + maxDistanceOutsideBox_2) {
    // The distance is OK.
    // We also want to check that the multipole degree needed is not too high.
    // For that we need max norms of subvectors for distrs of both branches.
    // First the case with distrs of 1 interacting with multipole of 2
    ergo_real r_1 = get_min_distance_from_point_to_box(boxCenterCoords_1, boxWidth / 2,
						       org_mm_2.data.multipole.centerCoords);
    int degreeNeeded_1 =
      integralInfo.GetMMLimitTable().get_minimum_multipole_degree_needed(r_1,
									 &org_mm_2.data.multipole,
									 MAX_MULTIPOLE_DEGREE_BASIC,
									 org_mm_1.data.maxMomentVectorNormForDistrsList,
									 threshold);
    if(degreeNeeded_1 < 0)
      return -1;
    // Now the case with distrs of 2 interacting with multipole of 1
    ergo_real r_2 = get_min_distance_from_point_to_box(boxCenterCoords_2, boxWidth / 2,
						       org_mm_1.data.multipole.centerCoords);
    int degreeNeeded_2 =
      integralInfo.GetMMLimitTable().get_minimum_multipole_degree_needed(r_2,
									 &org_mm_1.data.multipole,
									 MAX_MULTIPOLE_DEGREE_BASIC,
									 org_mm_2.data.maxMomentVectorNormForDistrsList,
									 threshold);
    if(degreeNeeded_2 < 0)
      return -1;
    // We need some margin compared to MAX_MULTIPOLE_DEGREE, because
    // in some cases the box multipole is alternating between
    // odd/even large/small elements. Therefore we increase
    // degreeNeeded_1 and degreeNeeded_2 by 1 here.
    degreeNeeded_1++;
    degreeNeeded_2++;
    if(degreeNeeded_1 < MAX_MULTIPOLE_DEGREE && degreeNeeded_2 < MAX_MULTIPOLE_DEGREE)
      useMultipoleDescription = 1;
  }

  return useMultipoleDescription;
}


int
create_list_of_multipoles_for_box(const IntegralInfo& integralInfo,
				  const distr_org_struct & org,
				  multipole_struct_small* result_multipoleList) {
  // create list of multipoles
  const batch_struct* batchList = &org.batchList[0];
  const cluster_struct* clusterList = &org.clusterList[0];
  const distr_group_struct* groupList = &org.groupList[0];
  const minimal_distr_struct* minimalDistrList = &org.minimalDistrList[0];
  int batchCount = org.batchList.size();
  int count_temp = 0;
  for(int batchIndex = 0; batchIndex < batchCount; batchIndex++) {
    int clusterCount = batchList[batchIndex].noOfClusters;
    int cluster_start = batchList[batchIndex].clusterStartIndex;
    for(int clusterIndex = cluster_start; clusterIndex < cluster_start + clusterCount; clusterIndex++) {
      int group_start = clusterList[clusterIndex].groupStartIndex;
      int group_end = group_start + clusterList[clusterIndex].noOfGroups;
      for(int groupIndex = group_start; groupIndex < group_end; groupIndex++) {
	const distr_group_struct* currGroup = &groupList[groupIndex];
	int distr_start = currGroup->startIndex;
	int distr_end = distr_start + currGroup->distrCount;
	for(int distrIndex = distr_start; distrIndex < distr_end; distrIndex++) {
	  int monomialIndex = minimalDistrList[distrIndex].monomialIndex;
	  ergo_real coeff = minimalDistrList[distrIndex].coeff;
	  // get monomialInts from monomialIndex
	  DistributionSpecStruct distr;
	  distr.monomialInts[0] = integralInfo.monomial_info.monomial_list[monomialIndex].ix;
	  distr.monomialInts[1] = integralInfo.monomial_info.monomial_list[monomialIndex].iy;
	  distr.monomialInts[2] = integralInfo.monomial_info.monomial_list[monomialIndex].iz;
	  distr.coeff = coeff;
	  distr.exponent = currGroup->exponent;
	  distr.centerCoords[0] = currGroup->centerCoords[0];
	  distr.centerCoords[1] = currGroup->centerCoords[1];
	  distr.centerCoords[2] = currGroup->centerCoords[2];
	  multipole_struct_small multipole;
	  if(compute_multipole_moments(integralInfo, &distr, &multipole) != 0) {
	    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_multipole_moments");
	    return -1;
	  }
	  result_multipoleList[count_temp] = multipole;
	  count_temp++;
	} // END FOR distrIndex
      } // END FOR groupIndex
    } // END FOR clusterIndex
  } // END FOR batchIndex
  return 0;
}
