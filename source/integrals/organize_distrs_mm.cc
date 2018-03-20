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

/** @file organize_distrs_mm.cc

    @brief Code for organizing a given set of primitive Gaussian
    distributions (typically coming from basis function products)
    regarding information related to multipole methods.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "organize_distrs_mm.h"
#include "serialization_tools.h"
#include <stdexcept>

/* distr_org_mm_struct functions */

distr_org_mm_struct::Data::Data() : chargeSum(0) {
  multipole.degree = -1;
  multipole.noOfMoments = 0;
  memset(multipole.momentList, 0, MAX_NO_OF_MOMENTS_PER_MULTIPOLE*sizeof(ergo_real));
  memset(&multipolePoint, 0x00, 3*sizeof(ergo_real));
  memset(maxMomentVectorNormForDistrsList, 0, (MAX_MULTIPOLE_DEGREE_BASIC+1)*sizeof(ergo_real));
}

void distr_org_mm_struct::writeToBuffer(char* dataBuffer, size_t const bufferSize) const {
  assert(bufferSize >= getSize());
  char* p = dataBuffer;
  memcpy(p, &data, sizeof(data));
  p += sizeof(data);
  std_vector_writeToBuffer_and_move_ptr(multipoleListForGroups, p);
  std_vector_writeToBuffer_and_move_ptr(multipoleListForDistrs, p);
}

size_t distr_org_mm_struct::getSize() const {
  size_t size = sizeof(distr_org_mm_struct::Data);
  size += std_vector_getSize(multipoleListForGroups);
  size += std_vector_getSize(multipoleListForDistrs);
  return size;
}

void distr_org_mm_struct::assignFromBuffer(char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  size_t remainingBytes = bufferSize;
  assert(remainingBytes >= sizeof(data));
  memcpy(&data, p, sizeof(data));
  p += sizeof(data);
  const char* bufEndPtr = &dataBuffer[bufferSize];
  std_vector_assignFromBuffer_and_move_ptr(multipoleListForGroups, p, bufEndPtr);
  std_vector_assignFromBuffer_and_move_ptr(multipoleListForDistrs, p, bufEndPtr);
}

/* distr_list_description_struct functions */

void distr_list_description_struct::writeToBuffer(char* dataBuffer, size_t const bufferSize) const {
  assert(bufferSize >= getSize());
  char* p = dataBuffer;
  org.writeToBuffer(p, org.getSize());
  p += org.getSize();
  org_mm.writeToBuffer(p, org_mm.getSize());
  p += org_mm.getSize();
  assert((size_t)(p - dataBuffer) <= bufferSize);
}

size_t distr_list_description_struct::getSize() const {
  return org.getSize() + org_mm.getSize();
}

void distr_list_description_struct::assignFromBuffer(char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  org.assignFromBuffer(p, bufferSize);
  p += org.getSize();
  org_mm.assignFromBuffer(p, bufferSize - org.getSize());
  p += org_mm.getSize();
  assert((size_t)(p - dataBuffer) == bufferSize);
}

/* ************************************************** */

int
generate_multipoles_for_groups(const IntegralInfo & integralInfo,
			       const distr_org_struct & org,
			       distr_org_mm_struct & result_org_mm,
			       ergo_real* averagePosList,
			       int & avgPosCounter
			       ) {
  ergo_real chargeSum = 0;

  const batch_struct* batchList = &org.batchList[0];
  const cluster_struct* clusterList = &org.clusterList[0];
  const distr_group_struct* groupList = &org.groupList[0];
  const minimal_distr_struct* minimalDistrList = &org.minimalDistrList[0];
  int batchCount = org.batchList.size();
  const basis_func_pair_struct* basisFuncPairList = &org.basisFuncPairList[0];

  ergo_real* maxMomentVectorNormForDistrsList = result_org_mm.data.maxMomentVectorNormForDistrsList;
  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
    maxMomentVectorNormForDistrsList[l] = 0;

  int groupCount = org.groupList.size();
  result_org_mm.multipoleListForGroups.resize(groupCount);

  for(int batchIndex = 0; batchIndex < batchCount; batchIndex++) {
    int clusterCount = batchList[batchIndex].noOfClusters;
    int cluster_start = batchList[batchIndex].clusterStartIndex;
    for(int clusterIndex = cluster_start; clusterIndex < cluster_start + clusterCount; clusterIndex++) {
      int group_start = clusterList[clusterIndex].groupStartIndex;
      int group_end = group_start + clusterList[clusterIndex].noOfGroups;
      for(int groupIndex = group_start; groupIndex < group_end; groupIndex++) {
	const distr_group_struct* currGroup = &groupList[groupIndex];

	// Now create a single multipole description of the density of this group.
	multipole_struct_small* multipoleCurrGroup = &result_org_mm.multipoleListForGroups[groupIndex];
	multipoleCurrGroup->degree = -1;
	multipoleCurrGroup->noOfMoments = 0;
	multipoleCurrGroup->centerCoords[0] = currGroup->centerCoords[0];
	multipoleCurrGroup->centerCoords[1] = currGroup->centerCoords[1];
	multipoleCurrGroup->centerCoords[2] = currGroup->centerCoords[2];
	memset(multipoleCurrGroup->momentList, 0, MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC*sizeof(ergo_real));

	int distr_start = currGroup->startIndex;
	int distr_end = distr_start + currGroup->distrCount;
	for(int distrIndex = distr_start; distrIndex < distr_end; distrIndex++) {
	  int basisFuncPairIndex = minimalDistrList[distrIndex].basisFuncPairIndex;
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

	  // add this multipole to multipole for group.
	  int a = basisFuncPairList[batchList[batchIndex].basisFuncPairListIndex+basisFuncPairIndex].index_1;
	  int b = basisFuncPairList[batchList[batchIndex].basisFuncPairListIndex+basisFuncPairIndex].index_2;
	  ergo_real factor = basisFuncPairList[batchList[batchIndex].basisFuncPairListIndex+basisFuncPairIndex].dmatElement;
	  if(a != b)
	    factor *= 2;

	  for(int l = 0; l <= multipole.degree; l++) {
	    int startIndex = l*l;
	    int endIndex = (l+1)*(l+1);
	    ergo_real sum = 0;
	    for(int A = startIndex; A < endIndex; A++)
	      sum += multipole.momentList[A]*multipole.momentList[A];
	    ergo_real subNorm = template_blas_sqrt(sum);
	    if(subNorm > maxMomentVectorNormForDistrsList[l])
	      maxMomentVectorNormForDistrsList[l] = subNorm;
	  }
			  
	  for(int kk = 0; kk < multipole.noOfMoments; kk++)
	    multipoleCurrGroup->momentList[kk] += factor * multipole.momentList[kk];
	  if(multipole.degree > multipoleCurrGroup->degree)
	    multipoleCurrGroup->degree = multipole.degree;
	  if(multipole.noOfMoments > multipoleCurrGroup->noOfMoments)
	    multipoleCurrGroup->noOfMoments = multipole.noOfMoments;
	} // END FOR distrIndex

	// OK, multipoleCurrGroup is complete.
	chargeSum += multipoleCurrGroup->momentList[0];
	for(int kk = 0; kk < 3; kk++)
	  averagePosList[kk] += multipoleCurrGroup->centerCoords[kk];
	avgPosCounter++;

      } // END FOR groupIndex
    } // END FOR clusterIndex
  } // END FOR batchIndex

  return 0;
}


int
get_multipole_pt_for_box(const ergo_real* boxCenterCoords,
			 ergo_real boxWidth,
			 const ergo_real* averagePosList,
			 int avgPosCounter,
			 ergo_real* resultMultipolePoint) {
  // use average position instead of center-of-charge, 
  // because center-of-charge is ill-defined when some charges are negative.
  if(avgPosCounter == 0) {
    for(int kk = 0; kk < 3; kk++)
      resultMultipolePoint[kk] = boxCenterCoords[kk];
  }
  else {
    for(int kk = 0; kk < 3; kk++)
      resultMultipolePoint[kk] = averagePosList[kk] / avgPosCounter;
  }
  // check that "resultMultipolePoint" is not too far from box center.
  ergo_real sumofsquares = 0;
  for(int kk = 0; kk < 3; kk++) {
    ergo_real dx = resultMultipolePoint[kk] - boxCenterCoords[kk];
    sumofsquares += dx*dx;
  }
  ergo_real distFromCenter = template_blas_sqrt(sumofsquares);
  if(distFromCenter > boxWidth) {
    do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_multipole_pt_for_box: (distFromCenter > boxWidth).");
    return -1;
  }
  return 0;
}


int
translate_multipoles_for_box(distr_org_mm_struct & result_org_mm,
			     const distr_org_struct & org,
			     const MMTranslator & translator
			     ) {

  ergo_real* multipolePointCoords = result_org_mm.data.multipolePoint;
  multipole_struct_large branchMultipole;
  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
    branchMultipole.momentList[A] = 0;
  for(int kk = 0; kk < 3; kk++)
    branchMultipole.centerCoords[kk] = multipolePointCoords[kk];
  branchMultipole.degree = MAX_MULTIPOLE_DEGREE;
  branchMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

  const batch_struct* batchList = &org.batchList[0];
  const cluster_struct* clusterList = &org.clusterList[0];
  const multipole_struct_small* multipoleListForGroups = &result_org_mm.multipoleListForGroups[0];

  int batchCount = org.batchList.size();
  for(int batchIndex = 0; batchIndex < batchCount; batchIndex++) {
    int clusterCount = batchList[batchIndex].noOfClusters;
    int cluster_start = batchList[batchIndex].clusterStartIndex;
    for(int clusterIndex = cluster_start; clusterIndex < cluster_start + clusterCount; clusterIndex++) {
      int group_start = clusterList[clusterIndex].groupStartIndex;
      int group_end = group_start + clusterList[clusterIndex].noOfGroups;
      for(int groupIndex = group_start; groupIndex < group_end; groupIndex++) {
	const multipole_struct_small & multipoleCurrGroup = multipoleListForGroups[groupIndex];

	// take multipole for this group, and translate it to center-of-charge point
	ergo_real dx = multipoleCurrGroup.centerCoords[0] - multipolePointCoords[0];
	ergo_real dy = multipoleCurrGroup.centerCoords[1] - multipolePointCoords[1];
	ergo_real dz = multipoleCurrGroup.centerCoords[2] - multipolePointCoords[2];

	ergo_real W[MAX_NO_OF_MOMENTS_PER_MULTIPOLE*MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
	translator.getTranslationMatrix
	  (dx, dy, dz, MAX_MULTIPOLE_DEGREE,
	   multipoleCurrGroup.degree, W);

	multipole_struct_large translatedMultipole;
	for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++) {
	  ergo_real sum = 0;
	  for(int B = 0; B < multipoleCurrGroup.noOfMoments; B++)
	    sum += W[A*multipoleCurrGroup.noOfMoments+B] * multipoleCurrGroup.momentList[B];
	  translatedMultipole.momentList[A] = sum;
	} // END FOR A
	for(int kk = 0; kk < 3; kk++)
	  translatedMultipole.centerCoords[kk] = multipolePointCoords[kk];
	translatedMultipole.degree = MAX_MULTIPOLE_DEGREE;
	translatedMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

	// add translated multipole to branch multipole
	for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
	  branchMultipole.momentList[A] += translatedMultipole.momentList[A];
      } // END FOR groupIndex
    } // END FOR clusterIndex
  } // END FOR batchIndex
  setup_multipole_maxAbsMomentList(&branchMultipole);
  result_org_mm.data.multipole = branchMultipole;

  return 0;
}



int
combine_mm_info_for_child_boxes(distr_list_description_struct & result_box_branch,
				const distr_list_description_struct** child_box_branches,
				int noOfChildren,
				const MMTranslator & translator) {
  multipole_struct_large & newMultipole = result_box_branch.org_mm.data.multipole;
  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
    newMultipole.momentList[A] = 0;

  // get average position of child multipoles
  ergo_real avgPosList[3];
  for(int kk = 0; kk < 3; kk++)
    avgPosList[kk] = 0;

  ergo_real* maxMomentVectorNormForDistrsList = result_box_branch.org_mm.data.maxMomentVectorNormForDistrsList;
  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
    maxMomentVectorNormForDistrsList[l] = 0;

  for(int childIndex = 0; childIndex < noOfChildren; childIndex++) {
    for(int kk = 0; kk < 3; kk++)
      avgPosList[kk] += child_box_branches[childIndex]->org_mm.data.multipole.centerCoords[kk];
  } // END FOR childIndex

  for(int kk = 0; kk < 3; kk++)
    newMultipole.centerCoords[kk] = avgPosList[kk] / noOfChildren;
  newMultipole.degree = MAX_MULTIPOLE_DEGREE;
  newMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

  // We also want to get maxExtent and maxDistanceOutsideBox for parent box (use largest values found among the children).
  ergo_real maxExtent = 0;
  ergo_real maxDistanceOutsideBox = 0;

  // Now translate child multipoles and add to parent multipole
  for(int childIndex = 0; childIndex < noOfChildren; childIndex++) {
    const multipole_struct_large* childMultipole = &child_box_branches[childIndex]->org_mm.data.multipole;

    if(child_box_branches[childIndex]->org.data.maxExtent > maxExtent)
      maxExtent = child_box_branches[childIndex]->org.data.maxExtent;

    if(child_box_branches[childIndex]->org.data.maxDistanceOutsideBox > maxDistanceOutsideBox)
      maxDistanceOutsideBox = child_box_branches[childIndex]->org.data.maxDistanceOutsideBox;

    ergo_real dx = childMultipole->centerCoords[0] - newMultipole.centerCoords[0];
    ergo_real dy = childMultipole->centerCoords[1] - newMultipole.centerCoords[1];
    ergo_real dz = childMultipole->centerCoords[2] - newMultipole.centerCoords[2];

    ergo_real W[MAX_NO_OF_MOMENTS_PER_MULTIPOLE*MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
    translator.getTranslationMatrix(dx, dy, dz,
				    MAX_MULTIPOLE_DEGREE,
				    MAX_MULTIPOLE_DEGREE, W);

    multipole_struct_large translatedMultipole;
    for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++) {
      ergo_real sum = 0;
      for(int B = 0; B < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; B++)
	sum += W[A*MAX_NO_OF_MOMENTS_PER_MULTIPOLE+B] * childMultipole->momentList[B];
      translatedMultipole.momentList[A] = sum;
    } // END FOR A
    for(int kk = 0; kk < 3; kk++)
      translatedMultipole.centerCoords[kk] = newMultipole.centerCoords[kk];
    translatedMultipole.degree = MAX_MULTIPOLE_DEGREE;
    translatedMultipole.noOfMoments = MAX_NO_OF_MOMENTS_PER_MULTIPOLE;

    // add translated multipole to parent multipole
    for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
      newMultipole.momentList[A] += translatedMultipole.momentList[A];

    for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++) {
      ergo_real childValue = child_box_branches[childIndex]->org_mm.data.maxMomentVectorNormForDistrsList[l];
      if(childValue > maxMomentVectorNormForDistrsList[l])
	maxMomentVectorNormForDistrsList[l] = childValue;
    }

  } // END FOR childIndex

  setup_multipole_maxAbsMomentList(&newMultipole);

  result_box_branch.org.data.maxExtent = maxExtent;
  result_box_branch.org.data.maxDistanceOutsideBox = maxDistanceOutsideBox;

  return 0;
}
