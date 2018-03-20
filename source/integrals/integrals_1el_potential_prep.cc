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

/** @file integrals_1el_potential_prep.cc

    @brief Code for 1-electron integrals, preparatory work for
    computation of electron-nuclear potential energy matrix V.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "integrals_1el_potential_prep.h"
#include "multipole.h"
#include <stdexcept>

SetOfDistrsForV::SetOfDistrsForV() { }

SetOfDistrsForV::SetOfDistrsForV(const SetOfDistrsForV & other) {
  distrList = other.distrList;
  multipoleList = other.multipoleList;
  groupList = other.groupList;
  maxMomentVectorNormList = other.maxMomentVectorNormList;
  info = other.info;
}

static void copy_data_and_advance_dest_ptr(char** destPtr, const char* srcPtr, size_t nBytes) {
  char* p = *destPtr;
  memcpy(p, srcPtr, nBytes);
  p += nBytes;
  *destPtr = p;
}

/** Function needed for Chunks and Tasks usage. */
void SetOfDistrsForV::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error: bufferSize too small in SetOfDistrsForV::write_to_buffer.");
  int nDistrs = distrList.size();
  int nGroups = groupList.size();
  // nDistrs
  copy_data_and_advance_dest_ptr(&p, (const char*)&nDistrs, sizeof(int));
  // nGroups
  copy_data_and_advance_dest_ptr(&p, (const char*)&nGroups, sizeof(int));
  // info
  copy_data_and_advance_dest_ptr(&p, (const char*)&info, sizeof(SetOfDistrsForVInfo));
  // distrList
  copy_data_and_advance_dest_ptr(&p, (const char*)&distrList[0], nDistrs*sizeof(DistributionSpecStructWithIndexes2));
  // multipoleList
  copy_data_and_advance_dest_ptr(&p, (const char*)&multipoleList[0], nDistrs*sizeof(multipole_struct_small));
  // groupList
  copy_data_and_advance_dest_ptr(&p, (const char*)&groupList[0], nGroups*sizeof(group_struct));
  // maxMomentVectorNormList
  copy_data_and_advance_dest_ptr(&p, (const char*) &maxMomentVectorNormList[0], nGroups*sizeof(maxMomentVectorNormStruct));
  // DONE!
}

/*
  std::vector<DistributionSpecStructWithIndexes2> distrList;
  std::vector<multipole_struct_small> multipoleList; // same size as distrList
  std::vector<group_struct> groupList;
  std::vector<maxMomentVectorNormStruct> maxMomentVectorNormList; // size same as groupList
  SetOfDistrsForVInfo info;
*/

/** Function needed for Chunks and Tasks usage. */
size_t SetOfDistrsForV::get_size() const {
  int nDistrs = distrList.size();
  int nGroups = groupList.size();
  return
    2 * sizeof(int) +
    sizeof(SetOfDistrsForVInfo) +
    nDistrs*sizeof(DistributionSpecStructWithIndexes2) +
    nDistrs*sizeof(multipole_struct_small) +
    nGroups*sizeof(group_struct) +
    nGroups*sizeof(maxMomentVectorNormStruct);
}

static void copy_data_and_advance_src_ptr(char* destPtr, const char** srcPtr, size_t nBytes) {
  const char* p = *srcPtr;
  memcpy(destPtr, p, nBytes);
  p += nBytes;
  *srcPtr = p;
}

/** Function needed for Chunks and Tasks usage. */
void SetOfDistrsForV::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  if(bufferSize < 2*sizeof(int))
    throw std::runtime_error("Error: bufferSize too small in SetOfDistrsForV::assign_from_buffer.");
  const char* p = dataBuffer;
  // nDistrs
  int nDistrs;
  copy_data_and_advance_src_ptr((char*)&nDistrs, &p, sizeof(int));
  // nGroups
  int nGroups;
  copy_data_and_advance_src_ptr((char*)&nGroups, &p, sizeof(int));
  distrList.resize(nDistrs);
  multipoleList.resize(nDistrs);
  groupList.resize(nGroups);
  maxMomentVectorNormList.resize(nGroups);
  assert(bufferSize >= get_size());
  // info
  copy_data_and_advance_src_ptr((char*)&info, &p, sizeof(SetOfDistrsForVInfo));
  // distrList
  copy_data_and_advance_src_ptr((char*)&distrList[0], &p, nDistrs*sizeof(DistributionSpecStructWithIndexes2));
  // multipoleList
  copy_data_and_advance_src_ptr((char*)&multipoleList[0], &p, nDistrs*sizeof(multipole_struct_small));
  // groupList
  copy_data_and_advance_src_ptr((char*)&groupList[0], &p, nGroups*sizeof(group_struct));
  // maxMomentVectorNormList
  copy_data_and_advance_src_ptr((char*)&maxMomentVectorNormList[0], &p, nGroups*sizeof(maxMomentVectorNormStruct));
  // DONE!
}

void
organize_distrs_for_V(const IntegralInfo & integralInfo,
		      SetOfDistrsForV & setOfDistrsForV,
		      const std::vector<DistributionSpecStructWithIndexes2> & inputList,
		      ergo_real threshold,
		      ergo_real maxCharge) {
  int nDistrs = inputList.size();
  setOfDistrsForV.distrList.resize(nDistrs);
  for(int i = 0; i < nDistrs; i++)
    setOfDistrsForV.distrList[i] = inputList[i];

  // Sort list of distrs by x, y, z, exponent.
  // The point of this is to group together distrs that have same center and same exponent.
  sort_distr_list(&setOfDistrsForV.distrList[0], nDistrs);

  // identify groups of distrs that have same center and same exponent.
  // Allocate according to worst case, each distr being a separate group.
  setOfDistrsForV.groupList.resize(nDistrs);
  int ind = 0;
  int currGroupInd = 0;
  int groupCount = 0;
  int maxNDistrsPerGroup = 0;
  while(ind < nDistrs) {
    ind++;
    if(ind < nDistrs) {
      if(compare_distrs<DistributionSpecStructWithIndexes2>(&setOfDistrsForV.distrList[ind], &setOfDistrsForV.distrList[currGroupInd]) == 0)
	continue;
    }
    // define new group
    setOfDistrsForV.groupList[groupCount].startIndex = currGroupInd;
    setOfDistrsForV.groupList[groupCount].count = ind - currGroupInd;
    if (setOfDistrsForV.groupList[groupCount].count > maxNDistrsPerGroup)
      maxNDistrsPerGroup = setOfDistrsForV.groupList[groupCount].count;
    groupCount++;
    // start next group
    currGroupInd = ind;
  }
  setOfDistrsForV.groupList.resize(groupCount);
  setOfDistrsForV.maxMomentVectorNormList.resize(groupCount);

  // Create multipoles for all distrs.
  setOfDistrsForV.multipoleList.resize(nDistrs);
  memset(&setOfDistrsForV.multipoleList[0], 0, nDistrs*sizeof(multipole_struct_small));
  for(int i = 0; i < nDistrs; i++)
    compute_multipole_moments(integralInfo, &setOfDistrsForV.distrList[i].distr, &setOfDistrsForV.multipoleList[i]);

  // Determine min and max coords for all distrs, to set boundingCubeCenterCoords and boundingCubeWidth.
  ergo_real minCoords[3];
  ergo_real maxCoords[3];
  for(int i = 0; i < nDistrs; i++) {
    for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
      ergo_real coord = setOfDistrsForV.distrList[i].distr.centerCoords[coordIdx];
      if(i == 0 || coord < minCoords[coordIdx])
	minCoords[coordIdx] = coord;
      if(i == 0 || coord > maxCoords[coordIdx])
	maxCoords[coordIdx] = coord;
    }
  }
  ergo_real maxWidth = 0;
  for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
    ergo_real centerCoord = (minCoords[coordIdx] + maxCoords[coordIdx]) / 2;
    setOfDistrsForV.info.boundingCubeCenterCoords[coordIdx] = centerCoord;
    ergo_real width = maxCoords[coordIdx] - minCoords[coordIdx];
    if(width > maxWidth)
      maxWidth = width;
  }
  setOfDistrsForV.info.boundingCubeWidth = maxWidth;

  setOfDistrsForV.info.maxExtentForAll = 0;
  // Get maxMomentVectorNorm info for each group
  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
    setOfDistrsForV.info.maxMomentVectorNormForAll.maxMomentVectorNormList[l] = 0;
  for(int groupIndex = 0; groupIndex < groupCount; groupIndex++) {
    int groupStartIdx = setOfDistrsForV.groupList[groupIndex].startIndex;
    multipole_struct_small* currMultipoleList = &setOfDistrsForV.multipoleList[groupStartIdx];
    int nDistrsCurrGroup = setOfDistrsForV.groupList[groupIndex].count;
    int maxNoOfMoments = 0;
    int maxDegree = 0;
    ergo_real maxExtentForGroup = 0;
    for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
      setOfDistrsForV.maxMomentVectorNormList[groupIndex].maxMomentVectorNormList[l] = 0;
    for(int i = 0; i < nDistrsCurrGroup; i++) {
      if(currMultipoleList[i].noOfMoments > maxNoOfMoments)
	maxNoOfMoments = currMultipoleList[i].noOfMoments;
      if(currMultipoleList[i].degree > maxDegree)
	maxDegree = currMultipoleList[i].degree;
      const multipole_struct_small* distrMultipole = &currMultipoleList[i];
      for(int l = 0; l <= distrMultipole->degree; l++) {
	int startIndex = l*l;
	int endIndex = (l+1)*(l+1);
	ergo_real sum = 0;
	for(int A = startIndex; A < endIndex; A++)
	  sum += distrMultipole->momentList[A]*distrMultipole->momentList[A];
	ergo_real subNorm = template_blas_sqrt(sum);
	if(subNorm > setOfDistrsForV.maxMomentVectorNormList[groupIndex].maxMomentVectorNormList[l])
	  setOfDistrsForV.maxMomentVectorNormList[groupIndex].maxMomentVectorNormList[l] = subNorm;
      }
      // Get extent
      // Here we use an extent such that beyond the extent the abs
      // value of any distr is smaller than threshold/maxCharge.
      ergo_real abscoeff = template_blas_fabs(setOfDistrsForV.distrList[groupStartIdx+i].distr.coeff);
      ergo_real exponent = setOfDistrsForV.distrList[groupStartIdx+i].distr.exponent;
      ergo_real R2 = -1 * (1/exponent) * template_blas_log(threshold/(abscoeff*maxCharge));
      ergo_real extent = 0;
      if(R2 > 0) // R2 can become negative, e.g. if abscoeff is very small, in such cases we let extent be zero.
	extent = template_blas_sqrt(R2);
      if(extent > maxExtentForGroup)
	maxExtentForGroup = extent;
    } // end for i
    setOfDistrsForV.groupList[groupIndex].maxExtent = maxExtentForGroup;
    setOfDistrsForV.groupList[groupIndex].maxNoOfMoments = maxNoOfMoments;
    setOfDistrsForV.groupList[groupIndex].maxDegree = maxDegree;
    if(maxExtentForGroup > setOfDistrsForV.info.maxExtentForAll)
      setOfDistrsForV.info.maxExtentForAll = maxExtentForGroup;
    // Update maxMomentVectorNormForAll
    for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++) {
      ergo_real currValue = setOfDistrsForV.maxMomentVectorNormList[groupIndex].maxMomentVectorNormList[l];
      if(currValue > setOfDistrsForV.info.maxMomentVectorNormForAll.maxMomentVectorNormList[l])
	setOfDistrsForV.info.maxMomentVectorNormForAll.maxMomentVectorNormList[l] = currValue;
    }
  } // end for groupIndex
}

