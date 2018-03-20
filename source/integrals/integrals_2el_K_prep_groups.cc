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

/** @file integrals_2el_K_prep_groups.cc

    \brief Code for preparing basis function group information to be
    used for computing the Hartree-Fock exchange matrix K.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include "integrals_2el_K_prep_groups.h"

int prep_info_for_K(int maxCount,
		    distr_org_struct & org,
		    int distrCountCurrBox,
		    const ergo_real* multipoleNormVectorList,
		    const int* multipoleDegreeList,
		    const ergo_real* limitingFactorList,
		    const int* basisFuncGroupList1,
		    const int* basisFuncGroupList2) {
  // go through all distrs of this box, and update basisFuncGroupInfoList accordingly.
  std::vector<basis_func_group_info_for_box> basisFuncGroupInfoListForK_tmp(maxCount);
  int count = 0;
  for(int jjj = 0; jjj < distrCountCurrBox; jjj++) {
    const ergo_real* multipoleNormVectorList_curr = NULL;
    if(multipoleNormVectorList)
      multipoleNormVectorList_curr = &multipoleNormVectorList[jjj*(MAX_MULTIPOLE_DEGREE_BASIC+1)];
    int multipoleDegree_curr = 0;
    if(multipoleDegreeList)
      multipoleDegree_curr = multipoleDegreeList[jjj];
    int basisFuncGroup_1 = basisFuncGroupList1[jjj];
    int basisFuncGroup_2 = basisFuncGroupList2[jjj];
    ergo_real CS_factor = limitingFactorList[jjj];
    // check if basisFuncGroup_1 and/or basisFuncGroup_2 is already present
    int foundIndex_1 = -1;
    int foundIndex_2 = -1;
    if(count > maxCount) {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (count > maxCount)");
      return -1;
    }
    for(int k = 0; k < count; k++) {
      if(basisFuncGroupInfoListForK_tmp[k].basisFuncGroupIndex == basisFuncGroup_1)
	foundIndex_1 = k;
      if(basisFuncGroupInfoListForK_tmp[k].basisFuncGroupIndex == basisFuncGroup_2)
	foundIndex_2 = k;
    }
    if(foundIndex_1 >= 0) {
      // check if max_CS_factor needs updating
      if(CS_factor > basisFuncGroupInfoListForK_tmp[foundIndex_1].max_CS_factor)
	basisFuncGroupInfoListForK_tmp[foundIndex_1].max_CS_factor = CS_factor;
      if(multipoleDegreeList) {
	// modfy maxMomentVectorNormList if needed.
	if(multipoleDegree_curr > basisFuncGroupInfoListForK_tmp[foundIndex_1].maxMultipoleDegree)
	  basisFuncGroupInfoListForK_tmp[foundIndex_1].maxMultipoleDegree = multipoleDegree_curr;
	if(multipoleNormVectorList_curr) {
	  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++) {
	    if(multipoleNormVectorList_curr[l] > basisFuncGroupInfoListForK_tmp[foundIndex_1].maxMomentVectorNormList[l])
	      basisFuncGroupInfoListForK_tmp[foundIndex_1].maxMomentVectorNormList[l] = multipoleNormVectorList_curr[l];
	  }
	}
      }
    }
    else {
      // add new entry for basisFuncGroup_1
      if(count >= maxCount) {
	do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (count >= maxCount)");
	return -1;
      }
      basisFuncGroupInfoListForK_tmp[count].basisFuncGroupIndex = basisFuncGroup_1;
      basisFuncGroupInfoListForK_tmp[count].max_CS_factor = CS_factor;
      basisFuncGroupInfoListForK_tmp[count].maxMultipoleDegree = multipoleDegree_curr;
      if(multipoleNormVectorList_curr) {
	for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
	  basisFuncGroupInfoListForK_tmp[count].maxMomentVectorNormList[l] = multipoleNormVectorList_curr[l];
      }
      count++;
    }
    if(basisFuncGroup_2 != basisFuncGroup_1) {
      if(foundIndex_2 >= 0) {
	// check if maxSize needs updating
	if(CS_factor > basisFuncGroupInfoListForK_tmp[foundIndex_2].max_CS_factor)
	  basisFuncGroupInfoListForK_tmp[foundIndex_2].max_CS_factor = CS_factor;
	if(multipoleDegreeList) {
	  // modfy maxMomentVectorNormList if needed.
	  if(multipoleDegree_curr > basisFuncGroupInfoListForK_tmp[foundIndex_2].maxMultipoleDegree)
	    basisFuncGroupInfoListForK_tmp[foundIndex_2].maxMultipoleDegree = multipoleDegree_curr;
	  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++) {
	    if(multipoleNormVectorList_curr[l] > basisFuncGroupInfoListForK_tmp[foundIndex_2].maxMomentVectorNormList[l])
	      basisFuncGroupInfoListForK_tmp[foundIndex_2].maxMomentVectorNormList[l] = multipoleNormVectorList_curr[l];
	  }
	}
      }
      else {
	// add new entry for basisFuncGroup_2
	if(count >= maxCount) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error doing basisFuncGroupInfoList: (count >= maxCount)");
	  return -1;
	}
	basisFuncGroupInfoListForK_tmp[count].basisFuncGroupIndex = basisFuncGroup_2;
	basisFuncGroupInfoListForK_tmp[count].max_CS_factor = CS_factor;
	basisFuncGroupInfoListForK_tmp[count].maxMultipoleDegree = multipoleDegree_curr;
	if(multipoleNormVectorList_curr) {
	  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE_BASIC; l++)
	    basisFuncGroupInfoListForK_tmp[count].maxMomentVectorNormList[l] = multipoleNormVectorList_curr[l];
	}
	count++;
      }
    }
    // OK, distr j done
  } // END FOR j
  org.basisFuncGroupInfoListForK.resize(count);
  for(int i = 0; i < count; i++)
    org.basisFuncGroupInfoListForK[i] = basisFuncGroupInfoListForK_tmp[i];
  return 0;
}
