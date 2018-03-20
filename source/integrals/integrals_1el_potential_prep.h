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

/** @file integrals_1el_potential_prep.h

    @brief Code for 1-electron integrals, preparatory work for
    computation of electron-nuclear potential energy matrix V.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef INTEGRALS_1EL_POTENTIAL_PREP_HEADER
#define INTEGRALS_1EL_POTENTIAL_PREP_HEADER

#include "basisinfo.h"
#include <algorithm>    // std::sort

struct DistributionSpecStructWithIndexes2 {
  DistributionSpecStruct distr;
  int basisFuncIdx1;
  int basisFuncIdx2;
};

struct group_struct {
  int startIndex;
  int count;
  int maxNoOfMoments;
  int maxDegree;
  ergo_real maxExtent;
};

struct maxMomentVectorNormStruct {
  ergo_real maxMomentVectorNormList[MAX_MULTIPOLE_DEGREE_BASIC+1];
};

struct SetOfDistrsForVInfo {
  ergo_real maxExtentForAll;
  maxMomentVectorNormStruct maxMomentVectorNormForAll;
  ergo_real boundingCubeCenterCoords[3];
  ergo_real boundingCubeWidth;
};

struct SetOfDistrsForV {
  std::vector<DistributionSpecStructWithIndexes2> distrList;
  std::vector<multipole_struct_small> multipoleList; // same size as distrList
  std::vector<group_struct> groupList;
  std::vector<maxMomentVectorNormStruct> maxMomentVectorNormList; // size same as groupList
  SetOfDistrsForVInfo info;
  // Stuff needed for Chunks&Tasks usage
  SetOfDistrsForV();
  SetOfDistrsForV(const SetOfDistrsForV & other);
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};

void
organize_distrs_for_V(const IntegralInfo & integralInfo,
		      SetOfDistrsForV & setOfDistrsForV,
		      const std::vector<DistributionSpecStructWithIndexes2> & inputList,
		      ergo_real threshold,
		      ergo_real maxCharge);

template <typename DistributionSpecStructType>
int
compare_distrs(const void* p1, const void* p2) {
  DistributionSpecStructType* d1 = (DistributionSpecStructType*)p1;
  DistributionSpecStructType* d2 = (DistributionSpecStructType*)p2;
  /* FIXME: Not nice to have these two hardcoded values here.  */
  const ergo_real tolernance_dist = 1e-10;
  const ergo_real tolernance_exponent = 1e-11;
  ergo_real dx = d1->distr.centerCoords[0] - d2->distr.centerCoords[0];
  if(dx > tolernance_dist)
    return 1;
  if(dx < -tolernance_dist)
    return -1;
  ergo_real dy = d1->distr.centerCoords[1] - d2->distr.centerCoords[1];
  if(dy > tolernance_dist)
    return 1;
  if(dy < -tolernance_dist)
    return -1;
  ergo_real dz = d1->distr.centerCoords[2] - d2->distr.centerCoords[2];
  if(dz > tolernance_dist)
    return 1;
  if(dz < -tolernance_dist)
    return -1;
  ergo_real de = d1->distr.exponent - d2->distr.exponent;
  if(de > tolernance_exponent)
    return 1;
  if(de < -tolernance_exponent)
    return -1; 
  return 0;
}

template <typename DistributionSpecStructType>
bool
compare_distrs_bool(const DistributionSpecStructType & p1, const DistributionSpecStructType & p2) {
  int i = compare_distrs<DistributionSpecStructType>(&p1, &p2);
  return (i == 1);
}

template <typename DistributionSpecStructType>
int
sort_distr_list(DistributionSpecStructType* list, int n) {
  std::sort(&list[0], &list[n], compare_distrs_bool<DistributionSpecStructType>);
  return 0;
}

#endif
