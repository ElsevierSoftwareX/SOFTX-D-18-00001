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

/** @file organize_distrs.h

    @brief Code for organizing a given set of primitive Gaussian
    distributions (typically coming from basis function products); the
    distributions are grouped according to their location in space,
    their exponents, etc.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef ORGANIZE_DISTRS_HEADER
#define ORGANIZE_DISTRS_HEADER

#include "output.h"
#include "multipole.h"
#include "simple_sparse_mat.h"

#include <vector>


typedef struct
{
  int startIndex;
  int distrCount;
  int nmax;
  ergo_real centerCoords[3];
  ergo_real exponent;
  ergo_real maxSizeGroup;
  ergo_real maxExtentGroup;
  ergo_real maxLimitingFactorGroup;
  ergo_real maxAbsDmatElementGroup;
  ergo_real multipoleEuclNormListForK[MAX_MULTIPOLE_DEGREE_BASIC+1];
} distr_group_struct;

typedef struct
{
  int basisFuncPairIndex;
  int monomialIndex;
  ergo_real coeff;
} minimal_distr_struct;

typedef struct
{
  int nmax;
  ergo_real exponent;
  int groupStartIndex;
  int noOfGroups;
  ergo_real maxLimitingFactorForCluster;
  ergo_real multipoleEuclNormListForK[MAX_MULTIPOLE_DEGREE_BASIC+1];
} cluster_struct;

typedef struct
{
  int index_1;
  int index_2;
  int index_1_mod;
  int index_2_mod;
  int index_inbox_1;
  int index_inbox_2;
  int pairIndex;
  ergo_real dmatElement;
} basis_func_pair_struct;

#ifndef BASIS_FUNC_POLY_MAX_DEGREE
#error The constant BASIS_FUNC_POLY_MAX_DEGREE must be defined.
#endif
#if BASIS_FUNC_POLY_MAX_DEGREE<6
#define MAX_NO_OF_BASIS_FUNC_PAIRS_PER_BATCH 1000
#else
#define MAX_NO_OF_BASIS_FUNC_PAIRS_PER_BATCH 10000
#endif

typedef struct
{
  int clusterStartIndex;
  int noOfClusters;
  int noOfBasisFuncPairs;
  int basisFuncPairListIndex;
  int basisFuncForBatchsIndex;
  int basisFuncForBatchCount;
  int global_debug_id;
} batch_struct;

struct basis_func_group_info_for_box {
  int basisFuncGroupIndex;
  ergo_real max_CS_factor;
  ergo_real maxMomentVectorNormList[MAX_MULTIPOLE_DEGREE_BASIC+1];
  int maxMultipoleDegree;
};


struct distr_org_struct {
  std::vector<minimal_distr_struct> minimalDistrList;
  std::vector<distr_group_struct> groupList;
  std::vector<cluster_struct> clusterList;
  std::vector<batch_struct> batchList;
  std::vector<basis_func_pair_struct> basisFuncPairList;
  std::vector<int> basisFuncListForBatchs;
  std::vector<int> basisFuncListForBatchs_map;
  std::vector<int> basisFuncList;
  std::vector<i_j_val_struct> spMatElementList;
  std::vector<int> spMatCountList;
  std::vector<int> spMatIdxList;
  std::vector<basis_func_group_info_for_box> basisFuncGroupInfoListForK;
  struct Data {
    ergo_real maxExtent;
    ergo_real maxDistanceOutsideBox;
    int maxNoOfMonomials;
    Data();
  };
  Data data;
  // Functions needed for CHT usage
  void writeToBuffer(char* dataBuffer, size_t const bufferSize) const;
  size_t getSize() const;
  void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
};


int
organize_distributions(const IntegralInfo & integralInfo,
		       DistributionSpecStructLabeled* distrList_in, 
		       int distrCount, 
		       distr_org_struct* result,
		       const ergo_real* boxCenterCoords,
		       ergo_real boxWidth);

#endif
