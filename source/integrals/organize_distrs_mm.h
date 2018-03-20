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

/** @file organize_distrs_mm.h

    @brief Code for organizing a given set of primitive Gaussian
    distributions (typically coming from basis function products)
    regarding information related to multipole methods.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef ORGANIZE_DISTRS_MM_HEADER
#define ORGANIZE_DISTRS_MM_HEADER

#include "organize_distrs.h"
#include "multipole.h"
#include <vector>

struct distr_org_mm_struct {
  std::vector<multipole_struct_small> multipoleListForGroups;
  std::vector<multipole_struct_small> multipoleListForDistrs; // For CHT usage
  struct Data {
    ergo_real multipolePoint[3];
    multipole_struct_large multipole;
    ergo_real maxMomentVectorNormForDistrsList[MAX_MULTIPOLE_DEGREE_BASIC+1];
    ergo_real chargeSum;
    Data();
  };
  Data data;
  // Functions needed for CHT usage
  void writeToBuffer(char* dataBuffer, size_t const bufferSize) const;
  size_t getSize() const;
  void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
};

struct distr_list_description_struct {
  distr_org_struct org;
  distr_org_mm_struct org_mm;
  // Functions needed for CHT usage
  void writeToBuffer(char* dataBuffer, size_t const bufferSize) const;
  size_t getSize() const;
  void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
};

int
generate_multipoles_for_groups(const IntegralInfo & integralInfo,
			       const distr_org_struct & org,
			       distr_org_mm_struct & result_org_mm,
			       ergo_real* averagePosList,
			       int & avgPosCounter
			       );

int
get_multipole_pt_for_box(const ergo_real* boxCenterCoords,
			 ergo_real boxWidth,
			 const ergo_real* averagePosList,
			 int avgPosCounter,
			 ergo_real* resultMultipolePoint);

int
translate_multipoles_for_box(distr_org_mm_struct & result_org_mm,
			     const distr_org_struct & org,
			     const MMTranslator & translator
			     );

int
combine_mm_info_for_child_boxes(distr_list_description_struct & result_box_branch,
				const distr_list_description_struct** child_box_branches,
				int noOfChildren,
				const MMTranslator & translator);


#endif
