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

/** @file mm_limit_table.h

    @brief MMLimitTable class used to predict the magnitude of
    contributions when using truncated multipole expansions.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef MM_LIMIT_TABLE_HEADER
#define MM_LIMIT_TABLE_HEADER

#include "realtype.h"
#include "multipole_prep.h"


class MMLimitTable {
  typedef struct {
    ergo_real x[MAX_MULTIPOLE_DEGREE+1][MAX_MULTIPOLE_DEGREE_BASIC+1];
  } interaction_matrix_limit_struct;
  static const int NO_OF_STEPS_PER_RANGE = 5;
  static const int NO_OF_RANGES = 40;
  typedef struct {
    ergo_real startDistance;
    ergo_real maxDistance;
    ergo_real step;
    interaction_matrix_limit_struct list[NO_OF_STEPS_PER_RANGE];
  } interaction_matrix_limit_range_struct;
  const interaction_matrix_limit_struct & get_x_from_distance(ergo_real distance) const;
 public:
  MMLimitTable();
  ~MMLimitTable();
  void inittt(const MultipolePrepManager & multipolePrep);
  ergo_real get_max_abs_mm_contrib(int degree1,
				   const ergo_real* maxMomentVectorNormList1,
				   int degree2,
				   const ergo_real* maxMomentVectorNormList2,
				   ergo_real distance) const;
  int get_minimum_multipole_degree_needed(ergo_real distance,
					  const multipole_struct_large* boxMultipole, 
					  int maxDegreeForDistrs, 
					  const ergo_real* maxMomentVectorNormForDistrsList, 
					  ergo_real threshold) const;
  int noOfRangesUsed;
  interaction_matrix_limit_range_struct rangeList[NO_OF_RANGES];
  // Stuff needed for Chunks&Tasks usage
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};


#endif
