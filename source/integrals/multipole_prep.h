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

/** @file multipole_prep.h

    @brief This file contains preparatory stuff for computing
    multipole moments and related things.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef MULTIPOLE_PREP_HEADER
#define MULTIPOLE_PREP_HEADER

#include "realtype.h"
#include "polydegree.h"
#include <cstddef> /* size_t */

#define MAX_MULTIPOLE_DEGREE 15
#define MAX_NO_OF_MOMENTS_PER_MULTIPOLE ((MAX_MULTIPOLE_DEGREE+1)*(MAX_MULTIPOLE_DEGREE+1))

#define MAX_MULTIPOLE_DEGREE_BASIC BASIS_FUNC_POLY_MAX_DEGREE
#define MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC ((MAX_MULTIPOLE_DEGREE_BASIC+1)*(MAX_MULTIPOLE_DEGREE_BASIC+1))

typedef struct
{
  ergo_real centerCoords[3];
  int degree;
  int noOfMoments;
  ergo_real momentList[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
  ergo_real maxAbsMomentList[MAX_MULTIPOLE_DEGREE+1];
  ergo_real euclideanNormList[MAX_MULTIPOLE_DEGREE+1];
} multipole_struct_large;

typedef struct
{
  ergo_real centerCoords[3];
  int degree;
  int noOfMoments;
  ergo_real momentList[MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC];
} multipole_struct_small;

class MultipolePrepManager {
 public:
  typedef struct {
    int l;
    int m;
  } l_m_struct;
 private:
  int initialized_flag;
  ergo_real prepared_lm_factor_list[MAX_MULTIPOLE_DEGREE+1][MAX_MULTIPOLE_DEGREE+1];
  l_m_struct prepared_l_m_list[MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
 public:
  MultipolePrepManager();
  void init();
  bool is_initialized() const;
  const l_m_struct* get_l_m_list_ptr() const { return prepared_l_m_list; }
  ergo_real get_lm_factor(int l, int m) const;
  // Stuff needed for Chunks&Tasks usage
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};

#endif
