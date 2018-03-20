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

/** @file multipole_prep.cc

    @brief This file contains preparatory stuff for computing
    multipole moments and related things.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdexcept>

#include "multipole_prep.h"
#include "output.h"
#include "template_blas_common.h"

static ergo_real
slow_factorial(int n) {
  if(n == 0)
    return 1;
  return n * slow_factorial(n-1);
}

static ergo_real
get_lm_factor_slow(int l, int m) {
  return 1.0 / template_blas_sqrt(slow_factorial(l-m)*slow_factorial(l+m));
}

static void
get_l_m_from_index(int index, int* result_l, int* result_m) {
  int l = 0;
  int count = 0;
  while(count < index) {
    l++;
    count += 2 * l + 1;
  }
  int m = l + index - count;
  *result_l = l;
  *result_m = m;
}

MultipolePrepManager::MultipolePrepManager() {
  initialized_flag = 0;
  memset(prepared_lm_factor_list, 0x00, sizeof(prepared_lm_factor_list));
  memset(prepared_l_m_list, 0x00, sizeof(prepared_l_m_list));
}

void MultipolePrepManager::init() {
  if(initialized_flag)
    return;
  for(int l = 0; l <= MAX_MULTIPOLE_DEGREE; l++)
    for(int m = 0; m <= l; m++)
      prepared_lm_factor_list[l][m] = get_lm_factor_slow(l, m);
  // initialize prepared_l_m_list
  for(int A = 0; A < MAX_NO_OF_MOMENTS_PER_MULTIPOLE; A++)
    get_l_m_from_index(A, &prepared_l_m_list[A].l, &prepared_l_m_list[A].m);
  initialized_flag = 1;
}

bool MultipolePrepManager::is_initialized() const {
  if(initialized_flag == 1)
    return true;
  else
    return false;
}

ergo_real MultipolePrepManager::get_lm_factor(int l, int m) const {
  if(m >= 0)
    return prepared_lm_factor_list[l][m];
  else
    return prepared_lm_factor_list[l][-m];
}

void MultipolePrepManager::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error in MultipolePrepManager::write_to_buffer: bufferSize too small.");
  memcpy(p, this, sizeof(MultipolePrepManager));
}

size_t MultipolePrepManager::get_size() const {
  return sizeof(MultipolePrepManager);
}

void MultipolePrepManager::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  if(bufferSize < sizeof(MultipolePrepManager))
    throw std::runtime_error("Error in MultipolePrepManager::assign_from_buffer: bufferSize too small.");
  memcpy(this, p, sizeof(MultipolePrepManager));
}
