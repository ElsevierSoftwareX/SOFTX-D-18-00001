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

/** @file monomial_info.h

    @brief Code used to organize monomials; contributions like x, y,
    z, xy, xz etc. that can be used to describe e.g. basis functions
    and basis function products.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef MONOMIAL_INFO_HEADER
#define MONOMIAL_INFO_HEADER

#include <cstring>
#include "polydegree.h"

/* We need a monomial degree that is 4 times the highest
   basisfunction polynomial degree, to handle two-electron integrals. 
*/
const int MONOMIAL_N_MAX = BASIS_FUNC_POLY_MAX_DEGREE*4;

typedef struct
{
  int ix;
  int iy;
  int iz;
} monomial_struct;

struct monomial_info_struct
{
  int noOfMonomialsTot;
  monomial_struct* monomial_list;
  int no_of_monomials_list[MONOMIAL_N_MAX+1];
  int monomial_index_list[MONOMIAL_N_MAX+1][MONOMIAL_N_MAX+1][MONOMIAL_N_MAX+1];
  void init();
  monomial_info_struct();
  ~monomial_info_struct();
  // Stuff needed for Chunks&Tasks usage
  monomial_info_struct(const monomial_info_struct & other);
  void write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const;
  size_t get_size() const;
  void assign_from_buffer ( char const * dataBuffer, size_t const bufferSize);
};

#endif
