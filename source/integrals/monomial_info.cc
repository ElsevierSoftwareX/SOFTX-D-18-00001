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

/** @file monomial_info.cc

    @brief Code used to organize monomials; contributions like x, y,
    z, xy, xz etc. that can be used to describe e.g. basis functions
    and basis function products.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "monomial_info.h"
#include <assert.h>
#include <stdexcept>
#include <cstdio>


void monomial_info_struct::init()
{
  // first get count
  int count = 0;
  for(int n1 = 0; n1 <= MONOMIAL_N_MAX; n1++)
    {
      for(int n1x = 0; n1x <= n1; n1x++)
	for(int n1y = 0; n1y <= n1; n1y++)
	  for(int n1z = 0; n1z <= n1; n1z++)
	    {
	      if(n1x+n1y+n1z != n1)
		continue;
	      count++;
	    } /* END FOR n1x n1y n1z */
    } /* END FOR n1 */
  
  noOfMonomialsTot = count;
  monomial_list = new monomial_struct[noOfMonomialsTot];
  
  count = 0;
  for(int n1 = 0; n1 <= MONOMIAL_N_MAX; n1++)
    {
      for(int n1x = 0; n1x <= n1; n1x++)
	for(int n1y = 0; n1y <= n1; n1y++)
	  for(int n1z = 0; n1z <= n1; n1z++)
	    {
	      if(n1x+n1y+n1z != n1)
		continue;
	      assert(count < noOfMonomialsTot);
	      monomial_list[count].ix = n1x;
	      monomial_list[count].iy = n1y;
	      monomial_list[count].iz = n1z;
	      monomial_index_list[n1x][n1y][n1z] = count;
	      count++;
	    } /* END FOR n1x n1y n1z */
      no_of_monomials_list[n1] = count;
    } /* END FOR n1 */
  assert(count == noOfMonomialsTot);
}

monomial_info_struct::monomial_info_struct() : noOfMonomialsTot(0), monomial_list(0) {
  for(int n1 = 0; n1 <= MONOMIAL_N_MAX; n1++)
    no_of_monomials_list[n1] = 0;
  for(int n1x = 0; n1x <= MONOMIAL_N_MAX; n1x++)
    for(int n1y = 0; n1y <= MONOMIAL_N_MAX; n1y++)
      for(int n1z = 0; n1z <= MONOMIAL_N_MAX; n1z++)
	monomial_index_list[n1x][n1y][n1z] = -1;
}

monomial_info_struct::~monomial_info_struct()
{
  delete []monomial_list;
}

/** Function needed for Chunks&Tasks usage. */
monomial_info_struct::monomial_info_struct(const monomial_info_struct & other) 
  : noOfMonomialsTot(other.noOfMonomialsTot)
{  
  monomial_list = new monomial_struct[noOfMonomialsTot];
  memcpy(monomial_list, other.monomial_list, noOfMonomialsTot*sizeof(monomial_struct));
  memcpy(no_of_monomials_list, other.no_of_monomials_list, sizeof(no_of_monomials_list));
  memcpy(monomial_index_list, other.monomial_index_list, sizeof(monomial_index_list));
}

/** Function needed for Chunks&Tasks usage. */
void monomial_info_struct::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error in monomial_info_struct::write_to_buffer: bufferSize too small.");
  // noOfMonomialsTot
  memcpy(p, &noOfMonomialsTot, sizeof(int));
  p += sizeof(int);
  // monomial_list
  memcpy(p, monomial_list, noOfMonomialsTot*sizeof(monomial_struct));
  p += noOfMonomialsTot*sizeof(monomial_struct);
  // no_of_monomials_list
  memcpy(p, no_of_monomials_list, sizeof(no_of_monomials_list));
  p += sizeof(no_of_monomials_list);
  // monomial_index_list
  memcpy(p, monomial_index_list, sizeof(monomial_index_list));
  p += sizeof(monomial_index_list);
  // DONE!
}

/** Function needed for Chunks&Tasks usage. */
size_t monomial_info_struct::get_size() const {
  return sizeof(int)
    + noOfMonomialsTot*sizeof(monomial_struct) 
    + sizeof(no_of_monomials_list) 
    + sizeof(monomial_index_list);
}

/** Function needed for Chunks&Tasks usage. */
void monomial_info_struct::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  // noOfMonomialsTot
  memcpy(&noOfMonomialsTot, p, sizeof(int));
  p += sizeof(int);
  // monomial_list
  monomial_list = new monomial_struct[noOfMonomialsTot];
  memcpy(monomial_list, p, noOfMonomialsTot*sizeof(monomial_struct));
  p += noOfMonomialsTot*sizeof(monomial_struct);
  // no_of_monomials_list
  memcpy(no_of_monomials_list, p, sizeof(no_of_monomials_list));
  p += sizeof(no_of_monomials_list);
  // monomial_index_list
  memcpy(monomial_index_list, p, sizeof(monomial_index_list));
  p += sizeof(monomial_index_list);
  // DONE!
  if(static_cast<size_t>(p-dataBuffer) > bufferSize)
    throw std::runtime_error("Error: (p > bufferSize).");  
}

