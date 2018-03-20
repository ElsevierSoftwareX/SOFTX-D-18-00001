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

/** @file serialization_tools.h

    \brief Tools to simplify serialization of e.g. std::vector
    objects, useful when writing serialization functions needed for
    Chunks and Tasks usage.

    @author: Elias Rudberg <em>responsible</em>.
*/

#ifndef SERIALIZATION_TOOLS_HEADER
#define SERIALIZATION_TOOLS_HEADER

#include <vector>

/* general std::vector serialization template code */

template <typename VectorType>
size_t std_vector_getSize(const VectorType & v) {
  // One size_t for storing the number of elements, plus the space needed for the elements themselves.
  size_t sizeOfOneElement = sizeof(typename VectorType::value_type);
  return sizeof(size_t) + v.size() * sizeOfOneElement;
}

typedef char* CharPtrType;

template <typename VectorType>
void std_vector_writeToBuffer_and_move_ptr(const VectorType & v, CharPtrType & p) {
  size_t nElements = v.size();
  memcpy(p, &nElements, sizeof(size_t));
  p += sizeof(size_t);
  size_t sizeOfOneElement = sizeof(typename VectorType::value_type);
  size_t sizeOfAllElements = v.size() * sizeOfOneElement;
  memcpy(p, &v[0], sizeOfAllElements);
  p += sizeOfAllElements;
}

typedef const char* ConstCharPtrType;

template <typename VectorType>
void std_vector_assignFromBuffer_and_move_ptr(VectorType & v, ConstCharPtrType & p, const char* bufEndPtr) {
  assert(bufEndPtr > p);
  size_t nBytesInBuffer = bufEndPtr - p;
  assert(nBytesInBuffer >= sizeof(size_t));
  size_t nElements;
  memcpy(&nElements, p, sizeof(size_t));
  p += sizeof(size_t);
  // Check that nElements has a reasonable value.
  size_t sizeOfOneElement = sizeof(typename VectorType::value_type);
  size_t sizeOfAllElements = nElements * sizeOfOneElement;
  assert(sizeOfAllElements + sizeof(size_t) <= nBytesInBuffer);
  v.resize(nElements);
  memcpy(&v[0], p, sizeOfAllElements);
  p += sizeOfAllElements;
}


#endif
