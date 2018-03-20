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

/** @file Allocator.h

    \brief Code for Allocator class used for memory
    allocation/deallocation in matrix library. The idea of the
    Allocator class is that it manages a list of buffers of identical
    size, keeping track of where the next available free buffer
    is. This functionality can be used for allocating submatrices in
    the matrix library, since almost all submatrices have identical
    size.

    @author Elias Rudberg
*/

#ifndef MAT_ALLOCATOR_HEADER
#define MAT_ALLOCATOR_HEADER

#include <stdexcept>

namespace mat {

template<class Treal>
class Allocator
{
 public:
  Allocator(int noOfRealsPerBuffer_,
	    int noOfBuffers_) : 
  noOfRealsPerBuffer(noOfRealsPerBuffer_),
    noOfBuffers(noOfBuffers_)
    {
      buffer = new Treal[noOfBuffers * noOfRealsPerBuffer];
      nextFreeIndexList = new int[noOfBuffers];
      // Initialize nextFreeIndexList to indicate that all slots are free.
      for(int i = 0; i < noOfBuffers-1; i++)
	nextFreeIndexList[i] = i + 1;
      nextFreeIndexList[noOfBuffers-1] = -1; // last one points to -1
      firstFreeIndex = 0;
      noOfOccupiedSlots = 0;
    }
  ~Allocator() 
    {
      delete [] buffer;
      delete [] nextFreeIndexList;
    }
  Treal* alloc() {
    if(firstFreeIndex < 0)
      throw std::runtime_error("Error in Allocator::alloc(): no free slots.");
    Treal* ptrToReturn = &buffer[firstFreeIndex*noOfRealsPerBuffer];
    int firstFreeIndex_new = nextFreeIndexList[firstFreeIndex];
    nextFreeIndexList[firstFreeIndex] = -1;
    firstFreeIndex = firstFreeIndex_new;
    noOfOccupiedSlots++;
    return ptrToReturn;
  }
  void free(Treal* ptr) {
    if(ptr < buffer || ptr >= &buffer[noOfBuffers * noOfRealsPerBuffer])
      throw std::runtime_error("Error in Allocator::free(): unknown ptr.");
    int count = ptr - buffer;
    if((count % noOfRealsPerBuffer) != 0)
      throw std::runtime_error("Error in Allocator::free(): bad ptr.");
    int bufferIdx = count / noOfRealsPerBuffer;
    if(nextFreeIndexList[bufferIdx] != -1)
      throw std::runtime_error("Error in Allocator::free(): -1 not found.");
    nextFreeIndexList[bufferIdx] = firstFreeIndex;
    firstFreeIndex = bufferIdx;
    noOfOccupiedSlots--;
  }
  bool isFull() {
    if(noOfOccupiedSlots == noOfBuffers)
      return true;
    return false;
  }
  bool isEmpty() {
    if(noOfOccupiedSlots == 0)
      return true;
    return false;
  }
  bool ownsPtr(Treal* ptr) {
    if(ptr < buffer || ptr >= &buffer[noOfBuffers * noOfRealsPerBuffer])
      return false;
    return true;
  }
  int getNoOfOccupiedSlots() {
    return noOfOccupiedSlots;
  }
 private:
  int noOfRealsPerBuffer;
  int noOfBuffers;
  Treal* buffer;
  int* nextFreeIndexList;
  int firstFreeIndex;
  int noOfOccupiedSlots;
}; // end class Allocator

} /* end namespace mat */

#endif
