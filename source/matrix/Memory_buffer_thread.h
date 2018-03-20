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

/** @file Memory_buffer_thread.h

    \brief Code for managing aligned memory buffers, used if SSE
    intrinsics enabled.
*/

#ifndef MEMORY_BUFFER_THREAD_HEADER
#define MEMORY_BUFFER_THREAD_HEADER

/* We need to include config.h to get the USE_SSE_INTRINSICS flag. */
#include "config.h"

/* This file is only used if USE_SSE_INTRINSICS is defined. */
#ifdef USE_SSE_INTRINSICS

#include <iostream> // for bad_alloc
#include <stdexcept>
#include <vector>
#include <emmintrin.h>
#ifdef __SSE3__
#include <pmmintrin.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mat{
  
  class Memory_aligned {
  protected:
    inline void * operator new (size_t s) {
      void * p = _mm_malloc (s, 16);
      if (p == NULL)
	throw std::bad_alloc();
      else
	return p;
    }
 
    inline void * operator new[] (size_t s) {
      void * p = _mm_malloc (s, 16);
      if (p == NULL)
	throw std::bad_alloc();
      else
	return p;
    }

    inline void operator delete (void * p) {
      _mm_free(p);
    }

    inline void operator delete[] (void * p) {
      _mm_free(p);
    }
  };

  class Memory_buffer_thread : public Memory_aligned {
    // Hidden stuff
  private:
    static void create();
    static Memory_buffer_thread* ptr_to_instance;
    // All values allowed for char - 
    // may be important for double checked locking pattern in instance()
    static volatile char ptr_to_instance_is_valid; 
    static unsigned int bufSize;
    Memory_buffer_thread(Memory_buffer_thread const &);
  protected:
    Memory_buffer_thread() {} // No instances of Memory_buffer_thread ever!
  
  public:
    static Memory_buffer_thread& instance();
    template<typename T>
      void get_buffer(size_t size, T* & buffer) const {
      if (sizeof(T) * size > bufSize)
	throw std::runtime_error("In Memory_buffer_thread::get_buffer : "
				 "Allocated buffer smaller than requested "
				 "buffer size");
      int threadID = 0;
#ifdef _OPENMP
      for (int ind = 0; ind <= omp_get_level(); ind++) {
	int tmp = omp_get_ancestor_thread_num(ind);
	threadID = threadID > tmp ? threadID : tmp;
      }
#endif
      if (buffers.empty())
	throw std::runtime_error("In Memory_buffer_thread::get_buffer : "
				 "buffers empty!");
      if ((unsigned)threadID >= buffers.size())
	throw std::runtime_error("In Memory_buffer_thread::get_buffer : "
				 "thread id larger than number of buffers");
      buffer = (T*)buffers[threadID];
    }
    void init_buffers(unsigned int const nThreads);
    ~Memory_buffer_thread();
  private:
    std::vector<char*> buffers;
  
  };
 

} // end namespace mat

#endif

#endif
