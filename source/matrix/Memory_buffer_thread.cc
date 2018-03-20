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

/** @file Memory_buffer_thread.cc

    \brief Code for managing aligned memory buffers, used if SSE
    intrinsics enabled.
*/

#include "Memory_buffer_thread.h"

/* This file is only used if USE_SSE_INTRINSICS is defined. */
#ifdef USE_SSE_INTRINSICS

namespace mat{
  
  // Initialization of static members
  Memory_buffer_thread* Memory_buffer_thread::ptr_to_instance = 0;
  volatile char Memory_buffer_thread::ptr_to_instance_is_valid = 0;
  unsigned int Memory_buffer_thread::bufSize = 1000000;

  void Memory_buffer_thread::create() {
    static Memory_buffer_thread theInstance;
    ptr_to_instance = &theInstance;
    theInstance.init_buffers(1); // Default number of threads is 1. 
  }

  Memory_buffer_thread& Memory_buffer_thread::instance() {
    if (!ptr_to_instance_is_valid) {
#ifdef _OPENMP
#pragma omp critical (mem_buf_thr_instance)
#endif
      {
	if (!ptr_to_instance_is_valid) {
	  create();
	  ptr_to_instance_is_valid = 1;
	}
      }
    }
    return *ptr_to_instance;
  }

  void Memory_buffer_thread::init_buffers(unsigned int const nThreads) {
#ifdef _OPENMP
#pragma omp critical (mem_buf_thr_init_buffers)
#endif
    {
      // First delete buffers if they were allocated before
      for (unsigned int ind = 0; ind < buffers.size(); ind++)
	delete[] buffers[ind];
      // Then resize
      buffers.resize(nThreads);
      for (unsigned int ind = 0; ind < buffers.size(); ind++)
	buffers[ind] = new char[bufSize];
    }
  }

  Memory_buffer_thread::~Memory_buffer_thread() {
    for (unsigned int ind = 0; ind < buffers.size(); ind++)
      delete[] buffers[ind];    
  }


  
} // end namespace mat

#endif

