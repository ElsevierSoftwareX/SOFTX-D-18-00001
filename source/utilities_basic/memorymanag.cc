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

/** @file memorymanag.cc

    @brief Memory allocation/deallocation routines.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <pthread.h>
#include <stdlib.h>

#include "memorymanag.h"
#include "output.h"

static double globalNoOfBytesAllocated = 0;
static int globalMallocCount = 0;
static int globalFreeCount = 0;
static pthread_mutex_t globalMemStatLock = PTHREAD_MUTEX_INITIALIZER;

void* 
ergo_malloc(size_t noOfBytes)
{
  void* res = malloc(noOfBytes);
  if(!res) {
    double noOfBytesAsDouble = (double)noOfBytes;
    double noOfMegaBytes = noOfBytesAsDouble / 1000000;
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
	      "error in ergo_malloc, noOfBytes = %12.0f  ( %7.2f MegaBytes )", noOfBytesAsDouble, noOfMegaBytes);
    do_output_time(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "ERROR. ergo_malloc() calling exit()");
    exit(0);
  }
  pthread_mutex_lock(&globalMemStatLock);
  globalMallocCount++;
  globalNoOfBytesAllocated += noOfBytes;
  pthread_mutex_unlock(&globalMemStatLock);
  return res;
}

void 
ergo_free(void* p)
{
  free(p);
  pthread_mutex_lock(&globalMemStatLock);
  globalFreeCount++;
  pthread_mutex_unlock(&globalMemStatLock);

}

void 
report_memory_status()
{
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "mallocCount - freeCount = %i",
	    globalMallocCount - globalFreeCount);
}

