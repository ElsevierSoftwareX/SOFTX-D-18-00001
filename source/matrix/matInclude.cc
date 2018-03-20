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

/** @file matInclude.cc

    \brief Some utilities related to hierarchical matrix library.
*/

#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdexcept>
#include <cstring>
#include <cstdio>

#include "matInclude.h"

namespace mat {

#ifdef _OPENMP
  unsigned int Params::nProcs = 0;
  unsigned int Params::matrixParallelLevel = 0;
#endif


  normType getNormType(const char* normStr) {
    if ( strcmp(normStr, "eucl") == 0)
      return mat::euclNorm;
    if ( strcmp(normStr, "frob") == 0)
      return mat::frobNorm;
    if ( strcmp(normStr, "mixed") == 0)
      return mat::mixedNorm;
    throw "Error in mat::getNormType: Unknown norm type string.";
  }

  std::string getNormTypeString(normType nType) {
    switch(nType) {
    case mat::euclNorm:  return "eucl";
    case mat::frobNorm:  return "frob";
    case mat::mixedNorm: return "mixed";
    }
    throw "Error in mat::getNormTypeString: Unknown norm type.";
  }


  // Class "Time" implementation starts here!

  double Time::get_wall_seconds() {
    struct timeval tv;
    if(gettimeofday(&tv, NULL) != 0)
      throw std::runtime_error("Error in get_wall_seconds(), in gettimeofday().");
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
  }

  Time::Time():ticTime(0) { }
  
  void Time::tic() {
    ticTime = get_wall_seconds();
  }

  float Time::toc() {
    double returnValue = get_wall_seconds() - ticTime;
    return (float)returnValue;
  }


  // Class "MemUsage" implementation starts here!

  int MemUsage::getNumberFromBuffer(const char* buffer, const char* s) {
    const char* p = buffer;
    int slen = strlen(s);
    while(1) {
      if(*p == '\0')
	return -11;
      /*  now p points to the beginning of a new line. */
      if(memcmp(p, s, slen) == 0)
	{
	  int number;
	  /*  string found! */
	  /*  skip until blank or tab */
	  while(*p != ' ' && *p != '\t' && *p != '\n' && *p != '\0')
	    p++;
	  /*  skip blanks and tabs */
	  while(*p == ' ' || *p == '\t')
	    p++;
	  /*  get number */
	  number = atoi(p);
	  /*  skip until blank or tab */
	  while(*p != ' ' && *p != '\t' && *p != '\n' && *p != '\0')
	    p++;
	  /*  skip blanks and tabs */
	  while(*p == ' ' || *p == '\t')
	    p++;
	  /*  now p should point to "kB" */
	  if(memcmp(p, "kB", 2) != 0)
	    return -22;
	  return number;
	}
      /*  skip to next line */
      while(*p != '\n' && *p != '\0')
	p++;
      p++;
    }
    return -33;
  }

  void MemUsage::getMemUsage(Values & values) {
    char fileName[888];
    const size_t PROCFILESIZE = 8888;
    char buffer[PROCFILESIZE];
    values.virt = 0;
    values.res  = 0;
    values.peak = 0;
    int pid = getpid();
    sprintf(fileName, "/proc/%i/status", pid);
    memset(buffer, 0, PROCFILESIZE);
    FILE* f = fopen(fileName, "rt");
    // Elias note 2011-01-19: Earlier an exception was thown here, but now we just return with zero result, to make it work on Mac.
    if(f == NULL)
      return;
    size_t noOfBytesRead = fread(buffer, 1, PROCFILESIZE, f);
    fclose(f);
    if(noOfBytesRead <= 0)
      throw std::runtime_error("Error reading proc status file to get mem usage.");
    if(noOfBytesRead >= PROCFILESIZE)
      throw std::runtime_error("Error reading proc status file to get mem usage: (noOfBytesRead >= PROCFILESIZE).");
    int VmSize_kB = getNumberFromBuffer(buffer, "VmSize:");
    int VmRSS_kB  = getNumberFromBuffer(buffer, "VmRSS:");
    int VmPeak_kB = getNumberFromBuffer(buffer, "VmPeak:");
    if(VmSize_kB <= 0 || VmRSS_kB <= 0)
      throw std::runtime_error("error getting VmSize_kB or VmRSS_kB.");
    values.virt   = (float)VmSize_kB / 1000000;
    values.res    = (float)VmRSS_kB  / 1000000;
    if(VmPeak_kB > 0)
      values.peak = (float)VmPeak_kB / 1000000;
    else
      values.peak = 0;
  }
  
} // end namespace mat
