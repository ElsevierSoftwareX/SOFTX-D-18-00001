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

/** @file utilities.h

    @brief Basic OS access utilities.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef UTILITIES_HEADER
#define UTILITIES_HEADER


#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>


#define MAX_HOST_NAME_LEN 100

typedef struct
{
  char s[MAX_HOST_NAME_LEN];
} host_name_struct;

#define MAX_WORKING_DIRECTORY_LEN 800

typedef struct
{
  char s[MAX_WORKING_DIRECTORY_LEN];
} working_directory_struct;

void get_host_name(host_name_struct* result);

void get_working_directory(working_directory_struct* result);

int get_memory_usage_by_ps(double* virtualMemoryGigaBytes, double* residentMemoryGigaBytes);

int get_memory_usage_by_procfile(double* virtualMemGigaBytes,
					  double* residentMemGigaBytes,
					  double* virtualMemPeakGigaBytes);

int generate_unique_random_filename(char* result, unsigned n);

long int get_file_size(const char* fileName);

#include <stdexcept>
#include "output.h"
#include "realtype.h"
namespace Util {
  /** Time-measuring class. Measures the time between the
      construction of the object and the call of the print method. */
  class TimeMeter {
  private:
    double startTimeCPU_sys;
    double startTimeCPU_usr;
    double startTimeWall;
  public:
    double get_start_time_wall_seconds() const {
      return startTimeWall;
    }
    static double get_wall_seconds() {
      struct timeval tv;
      if(gettimeofday(&tv, NULL) != 0)
        throw std::runtime_error("Error in get_wall_seconds(), in gettimeofday().");
      double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
      return seconds;
    }
    static void get_current_cpu_times(double & seconds_usr, double & seconds_sys) {
      struct rusage usage;
      if(getrusage (RUSAGE_SELF, &usage) != 0)
	throw std::runtime_error("Error in get_current_cpu_times(), in getrusage().");
      seconds_usr = usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1000000;
      seconds_sys = usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec / 1000000;
    }
    TimeMeter() {
      startTimeWall = get_wall_seconds();
      get_current_cpu_times(startTimeCPU_usr, startTimeCPU_sys);
    }
    double get_elapsed_wall_seconds() {
      double endTimeWall = get_wall_seconds();
      return endTimeWall - startTimeWall;
    }
    void print(int area, const char *routine) {
      double endTimeWall = get_wall_seconds();
      double secondsTakenWall = endTimeWall - startTimeWall;
      double seconds_usr, seconds_sys;
      get_current_cpu_times(seconds_usr, seconds_sys);
      double secondsTakenCPU_usr = seconds_usr - startTimeCPU_usr;
      double secondsTakenCPU_sys = seconds_sys - startTimeCPU_sys;
      do_output(LOG_CAT_TIMINGS, area, "%s took %9.2f usr cpu s  %9.2f sys cpu s  %9.2f wall s", 
		routine, secondsTakenCPU_usr, secondsTakenCPU_sys, secondsTakenWall);
    }
    

  };
}


#endif /* UTILITIES_HEADER */
