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

/** @file utilities.cc

    @brief Basic OS access utilities.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>

#include "utilities.h"
#include "output.h"



static int unique_random_filename_counter = 0;

int
generate_unique_random_filename(char* result, unsigned n)
{
  char s[888];
  unique_random_filename_counter++;
  sprintf(s, "random_unique_filename_%i_%i_%i", unique_random_filename_counter, getpid(), rand());
  if(strlen(s) >= n)
    {
      memset(result, 0, n);
      return -1;
    }
  strcpy(result, s);
  return 0;
}

long int get_file_size(const char* fileName) {
  FILE * fptr = fopen ( fileName , "rb" );
  if ( fptr == NULL )
    return -1;
  fseek( fptr, (long int)0, SEEK_END );
  long int sz = ftell( fptr );
  fclose ( fptr );
  return sz;
}

void 
get_host_name(host_name_struct* result)
{
  if(result == NULL)
    return;
  memset(result, 0, sizeof(host_name_struct));
  gethostname(result->s, sizeof(result->s)-1);
}

/* FIXME use getwd() system call instead. */
void 
get_working_directory(working_directory_struct* result)
{
#if 0
  int noOfBytesRead, i;
  char fileName[888];
  char cmdString[888];
  FILE* f;
  if(result == NULL)
    return;
  memset(result, 0, sizeof(working_directory_struct));
  generate_unique_random_filename(fileName, 888);
  sprintf(cmdString, "pwd > %s", fileName);
  if (system(cmdString) != 0)
    return;
  f = fopen(fileName, "rt");
  if(f == NULL)
    return;
  noOfBytesRead = (int)fread(result->s, 1, MAX_WORKING_DIRECTORY_LEN, f);
  fclose(f);
  if(noOfBytesRead <= 0)
    return;
  if(noOfBytesRead >= MAX_WORKING_DIRECTORY_LEN)
    return;
  /*  remove newline */
  for(i = 0; i < noOfBytesRead; i++)
    {
      if(result->s[i] == '\n')
	result->s[i] = '\0';
    }
  /*  remove temp file */
  sprintf(cmdString, "rm %s", fileName);
  if (system(cmdString) != 0)
    return;
#else
  if(result == NULL)
    return;
  if (getcwd(result->s, MAX_WORKING_DIRECTORY_LEN) == 0) 
    return;
#endif
}

int
get_memory_usage_by_ps(double* virtualMemoryGigaBytes, double* residentMemoryGigaBytes)
{
  char fileName[888];
  char cmdString[888];
  char* p;
  int pid;
  const unsigned BufSize = 888;
  char buf[BufSize];
  FILE* f;
  size_t noOfBytesRead;

  if(virtualMemoryGigaBytes == NULL || residentMemoryGigaBytes == NULL)
    return -1;
  *virtualMemoryGigaBytes  = 0;
  *residentMemoryGigaBytes = 0;

  generate_unique_random_filename(fileName, 888);

  pid = getpid();

  sprintf(cmdString, "ps -p %i -o vsize,rss | tail -n 1 > %s", pid, fileName);

  if (system(cmdString) != 0)
    return -1;
  f = fopen(fileName, "rt");
  if(f == NULL)
    return -1;

  memset(buf, 0, BufSize);
  noOfBytesRead = fread(buf, 1, BufSize, f);
  fclose(f);

  if(noOfBytesRead <= 0)
    return -1;
  if(noOfBytesRead >= BufSize)
    return -1;
  
  p = buf;

  /* printf("get_memory_usage, buf = '%s'\n", buf); */

  *virtualMemoryGigaBytes = atof(p) / 1000000;
  while(*p != ' ' && *p != '\0')
    p++;
  /*  skip blanks */
  while(*p == ' ')
    p++;
  *residentMemoryGigaBytes = atof(p) / 1000000;

  /*  remove temp file */
  sprintf(cmdString, "rm %s", fileName);
  if (system(cmdString) != -1)
    return -1;

  return 0;
}

static int getNumberFromBuffer(const char* buffer, const char* s)
{
  const char* p = buffer;
  int slen = strlen(s);
  while(1)
    {
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

#define PROCFILESIZE 8888
int 
get_memory_usage_by_procfile(double* virtualMemGigaBytes,
			     double* residentMemGigaBytes,
			     double* virtualMemPeakGigaBytes)
{
  char fileName[888];
  char buffer[PROCFILESIZE];
  int pid;
  size_t noOfBytesRead;
  FILE* f;
  int VmSize_kB, VmRSS_kB, VmPeak_kB;

  if(virtualMemGigaBytes == NULL || residentMemGigaBytes == NULL || virtualMemPeakGigaBytes == NULL)
    return -1;
  *virtualMemGigaBytes  = 0;
  *residentMemGigaBytes = 0;
  *virtualMemPeakGigaBytes = 0;

  pid = getpid();
  sprintf(fileName, "/proc/%i/status", pid);
  memset(buffer, 0, PROCFILESIZE);
  f = fopen(fileName, "rt");
  if(f == NULL)
    return -1;
  
  noOfBytesRead = fread(buffer, 1, PROCFILESIZE, f);
  fclose(f);
  
  if(noOfBytesRead <= 0)
    return -1;
  if(noOfBytesRead >= PROCFILESIZE)
    return -1;
  
  VmSize_kB = getNumberFromBuffer(buffer, "VmSize:");
  VmRSS_kB  = getNumberFromBuffer(buffer, "VmRSS:");
  VmPeak_kB = getNumberFromBuffer(buffer, "VmPeak:");
  
  if(VmSize_kB <= 0 || VmRSS_kB <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, 
		"error getting VmSize_kB or VmRSS_kB. Values returned: "
		"VmSize_kB = %i, VmRSS_kB = %i", 
		VmSize_kB, VmRSS_kB);
      return -1;
    }

  *virtualMemGigaBytes     = (double)VmSize_kB / 1000000;
  *residentMemGigaBytes    = (double)VmRSS_kB  / 1000000;
  if(VmPeak_kB > 0)
    *virtualMemPeakGigaBytes = (double)VmPeak_kB / 1000000;
  else
    *virtualMemPeakGigaBytes = 0;

  return 0;
}
#undef PROCFILESIZE
