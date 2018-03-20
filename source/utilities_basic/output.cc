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

/** @file output.cc

    @brief Functionality for writing output messages to a text file.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <time.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "output.h"
#include "utilities.h"

static int global_memory_usage_output_flag = 0; /* mem output disabled by default */

/* We choose to have output disabled by default, until enable_output() is called, to make sure output is disabled on worker processes when CHT is used. */ 
static int global_output_enabled_flag = 0; /* output disabled by default */

static int use_printf = 0; /* use printf */

/* general output routine  */
void 
do_output(int logCategory, int logArea, const char* format, ...)
{
  va_list a;
  va_start(a, format);
  if(use_printf)
    do_voutput_printf(logCategory, logArea, format, a);
  else
    do_voutput(logCategory, logArea, format, a);
  va_end(a);
}

/* output in stdout using printf */
void 
enable_printf_output()
{
  global_output_enabled_flag = 1;
  use_printf = 1;

  //global_memory_usage_output_flag = 1;
}


int
do_voutput_printf(int logCategory, int logArea, const char* format, va_list a)
{
  if(global_output_enabled_flag == 0)
    return 0; // Output is disabled; do nothing, just return.
  
  char ss[8888]; /* FIXME: Do something nicer here. */
  int r;

  memset(ss, 0, sizeof(ss));
  r = vsnprintf(ss, sizeof(ss), format, a);

  printf("%s\n", ss);

  return r;
}


int
do_voutput(int logCategory, int logArea, const char* format, va_list a)
{
  if(global_output_enabled_flag == 0)
    return 0; // Output is disabled; do nothing, just return.
  
  char ss[8888]; /* FIXME: Do something nicer here. */
  int r;

  memset(ss, 0, sizeof(ss));
  
  /* The output line should begin with two characters 
     specifying the log category, followed by two characters
     specifying the log area. */
  switch(logCategory)
    {
    case LOG_CAT_UNDEFINED: strcat(ss, "  "); break;
    case LOG_CAT_ERROR:     strcat(ss, "ER"); break;
    case LOG_CAT_WARNING:   strcat(ss, "WA"); break;
    case LOG_CAT_INFO:      strcat(ss, "IN"); break;
    case LOG_CAT_EXTRAINFO: strcat(ss, "EX"); break;
    case LOG_CAT_RESULTS:   strcat(ss, "RE"); break;
    case LOG_CAT_TIMINGS:   strcat(ss, "TI"); break;
    case LOG_CAT_MEMUSAGE:  strcat(ss, "ME"); break;
    default:                strcat(ss, "  "); break;
    }
  switch(logArea)
    {
    case LOG_AREA_UNDEFINED: strcat(ss, "  "); break;
    case LOG_AREA_MAIN:      strcat(ss, "MA"); break;
    case LOG_AREA_SCF:       strcat(ss, "SC"); break;
    case LOG_AREA_LR:        strcat(ss, "LR"); break;
    case LOG_AREA_INTEGRALS: strcat(ss, "IN"); break;
    case LOG_AREA_DENSFROMF: strcat(ss, "DE"); break;
    case LOG_AREA_DFT      : strcat(ss, "DF"); break;
    case LOG_AREA_LOWLEVEL:  strcat(ss, "LO"); break;
    case LOG_AREA_CI:        strcat(ss, "CI"); break;
    case LOG_AREA_ED:        strcat(ss, "ED"); break;
    case LOG_AREA_GS:        strcat(ss, "GS"); break;
    default:                 strcat(ss, "  "); break;
    }

  strcat(ss, " ");

  r = vsnprintf(ss+5, sizeof(ss)-5, format, a);

#if 0
  /* This needs to be protected with mutex. */
  if(global_output_file == NULL)
    global_output_file = fopen("ergoscf.out", "wt");
  fprintf(global_output_file, "%s\n", ss);
  fflush(global_output_file);
#else
  /* append to output file */
  {
    FILE *output_file = fopen("ergoscf.out", "at");
    fprintf(output_file, "%s\n", ss);
    fclose(output_file);
  }
#endif
  return r;
}



void
do_output_time(int logCategory, int logArea, const char* s)
{
  int len;
  char timeString[88];
  char ss[222];
  time_t rawtime;
  time(&rawtime);
  strcpy(timeString, ctime(&rawtime));
  len = (int)strlen(timeString);
  if(timeString[len-1] == '\n')
    timeString[len-1] = '\0';
  sprintf(ss, "%s %s", s, timeString);
  do_output(logCategory, logArea, ss);
}

void
enable_memory_usage_output()
{
  global_memory_usage_output_flag = 1;
}

void 
enable_output()
{
  global_output_enabled_flag = 1;
}

void output_current_memory_usage(int logArea, const char* contextString)
{
  double virt, res, virtPeak;

  if(global_memory_usage_output_flag == 0)
    return;

  if(get_memory_usage_by_procfile(&virt, &res, &virtPeak) != 0)
    {
      /* error getting memory usage */
      do_output(LOG_CAT_MEMUSAGE, logArea, "memory usage at '%s': virt ???? res ???? v peak ????",
		contextString);
    }
  else
    {
      do_output(LOG_CAT_MEMUSAGE, logArea, "memory usage at '%s': virt %6.3f G res %6.3f G v peak %6.3f G",
		contextString, virt, res, virtPeak);
    }
}


