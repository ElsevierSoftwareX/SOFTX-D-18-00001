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

/** @file boys_function_test.cc Tests the Boys function evaluation
    accuracy. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include <vector>

#include "integral_info.h"
#include "matInclude.h"
#include "template_blas_common.h"
#include "config.h" // Needed to get the PRECISION_SINGLE macro

static ergo_real rand_0_to_1() {
  int randomint = rand();
  ergo_real x = (ergo_real)randomint;
  return x / RAND_MAX;
}

static ergo_real rand_1_to_10() {
  ergo_real rand_0_to_9 = rand_0_to_1() * 9;
  ergo_real result = rand_0_to_9 + 1;
  return result;
}

static ergo_real BoysFuncAccurate(int n, ergo_real x, const IntegralInfo & integralInfo) {
  int noOfIntegrationIntervals = 40;
  ergo_real machine_epsilon = mat::getMachineEpsilon<ergo_real>();
  ergo_real diffLimit = template_blas_sqrt(template_blas_get_num_limit_min<ergo_real>());
  ergo_real reldiffLimit = template_blas_pow(machine_epsilon, (ergo_real)0.75);
  const int NMAX = 30;
  ergo_real list[NMAX];
  for(int i = 0; i < NMAX; i++) {
    list[i] = integralInfo.BoysFunction_expensive(n, x, noOfIntegrationIntervals);
    assert(list[i] >= 0);
    if(i >= 2) {
      // Check if difference between different results is small
      ergo_real diff1 = template_blas_fabs(list[i-0]-list[i-1]);
      ergo_real diff2 = template_blas_fabs(list[i-1]-list[i-2]);
      if(diff1 < diffLimit && diff2 < diffLimit)
	return list[i];
      ergo_real reldiff1 = diff1 / list[i];
      ergo_real reldiff2 = diff2 / list[i];
      if(reldiff1 < reldiffLimit && reldiff2 < reldiffLimit)
	return list[i];
    }
    noOfIntegrationIntervals *= 2;
  }
  printf("ERROR in BoysFuncAccurate(): failed to get accurate result.\n");
  return 0;
}

int main(int argc, char *argv[]) {
  IntegralInfo integralInfo(true);
  integralInfo.init();
  ergo_real machine_epsilon = mat::getMachineEpsilon<ergo_real>();
  ergo_real num_limit_min = template_blas_get_num_limit_min<ergo_real>();
  ergo_real min_value_for_meaningful_reldiff = template_blas_sqrt(num_limit_min);
  printf("boys_function_test start, machine_epsilon = %g\n", (double)machine_epsilon);
  // Different things to test: try different values of n and different values of x.
  // Define constant POW_OF_10_LIMIT that determines the range of sizes of x-values that are tested.
#ifdef PRECISION_SINGLE
  const int POW_OF_10_LIMIT = 4;
#else
  const int POW_OF_10_LIMIT = 10;
#endif
  int saved_n = -1;
  ergo_real saved_x = 0;
  ergo_real saved_absdiff = 0;
  ergo_real saved_value = 0;
  ergo_real saved_refvalue = 0;
  ergo_real maxreldiff = 0;
  // Try all values of n from 0 up to BOYS_N_MAX-1
  printf("Testing Boys function accuracy for 0 <= n <= %d\n", BOYS_N_MAX-1);
  for(int n = 0; n < BOYS_N_MAX; n++) {
    // Try different order of magnitude of x-values
    for(int i = -POW_OF_10_LIMIT; i < POW_OF_10_LIMIT; i++) {
      ergo_real x_base = template_blas_pow((ergo_real)10, (ergo_real)i);
      const int nTests = 10;
      for(int k = 0; k < nTests; k++) {
	// Get random x-value around x_base
	ergo_real x = x_base * rand_1_to_10();
	// Compare computed Boys func values for that x-value
	ergo_real boysFuncResultRef = BoysFuncAccurate(n, x, integralInfo);
	ergo_real boysFuncResult = integralInfo.BoysFunction(n, x);
	ergo_real absdiff = template_blas_fabs(boysFuncResult - boysFuncResultRef);
	ergo_real reldiff = 0;
	if(template_blas_fabs(boysFuncResultRef) > min_value_for_meaningful_reldiff)
	  reldiff = absdiff / template_blas_fabs(boysFuncResultRef);
	else {
	  // boysFuncResultRef is extremely small. In this case we just check that absdiff is extremely small also
	  if(absdiff > 2*min_value_for_meaningful_reldiff)
	    reldiff = 1;
	}
	if(reldiff > maxreldiff) {
	  maxreldiff = reldiff;
	  saved_n = n;
	  saved_x = x;
	  saved_absdiff = absdiff;
	  saved_refvalue = boysFuncResultRef;
	  saved_value = boysFuncResult;
	}
      }
    }
  }
  printf("maxreldiff = %g\n", (double)maxreldiff);
  printf("saved_n = %d\n", saved_n);
  printf("saved_x = %g\n", (double)saved_x);
  printf("saved_absdiff = %g\n", (double)saved_absdiff);
  printf("saved_refvalue = %g\n", (double)saved_refvalue);
  printf("saved_value    = %g\n", (double)saved_value);
  ergo_real tolerance = template_blas_pow(machine_epsilon, (ergo_real)0.75);
  if(maxreldiff > tolerance) {
    printf("Error in boys_function_test: (maxreldiff > tolerance). tolerance= %g\n", (double)tolerance);
    return -1;
  }
  printf("boys_function_test finished OK (maxreldiff was below tolerance %g).\n", (double)tolerance);
  return 0;
}
