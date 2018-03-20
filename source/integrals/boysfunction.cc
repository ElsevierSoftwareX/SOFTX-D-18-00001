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

/** @file boysfunction.cc

    @brief Code for Boys function evaluation.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <cassert>
#include "boysfunction.h"
#include "pi.h"
#include "output.h"
#include "utilities.h"
#include "mat_gblas.h"
#include "config.h" // Needed to get the PRECISION_SINGLE macro


static ergo_real
semiFactorial(int n)
{
  assert(n >= -1);
  if(n <= 1)
    return 1;
  return n * semiFactorial(n - 2);
}


static ergo_real
BoysFunction_raw_simpson(int n, ergo_real x, int noOfIntegrationIntervals)
{
  const int N = noOfIntegrationIntervals;
  ergo_real h = (ergo_real)0.5 / N;
  ergo_real sum = 0;
  for(int k = 0; k <= 2*N; k++)
    {
      ergo_real tk = (ergo_real)k / (2*N);
      // Compute f(tk) = exp(-x*tk*tk) * pow(tk, 2*n)
      ergo_real foftk = template_blas_exp(-x*tk*tk);
      if(n != 0)
	{
	  if(k != 0)
	    foftk *= template_blas_pow(tk, (ergo_real)(2*n));
	  else
	    foftk = 0;
	}
      // OK, foftk done, now add to sum.
      if(k == 0 || k == 2*N)
	{
	  sum += foftk;
	  continue;
	}
      if(k % 2 == 1)
	{
	  sum += 4 * foftk;
	  continue;
	}
      sum += 2 * foftk;
    }
  return (h/3) * sum;
}


static ergo_real BoysFunctionIntegrand(int n, ergo_real x, ergo_real t) {
  return template_blas_exp(-x*t*t) * template_blas_pow(t, (ergo_real)(2*n));
}

/* Numerical integration using Boole's rule */
static ergo_real
BoysFunction_raw_booles_rule(int n, ergo_real x, int noOfIntegrationIntervals)
{
  const int N = noOfIntegrationIntervals;
  ergo_real intervalWidth = (ergo_real)1 / N;
  ergo_real h = intervalWidth / 4;
  ergo_real sum = 0;
  for(int k = 0; k < N; k++) {
    ergo_real x1 = (ergo_real)k / N;
    ergo_real x2 = x1 + h;
    ergo_real x3 = x1 + 2*h;
    ergo_real x4 = x1 + 3*h;
    ergo_real x5 = x1 + 4*h;
    ergo_real f_of_x1 = BoysFunctionIntegrand(n, x, x1);
    ergo_real f_of_x2 = BoysFunctionIntegrand(n, x, x2);
    ergo_real f_of_x3 = BoysFunctionIntegrand(n, x, x3);
    ergo_real f_of_x4 = BoysFunctionIntegrand(n, x, x4);
    ergo_real f_of_x5 = BoysFunctionIntegrand(n, x, x5);
    sum +=  7 * f_of_x1 +
           32 * f_of_x2 +
           12 * f_of_x3 +
           32 * f_of_x4 +
            7 * f_of_x5;
  }
  return (2*h/45)*sum;
}

/* Numerical integration using 7-point Gauss-Lobatto rule */
static ergo_real
BoysFunction_raw_GaussLobatto(int n, ergo_real x, int noOfIntegrationIntervals, ergo_real endPt = 1)
{
  if(endPt == 0)
    return 0;
  // If integrand at endPt is almost zero we get better accuracy by integrating over a shorter interval.
  if(BoysFunctionIntegrand(n, x, endPt) < template_blas_get_num_limit_min<ergo_real>())
    return BoysFunction_raw_GaussLobatto(n, x, noOfIntegrationIntervals, endPt*0.5);
  const ergo_real c_5_11 = (ergo_real)5/11;
  const ergo_real c_2_11 = (ergo_real)2/11;
  const ergo_real c_5_3  = (ergo_real)5/3 ;
  const ergo_real sqrt15 = template_blas_sqrt((ergo_real)15);
  // points xi
  ergo_real x1 = 0;
  ergo_real x2 = template_blas_sqrt(c_5_11 - c_2_11 * template_blas_sqrt(c_5_3));
  ergo_real x3 = template_blas_sqrt(c_5_11 + c_2_11 * template_blas_sqrt(c_5_3));
  ergo_real x4 = 1;
  ergo_real x5 = -x2;
  ergo_real x6 = -x3;
  ergo_real x7 = -x4;
  // weights wi
  ergo_real w1 = (ergo_real)256/525;
  ergo_real w2 = ((ergo_real)124 + (ergo_real)7*sqrt15) / 350;
  ergo_real w3 = ((ergo_real)124 - (ergo_real)7*sqrt15) / 350;
  ergo_real w4 = (ergo_real)1 / 21;
  ergo_real w5 = w2;
  ergo_real w6 = w3;
  ergo_real w7 = w4;
  const int N = noOfIntegrationIntervals;
  ergo_real intervalWidth = (ergo_real)endPt / N;
  ergo_real sum = 0;
  for(int k = 0; k < N; k++) {
    ergo_real a = (ergo_real)k*endPt / N;
    ergo_real b = a + intervalWidth;
    ergo_real b_minus_a_over_2 = (b-a)/2;
    ergo_real a_plus_b_over_2 = (a+b)/2;
    ergo_real f1 = BoysFunctionIntegrand(n, x, b_minus_a_over_2*x1+a_plus_b_over_2);
    ergo_real f2 = BoysFunctionIntegrand(n, x, b_minus_a_over_2*x2+a_plus_b_over_2);
    ergo_real f3 = BoysFunctionIntegrand(n, x, b_minus_a_over_2*x3+a_plus_b_over_2);
    ergo_real f4 = BoysFunctionIntegrand(n, x, b_minus_a_over_2*x4+a_plus_b_over_2);
    ergo_real f5 = BoysFunctionIntegrand(n, x, b_minus_a_over_2*x5+a_plus_b_over_2);
    ergo_real f6 = BoysFunctionIntegrand(n, x, b_minus_a_over_2*x6+a_plus_b_over_2);
    ergo_real f7 = BoysFunctionIntegrand(n, x, b_minus_a_over_2*x7+a_plus_b_over_2);
    sum += b_minus_a_over_2 * (w1 * f1 + w2 * f2 + w3 * f3 + w4 * f4 + w5 * f5 + w6 * f6 + w7 * f7);
  }
  return sum;
}


void
BoysFunctionManager::init(void) {
  if(Boys_init_flag == 1)
    return;
  Util::TimeMeter timeMeter;
  ergo_real halfstep, kfactorial, BoysFuncRawResult, Ak, midx;
  // Prepare Boys_list
  halfstep = (ergo_real)BOYS_X_MAX / BOYS_NO_OF_INTERVALS * 0.5;
  for(int j = 0; j < BOYS_NO_OF_INTERVALS; j++) {
    midx = (ergo_real)BOYS_X_MAX * j / BOYS_NO_OF_INTERVALS + halfstep;
    const int highest_N_needed = BOYS_N_MAX+BOYS_TAB_DEGREE-2;
    ergo_real Boys_list_curr_midx[highest_N_needed+1];
    // Use downward recursion
    Boys_list_curr_midx[highest_N_needed] = BoysFunction_expensive(highest_N_needed, midx, 160); // FIXME DO NOT USE HARD-CODED VALUE HERE?
    for(int n = highest_N_needed-1; n >= 0; n--)
      Boys_list_curr_midx[n] = (2*midx*Boys_list_curr_midx[n+1] + template_blas_exp(-midx)) / (2*n+1);
    // Now we have the Boys_list_curr_midx list so we can use that below
    for(int n = 0; n < BOYS_N_MAX; n++) {
      Boys_list[n].list[j].midx = midx;
      kfactorial = 1;
      int minusOneToPowk = 1;
      for(int k = 0; k < BOYS_TAB_DEGREE; k++) {
	BoysFuncRawResult = Boys_list_curr_midx[n+k];
	Ak = minusOneToPowk * BoysFuncRawResult / kfactorial;
	Boys_list[n].list[j].A[k] = Ak;
	kfactorial *= k+1;
	minusOneToPowk *= -1;
      } /* END FOR k */
    } /* END FOR j */
  } /* END FOR n */
  // Also prepare SavedPrefactor_list
  for(int n = 0; n < BOYS_N_MAX; n++)
    SavedPrefactor_list[n] = (semiFactorial(2*n-1) * template_blas_sqrt(pi) / template_blas_pow((ergo_real)2, (ergo_real)(n+1)));
  // Set Boys_init_flag to indicate that initialization is done
  Boys_init_flag = 1;
  timeMeter.print(LOG_AREA_INTEGRALS, "BoysFunctionManager::init");
}

ergo_real
BoysFunctionManager::BoysFunction_pretabulated(int n, ergo_real x) const
{
  const BoysFuncIntervalStruct* interval;
  ergo_real intervalWidth, count, sum, deltax, deltaxtopowk;
  int intervalIndex, k;
  if(Boys_init_flag != 1)
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: (Boys_init_flag != 1).");
  if(x < 0)
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: (x < 0).");
  if(n >= BOYS_N_MAX)
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: (n >= BOYS_N_MAX).");
  if(n < 0)
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: (n < 0).");
  if(x >= BOYS_X_MAX) {
    /* use "large x formula" */
    return SavedPrefactor_list[n] / template_blas_pow((ergo_real)x, ((ergo_real)(2*n+1))/2);
  }
  /* choose which interval to use */
  intervalWidth = (ergo_real)BOYS_X_MAX / BOYS_NO_OF_INTERVALS;
  count = x / intervalWidth;
  intervalIndex = (int)std::floor((long double)count); // FIXME: ACCURACY PROBLEM HERE FOR QUAD PRECISION?
  if((intervalIndex < 0) || (intervalIndex >= BOYS_NO_OF_INTERVALS))
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_pretabulated: bad intervalIndex.");
  interval = &Boys_list[n].list[intervalIndex];
  sum = 0;
  deltax = x - interval->midx;
  deltaxtopowk = 1;
  for(k = 0; k < BOYS_TAB_DEGREE; k++)
    {
      ergo_real Ak = interval->A[k];
      sum += Ak * deltaxtopowk;
      deltaxtopowk *= deltax;
    }
  return sum;
}

ergo_real
BoysFunctionManager::BoysFunction(int n, ergo_real x) const {
  return BoysFunction_pretabulated(n, x);
}

ergo_real BoysFunctionManager::BoysFunction_expensive(int n, ergo_real x, int noOfIntegrationIntervals, int method) const {
  if(method == 0) {
    // default case
#ifdef PRECISION_SINGLE
    return BoysFunction_raw_simpson(n, x, noOfIntegrationIntervals);
#else
    return BoysFunction_raw_GaussLobatto(n, x, noOfIntegrationIntervals);
#endif
  }
  else if(method == 1)
    return BoysFunction_raw_simpson(n, x, noOfIntegrationIntervals);
  else if(method == 2)
    return BoysFunction_raw_booles_rule(n, x, noOfIntegrationIntervals);
  else if(method == 3)
    return BoysFunction_raw_GaussLobatto(n, x, noOfIntegrationIntervals);
  else
    throw std::runtime_error("Error in BoysFunctionManager::BoysFunction_expensive: bad mthod value.");
}

BoysFunctionManager::BoysFunctionManager() : Boys_list(BOYS_N_MAX), Boys_init_flag(0) {
  for(int i = 0; i < BOYS_N_MAX; i++)
    SavedPrefactor_list[i] = 0;
}

/** Function needed for Chunks&Tasks usage. */
void BoysFunctionManager::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error in BoysFunctionManager::write_to_buffer: bufferSize too small.");
  // Boys_list
  memcpy(p, &Boys_list[0], BOYS_N_MAX*sizeof(BoysFuncIntervalSetStruct));
  p += BOYS_N_MAX*sizeof(BoysFuncIntervalSetStruct);
  // SavedPrefactor_list
  memcpy(p, SavedPrefactor_list, sizeof(SavedPrefactor_list));
  p += sizeof(SavedPrefactor_list);
  // Boys_init_flag
  memcpy(p, &Boys_init_flag, sizeof(int));
  p += sizeof(int);
  // DONE!  
}

/** Function needed for Chunks&Tasks usage. */
size_t BoysFunctionManager::get_size() const {
  return BOYS_N_MAX*sizeof(BoysFuncIntervalSetStruct) + sizeof(SavedPrefactor_list) + sizeof(int);
}

/** Function needed for Chunks&Tasks usage. */
void BoysFunctionManager::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  // Boys_list
  memcpy(&Boys_list[0], p, BOYS_N_MAX*sizeof(BoysFuncIntervalSetStruct));
  p += BOYS_N_MAX*sizeof(BoysFuncIntervalSetStruct);
  // SavedPrefactor_list
  memcpy(SavedPrefactor_list, p, sizeof(SavedPrefactor_list));
  p += sizeof(SavedPrefactor_list);
  // Boys_init_flag
  memcpy(&Boys_init_flag, p, sizeof(int));
  p += sizeof(int);
  // DONE!
  if(static_cast<size_t>(p-dataBuffer) > bufferSize)
    throw std::runtime_error("Error: (p > bufferSize).");  
}
