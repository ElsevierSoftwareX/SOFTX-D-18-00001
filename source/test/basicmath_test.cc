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

/** @file basicmath_test.cc Tests some basic math functions
    such as template_blas_sqrt() template_blas_log() etc to
    see that they are working properly and that they deliver
    the expected accuracy. */

#include <stdio.h>
#include <stdlib.h>
#include <limits>

#include "realtype.h"
#include "template_blas_common.h"
#include "matInclude.h"

static ergo_real find_approx_smallest_number() {
  ergo_real x = 1;
  int count = 0;
  while(1) {
    ergo_real next = x * 0.9;
    if(!(next > 0))
      break;
    if(!(next < x))
      break;
    x = next;
    count++;
    if(count > 1000000)
      throw std::runtime_error("Error in basicmath_test.cc: max loop interations reached in find_approx_smallest_number().");
  }
  return x;
}

static ergo_real find_approx_largest_number() {
  ergo_real x = 1;
  int count = 0;
  while(1) {
    ergo_real next = x * 1.1;
    ergo_real next2 = next * 1.1;
    if(!(next2>next)) // here we check if next2 is inf
      break;
    if(!(next > 0))
      break;
    if(!(next > x))
      break;
    x = next;
    count++;
    if(count > 1000000)
      throw std::runtime_error("Error in basicmath_test.cc: max loop interations reached in find_approx_largest_number().");
  }
  return x;
}

static void verify_within_bounds(ergo_real x, ergo_real xmin, ergo_real xmax) {
  if(x >= xmin && x <= xmax)
    return;
  throw std::runtime_error("Error in basicmath_test.cc: verify_within_bounds failed. inf or nan number encountered?");
}

int main(int argc, char *argv[])
{
  int failed = 0;
  int verbose = getenv("VERBOSE") != NULL;
  ergo_real machine_epsilon = mat::getMachineEpsilon<ergo_real>();
  
  printf("machine_epsilon = %g Run with env VERBOSE for more info.\n",
         (double)machine_epsilon);

  ergo_real approx_smallest_number = find_approx_smallest_number();
  ergo_real approx_largest_number = find_approx_largest_number();
  if(verbose) {
    int ten_to_approx_what_power = -1;
    while(template_blas_pow((ergo_real)10, (ergo_real)ten_to_approx_what_power) > approx_smallest_number)
      ten_to_approx_what_power--;
    printf("approx_smallest_number is approx 10 to power %d\n", ten_to_approx_what_power);
    ten_to_approx_what_power = 1;
    while(template_blas_pow((ergo_real)10, (ergo_real)ten_to_approx_what_power) < approx_largest_number)
      ten_to_approx_what_power++;
    printf("approx_largest_number is approx 10 to power %d\n", ten_to_approx_what_power);
  }

  /* Test sqrt function for a set of random numbers. */
  ergo_real maxabsdiff_sqrt = 0;
  for(int i = 0; i < 7777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX;
      ergo_real x2 = x * x;
      ergo_real sqrt_of_x2 = template_blas_sqrt(x2);
      ergo_real diff = sqrt_of_x2 - x;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_sqrt)
	maxabsdiff_sqrt = absdiff;
    }
  if(verbose)
    printf("maxabsdiff for template_blas_sqrt: %g\n",
           (double)maxabsdiff_sqrt);
  ergo_real maxabsdiff_sqrt_requested = machine_epsilon;
  if(maxabsdiff_sqrt > maxabsdiff_sqrt_requested)
    {
      printf("ERROR: template_blas_sqrt() not accurate enough!\n");
      printf("maxabsdiff_sqrt: %g, requested: %g\n", (double)maxabsdiff_sqrt, (double)maxabsdiff_sqrt_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_sqrt() accuracy OK.\n");
  }




  /* Test exp function by computing exp(a)*exp(b) and comparing to exp(a+b) for a list of random pairs (a,b) */
  ergo_real maxabsdiff_exp = 0;
  for(int i = 0; i < 7777; i++)
    {
      ergo_real a = (ergo_real)rand() / RAND_MAX;
      ergo_real b = (ergo_real)rand() / RAND_MAX;
      ergo_real product_of_exps = template_blas_exp(a) * template_blas_exp(b);
      ergo_real exp_of_sum = template_blas_exp(a + b);
      ergo_real diff = product_of_exps - exp_of_sum;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_exp)
	maxabsdiff_exp = absdiff;
    }
  if(verbose)
    printf("maxabsdiff for template_blas_exp: %g\n", (double)maxabsdiff_exp);
  ergo_real maxabsdiff_exp_requested = machine_epsilon * 15;
  if(maxabsdiff_exp > maxabsdiff_exp_requested)
    {
      printf("ERROR: template_blas_exp() not accurate enough!\n");
      printf("maxabsdiff_exp: %g, requested: %g\n", (double)maxabsdiff_exp, (double)maxabsdiff_exp_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_exp() accuracy OK.\n");
  }





  /* Test log function by computing log(a) + log(b) and comparing to log(a*b) for a list of random pairs (a,b) */
  ergo_real maxabsdiff_log = 0;
  for(int i = 0; i < 7777; i++)
    {
      ergo_real a = (ergo_real)rand() / RAND_MAX + (ergo_real)0.1;
      ergo_real b = (ergo_real)rand() / RAND_MAX + (ergo_real)0.1;
      ergo_real sum_of_logs = template_blas_log(a) + template_blas_log(b);
      ergo_real log_of_product = template_blas_log(a * b);
      ergo_real diff = sum_of_logs - log_of_product;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_log)
	maxabsdiff_log = absdiff;
    }
  if(verbose)
    printf("maxabsdiff for template_blas_log: %g\n", (double)maxabsdiff_log);
  ergo_real maxabsdiff_log_requested = machine_epsilon * 10;
  if(maxabsdiff_log > maxabsdiff_log_requested)
    {
      printf("ERROR: template_blas_log() not accurate enough!\n");
      printf("maxabsdiff_log: %g, requested: %g\n", (double)maxabsdiff_log, (double)maxabsdiff_log_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_log() accuracy OK.\n");
  }




  /* Test log10 function by computing log10(a) + log10(b) and comparing to log10(a*b) for a list of random pairs (a,b) */
  ergo_real maxabsdiff_log10 = 0;
  for(int i = 0; i < 7777; i++)
    {
      ergo_real a = (ergo_real)rand() / RAND_MAX;
      ergo_real b = (ergo_real)rand() / RAND_MAX;
      ergo_real sum_of_log10s = template_blas_log10(a) + template_blas_log10(b);
      ergo_real log10_of_product = template_blas_log10(a * b);
      ergo_real diff = sum_of_log10s - log10_of_product;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_log10)
	maxabsdiff_log10 = absdiff;
    }
  if(verbose)
    printf("maxabsdiff for template_blas_log10: %g\n", (double)maxabsdiff_log10);
  ergo_real maxabsdiff_log10_requested = machine_epsilon * 10;
  if(maxabsdiff_log10 > maxabsdiff_log10_requested)
    {
      printf("ERROR: template_blas_log10() not accurate enough!\n");
      printf("maxabsdiff_log10: %g, requested: %g\n", (double)maxabsdiff_log10, (double)maxabsdiff_log10_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_log10() accuracy OK.\n");
  }




  /* Test erf function by comparing with a series expression */
  ergo_real piBBP = template_blas_compute_pi_BBP((ergo_real)0);
  ergo_real maxabsdiff_erf = 0;
  for(int i = 0; i < 777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX;
      ergo_real minus_1_to_pow_n = 1;
      ergo_real n_factorial = 1;
      ergo_real x_to_pow_2n_plus_1 = x;
      int n = 0;
      ergo_real sum = 0;
      while(((ergo_real)1 / n_factorial) > machine_epsilon)
	{
	  sum += (minus_1_to_pow_n / ( n_factorial * (ergo_real)( 2 * n + 1) )) * x_to_pow_2n_plus_1;
	  n++;
	  minus_1_to_pow_n *= -1;
	  n_factorial *= n;
	  x_to_pow_2n_plus_1 *= x * x;
	}
      ergo_real series_result = ((ergo_real)2 / template_blas_sqrt(piBBP)) * sum;
      ergo_real erf_value = template_blas_erf(x);
      ergo_real diff = series_result - erf_value;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_erf)
	maxabsdiff_erf = absdiff;
    } // END FOR i
  if(verbose)
    printf("maxabsdiff for template_blas_erf: %g\n", (double)maxabsdiff_erf);
  ergo_real maxabsdiff_erf_requested = machine_epsilon * 10;
  if(maxabsdiff_erf > maxabsdiff_erf_requested)
    {
      printf("ERROR: template_blas_erf() not accurate enough!\n");
      printf("maxabsdiff_erf: %g, requested: %g\n", (double)maxabsdiff_erf, (double)maxabsdiff_erf_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_erf() accuracy OK.\n");
  }




  /* Test erfc function by computing erf(x) + erfc(x) and comparing to 1 */
  ergo_real maxabsdiff_erfc = 0;
  for(int i = 0; i < 777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX;
      ergo_real erf_of_x = template_blas_erf(x);
      ergo_real erfc_of_x = template_blas_erfc(x);
      ergo_real sum = erf_of_x + erfc_of_x;
      ergo_real diff = sum - (ergo_real)1.0;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_erfc)
	maxabsdiff_erfc = absdiff;
    } // END FOR i
  if(verbose)
    printf("maxabsdiff for template_blas_erfc: %g\n", (double)maxabsdiff_erfc);
  ergo_real maxabsdiff_erfc_requested = machine_epsilon * 1;
  if(maxabsdiff_erfc > maxabsdiff_erfc_requested)
    {
      printf("ERROR: template_blas_erfc() not accurate enough!\n");
      printf("maxabsdiff_erfc: %g, requested: %g\n", (double)maxabsdiff_erfc, (double)maxabsdiff_erfc_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_erfc() accuracy OK.\n");
  }




  /* Test sin function by comparing with a series expression */
  ergo_real maxabsdiff_sin = 0;
  for(int i = 0; i < 777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX;
      ergo_real minus_1_to_pow_n = 1;
      ergo_real two_n_plus1_factorial = 1;
      ergo_real x_to_pow_2n_plus_1 = x;
      int n = 0;
      ergo_real sum = 0;
      while(((ergo_real)1 / two_n_plus1_factorial) > machine_epsilon)
	{
	  sum += (minus_1_to_pow_n / ( two_n_plus1_factorial )) * x_to_pow_2n_plus_1;
	  n++;
	  minus_1_to_pow_n *= -1;
	  two_n_plus1_factorial *= 2*n * (2*n+1);
	  x_to_pow_2n_plus_1 *= x * x;
	}
      ergo_real series_result = sum;
      ergo_real sin_value = template_blas_sin(x);
      ergo_real diff = series_result - sin_value;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_sin)
	maxabsdiff_sin = absdiff;
    } // END FOR i
  if(verbose)
    printf("maxabsdiff for template_blas_sin: %g\n", (double)maxabsdiff_sin);
  ergo_real maxabsdiff_sin_requested = machine_epsilon * 5;
  if(maxabsdiff_sin > maxabsdiff_sin_requested)
    {
      printf("ERROR: template_blas_sin() not accurate enough!\n");
      printf("maxabsdiff_sin: %g, requested: %g\n", (double)maxabsdiff_sin, (double)maxabsdiff_sin_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_sin() accuracy OK.\n");
  }




  /* Test cos function by computing cos(x) and comparing to sin(x+pi/2) */
  ergo_real maxabsdiff_cos = 0;
  for(int i = 0; i < 777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX;
      ergo_real cos_of_x = template_blas_cos(x);
      ergo_real sin_of_x_plus_pihalf = template_blas_sin(x + piBBP/2);
      ergo_real diff = cos_of_x - sin_of_x_plus_pihalf;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_cos)
	maxabsdiff_cos = absdiff;
    } // END FOR i
  if(verbose)
    printf("maxabsdiff for template_blas_cos: %g\n", (double)maxabsdiff_cos);
  ergo_real maxabsdiff_cos_requested = machine_epsilon * 3;
  if(maxabsdiff_cos > maxabsdiff_cos_requested)
    {
      printf("ERROR: template_blas_cos() not accurate enough!\n");
      printf("maxabsdiff_cos: %g, requested: %g\n", (double)maxabsdiff_cos, (double)maxabsdiff_cos_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_cos() accuracy OK.\n");
  }




  /* Test fabs function by computing x + fabs(x) for some negative numbers and comparing to 0 */
  ergo_real maxabsdiff_fabs = 0;
  for(int i = 0; i < 777; i++)
    {
      ergo_real x = (ergo_real)rand() / RAND_MAX - (ergo_real)1;
      ergo_real fabs_of_x = template_blas_fabs(x);
      ergo_real sum = x + fabs_of_x;
      ergo_real diff = sum - (ergo_real)0;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_fabs)
	maxabsdiff_fabs = absdiff;
    } // END FOR i
  if(verbose)
    printf("maxabsdiff for template_blas_fabs: %g\n", (double)maxabsdiff_fabs);
  ergo_real maxabsdiff_fabs_requested = machine_epsilon * 1;
  if(maxabsdiff_fabs > maxabsdiff_fabs_requested)
    {
      printf("ERROR: template_blas_fabs() not accurate enough!\n");
      printf("maxabsdiff_fabs: %g, requested: %g\n", (double)maxabsdiff_fabs, (double)maxabsdiff_fabs_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_fabs() accuracy OK.\n");
  }




  /* Test pow function by computing pow(a,b) and comparing to exp(b*log(a)) */
  ergo_real maxabsdiff_pow = 0;
  for(int i = 0; i < 777; i++)
    {
      ergo_real a = (ergo_real)rand() / RAND_MAX;
      ergo_real b = (ergo_real)rand() / RAND_MAX;
      ergo_real pow_ab = template_blas_pow(a, b);
      ergo_real exp_b_log_a = template_blas_exp(b*template_blas_log(a));
      ergo_real diff = pow_ab - exp_b_log_a;
      ergo_real absdiff = diff;
      if(absdiff < 0)
	absdiff *= -1;
      verify_within_bounds(absdiff, 0, approx_largest_number);
      if(absdiff > maxabsdiff_pow)
	maxabsdiff_pow = absdiff;
    } // END FOR i
  if(verbose)
    printf("maxabsdiff for template_blas_pow: %g\n", (double)maxabsdiff_pow);
  ergo_real maxabsdiff_pow_requested = machine_epsilon * 2;
  if(maxabsdiff_pow > maxabsdiff_pow_requested)
    {
      printf("ERROR: template_blas_pow() not accurate enough!\n");
      printf("maxabsdiff_pow: %g, requested: %g\n", (double)maxabsdiff_pow, (double)maxabsdiff_pow_requested);
      failed++;
    }
  else {
    if(verbose)
      printf("template_blas_pow() accuracy OK.\n");
  }




  if (!failed)
    puts("Basic math functions test succeeded.");
  else
    puts("Basic math functions test FAILED.");

  return failed;
}
