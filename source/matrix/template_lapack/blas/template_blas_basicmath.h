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
 
 /* This file belongs to the template_lapack part of the Ergo source 
  * code. The source files in the template_lapack directory are modified
  * versions of files originally distributed as CLAPACK, see the
  * Copyright/license notice in the file template_lapack/COPYING.
  */
 

#ifndef TEMPLATE_BLAS_BASICMATH_HEADER
#define TEMPLATE_BLAS_BASICMATH_HEADER

#include "template_blas_num_limits.h"

template<class Treal>
Treal template_blas_fabs(Treal x);

template<class Treal>
Treal template_blas_sqrt(Treal x);

template<class Treal>
Treal template_blas_exp(Treal x);

template<class Treal>
Treal template_blas_log(Treal x);

template<class Treal>
Treal template_blas_log10(Treal x);

template<class Treal>
Treal template_blas_erf(Treal x);

template<class Treal>
Treal template_blas_erfc(Treal x);

template<class Treal>
Treal template_blas_sin(Treal x);

template<class Treal>
Treal template_blas_cos(Treal x);

template<class Treal>
Treal template_blas_pow(Treal x, Treal y);


/* template_blas_compute_pi_BBP
   This routine computes the number pi up to the precision of Treal
   using the BBP formula. */
template<class Treal>
Treal template_blas_compute_pi_BBP(Treal dummy)
{
  Treal epsilon = template_blas_get_machine_epsilon<Treal>();
  Treal one_over_16 = (Treal)1 / (Treal)16;
  Treal one_over_16_to_pow_k = 1;
  Treal sum = 0;
  int k = 0;
  do
    {
      Treal factor = 
	(Treal)4 / (Treal)(8*k + 1) - 
	(Treal)2 / (Treal)(8*k + 4) - 
	(Treal)1 / (Treal)(8*k + 5) - 
	(Treal)1 / (Treal)(8*k + 6);
      sum += one_over_16_to_pow_k * factor;
      k++;
      one_over_16_to_pow_k *= one_over_16;
    }
  while(one_over_16_to_pow_k > epsilon);
  return sum;
}


#endif
