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


#ifndef TEMPLATE_BLAS_NUM_LIMITS_HEADER
#define TEMPLATE_BLAS_NUM_LIMITS_HEADER

#include <limits>

/* We need to include config.h to get macro PRECISION_QUAD_FLT128 */
#include "config.h"

#ifdef PRECISION_QUAD_FLT128
#include <quadmath.h>
#endif


/* template_blas_get_machine_epsilon(): function for getting the
   machine epsilon (the difference between 1 and the least value
   greater than 1 that is representable) for the given
   floating-point type.  */
template<typename Treal>
inline static Treal template_blas_get_machine_epsilon() {
  return std::numeric_limits<Treal>::epsilon();
}

#ifdef PRECISION_QUAD_FLT128
template<>
inline __float128 template_blas_get_machine_epsilon<__float128>() {
  return FLT128_EPSILON;
}
#endif


/* template_blas_get_num_limit_min(): function for getting the minimum
   positive normalized value for the given floating-point type.  */
template<typename Treal>
inline static Treal template_blas_get_num_limit_min() {
  return std::numeric_limits<Treal>::min();
}

#ifdef PRECISION_QUAD_FLT128
template<>
inline __float128 template_blas_get_num_limit_min<__float128>() {
  return FLT128_MIN; // FLT128_MIN: smallest positive number with full precision
}
#endif


/* template_blas_get_num_limit_max(): function for getting the maximum
   finite value for the given floating-point type.  */
template<typename Treal>
inline static Treal template_blas_get_num_limit_max() {
  return std::numeric_limits<Treal>::max();
}

#ifdef PRECISION_QUAD_FLT128
template<>
inline __float128 template_blas_get_num_limit_max<__float128>() {
  return FLT128_MAX; // FLT128_MAX: largest finite number
}
#endif


#endif
