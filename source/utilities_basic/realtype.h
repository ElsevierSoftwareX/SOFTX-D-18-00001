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

/** @file realtype.h

    @brief Definition of the main floating-point datatype used; the
    ergo_real type.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef REALTYPEHEADER
#define REALTYPEHEADER

#include "config.h"

#ifdef PRECISION_SINGLE
typedef float ergo_real;
typedef double ergo_long_real;
#define REALTYPE_DEFINED_OK
#endif

#ifdef PRECISION_DOUBLE
typedef double ergo_real;
typedef double ergo_long_real;
#define REALTYPE_DEFINED_OK
#endif

#ifdef PRECISION_LONG_DOUBLE
typedef long double ergo_real;
typedef long double ergo_long_real;
#define REALTYPE_DEFINED_OK
#endif

#ifdef PRECISION_QUAD_FLT128
typedef __float128 ergo_real;
typedef __float128 ergo_long_real;
#define REALTYPE_DEFINED_OK
#endif

/* if precision not specified, use double as default */
#ifndef REALTYPE_DEFINED_OK
typedef double ergo_real;
typedef double ergo_long_real;
#endif

#endif
