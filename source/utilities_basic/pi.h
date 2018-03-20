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

/** @file pi.h

    @brief Constants for the number pi and some related numbers like
    sqrt(pi).

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef PI_HEADER
#define PI_HEADER

#include "config.h" // need config.h to get PRECISION_QUAD_FLT128 macro

/* ELIAS NOTE 2016-08-30: the numbers below were produced using the
   quadmath library, by running the program in
   source/standalone/pi_quadmath.c */

#ifdef PRECISION_QUAD_FLT128

/* __float128 case: use suffix 'Q' */
#define pi         3.1415926535897932384626433832795028Q
#define sqrtpi     1.7724538509055160272981674833411451Q
#define pitopow52 17.4934183276248628462628216798715544Q

#else

/* regular, non-__float128 case: use suffix 'L' */
#define pi         3.1415926535897932384626433832795028L
#define sqrtpi     1.7724538509055160272981674833411451L
#define pitopow52 17.4934183276248628462628216798715544L

#endif

#endif
