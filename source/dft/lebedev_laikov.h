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

#if !defined(_LEBEDEV_LAIKOV_H_)
#define _LEBEDEV_LAIKOV_H_ 1

/** @file lebedev_laikov.h Headers of lebedev_laikov.c. 
   Based on V.I. Lebedev, and D.N. Laikov "A quadrature formula for
       the sphere of the 131st algebraic order of accuracy" Doklady
       Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
*/  

#include "realtype.h"

typedef ergo_real real;

#if defined(__cplusplus)
extern "C" {
#endif

/**
** ll_npoint returns number of angular grid points for given L-angular
** polynomial integration accuracy.
**
** @param lvalue : grid complete through this value of angular momentum
** quantum number l.
**
** @return value : number of points in sought Lebedev-Laikov grid.
**
*/

int ll_npoint(int lvalue);


/** ll_order returns order of the smallest angular grid that has at
    least that many grid points as specified. */
int ll_order(int npoint);

/** ll_sphere fills in arrays X, Y, Z and W with the cartesian
    coordinates and weights of the grid points. 

    @param N one of the possible values returned by ll_npoint().
    @param X x cartesian coordinates of the grid points.
    @param Y y cartesian coordinates of the grid points.
    @param Z z cartesian coordinates of the grid points.
    @param W associated weights.
    @return number of actually generated points (0 for unknown value
    of N).
*/

int ll_sphere(int N, real *X, real *Y, real *Z, real *W);

#if defined(__cplusplus)
}
#endif

#endif /* _LEBEDEV_LAIKOV_H_ */
