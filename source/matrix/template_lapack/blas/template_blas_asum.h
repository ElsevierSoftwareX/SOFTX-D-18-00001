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
 

#ifndef TEMPLATE_BLAS_ASUM_HEADER
#define TEMPLATE_BLAS_ASUM_HEADER


template<class Treal>
Treal template_blas_asum(const integer *n, const Treal *dx, const integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    Treal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;
    /* Local variables */
     integer i__, m;
     Treal dtemp;
     integer nincx, mp1;
/*     takes the sum of the absolute values.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --dx;
    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }
/*        code for increment not equal to 1 */
    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dtemp += (d__1 = dx[i__], absMACRO(d__1));
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;
/*        code for increment equal to 1   
          clean-up loop */
L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dtemp += (d__1 = dx[i__], absMACRO(d__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
	dtemp = dtemp + (d__1 = dx[i__], absMACRO(d__1)) + (d__2 = dx[i__ + 1], 
		absMACRO(d__2)) + (d__3 = dx[i__ + 2], absMACRO(d__3)) + (d__4 = dx[i__ 
		+ 3], absMACRO(d__4)) + (d__5 = dx[i__ + 4], absMACRO(d__5)) + (d__6 = 
		dx[i__ + 5], absMACRO(d__6));
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* dasum_ */

#endif
