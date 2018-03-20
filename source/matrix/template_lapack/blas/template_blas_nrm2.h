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
 

#ifndef TEMPLATE_BLAS_NRM2_HEADER
#define TEMPLATE_BLAS_NRM2_HEADER


template<class Treal>
Treal template_blas_nrm2(const integer *n, const Treal *x, const integer *incx)
{
/*        The following loop is equivalent to this call to the LAPACK   
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */
    /* System generated locals */
    integer i__1, i__2;
    Treal ret_val, d__1;
    /* Local variables */
     Treal norm, scale, absxi;
     integer ix;
     Treal ssq;
/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   
       DNRM2 := sqrt( x'*x )   
    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to DLASSQ.   
       Sven Hammarling, Nag Ltd.   
       Parameter adjustments */
    --x;
    /* Function Body */
    if (*n < 1 || *incx < 1) {
	norm = 0.;
    } else if (*n == 1) {
	norm = absMACRO(x[1]);
    } else {
	scale = 0.;
	ssq = 1.;


	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != 0.) {
		absxi = (d__1 = x[ix], absMACRO(d__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    d__1 = scale / absxi;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * template_blas_sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of DNRM2. */

} /* dnrm2_ */

#endif
