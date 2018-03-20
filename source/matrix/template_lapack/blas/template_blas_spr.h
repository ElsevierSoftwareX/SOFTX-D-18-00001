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
 

#ifndef TEMPLATE_BLAS_SPR_HEADER
#define TEMPLATE_BLAS_SPR_HEADER

#include "template_blas_common.h"

template<class Treal>
int template_blas_spr(const char *uplo, integer *n, Treal *alpha, 
	Treal *x, integer *incx, Treal *ap)
{
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
     integer info;
     Treal temp;
     integer i__, j, k;
     integer kk, ix, jx, kx;
/*  Purpose   
    =======   
    DSPR    performs the symmetric rank 1 operation   
       A := alpha*x*x' + A,   
    where alpha is a real scalar, x is an n element vector and A is an   
    n by n symmetric matrix, supplied in packed form.   
    Parameters   
    ==========   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the matrix A is supplied in the packed   
             array AP as follows:   
                UPLO = 'U' or 'u'   The upper triangular part of A is   
                                    supplied in AP.   
                UPLO = 'L' or 'l'   The lower triangular part of A is   
                                    supplied in AP.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    AP     - DOUBLE PRECISION array of DIMENSION at least   
             ( ( n*( n + 1 ) )/2 ).   
             Before entry with  UPLO = 'U' or 'u', the array AP must   
             contain the upper triangular part of the symmetric matrix   
             packed sequentially, column by column, so that AP( 1 )   
             contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )   
             and a( 2, 2 ) respectively, and so on. On exit, the array   
             AP is overwritten by the upper triangular part of the   
             updated matrix.   
             Before entry with UPLO = 'L' or 'l', the array AP must   
             contain the lower triangular part of the symmetric matrix   
             packed sequentially, column by column, so that AP( 1 )   
             contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )   
             and a( 3, 1 ) respectively, and so on. On exit, the array   
             AP is overwritten by the lower triangular part of the   
             updated matrix.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    --ap;
    --x;
    /* Initialization added by Elias to get rid of compiler warnings. */
    kx = 0;
    /* Function Body */
    info = 0;
    if (! template_blas_lsame(uplo, "U") && ! template_blas_lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    }
    if (info != 0) {
	template_blas_erbla("SPR   ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0 || *alpha == 0.) {
	return 0;
    }
/*     Set the start point in X if the increment is not unity. */
    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }
/*     Start the operations. In this version the elements of the array AP   
       are accessed sequentially with one pass through AP. */
    kk = 1;
    if (template_blas_lsame(uplo, "U")) {
/*        Form  A  when upper triangle is stored in AP. */
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0.) {
		    temp = *alpha * x[j];
		    k = kk;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			ap[k] += x[i__] * temp;
			++k;
/* L10: */
		    }
		}
		kk += j;
/* L20: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    ix = kx;
		    i__2 = kk + j - 1;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] += x[ix] * temp;
			ix += *incx;
/* L30: */
		    }
		}
		jx += *incx;
		kk += j;
/* L40: */
	    }
	}
    } else {
/*        Form  A  when lower triangle is stored in AP. */
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0.) {
		    temp = *alpha * x[j];
		    k = kk;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			ap[k] += x[i__] * temp;
			++k;
/* L50: */
		    }
		}
		kk = kk + *n - j + 1;
/* L60: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    ix = jx;
		    i__2 = kk + *n - j;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] += x[ix] * temp;
			ix += *incx;
/* L70: */
		    }
		}
		jx += *incx;
		kk = kk + *n - j + 1;
/* L80: */
	    }
	}
    }
    return 0;
/*     End of DSPR  . */
} /* dspr_ */

#endif
