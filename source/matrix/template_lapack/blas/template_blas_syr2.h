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
 

#ifndef TEMPLATE_BLAS_SYR2_HEADER
#define TEMPLATE_BLAS_SYR2_HEADER


template<class Treal>
int template_blas_syr2(const char *uplo, const integer *n, const Treal *alpha, 
	const Treal *x, const integer *incx, const Treal *y, const integer *incy, 
	Treal *a, const integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
     integer info;
     Treal temp1, temp2;
     integer i__, j;
     integer ix, iy, jx, jy, kx, ky;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    DSYR2  performs the symmetric rank 2 operation   
       A := alpha*x*y' + alpha*y*x' + A,   
    where alpha is a scalar, x and y are n element vectors and A is an n   
    by n symmetric matrix.   
    Parameters   
    ==========   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   
                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   
                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   
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
    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   
    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper   
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower   
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    /* Initialization added by Elias to get rid of compiler warnings. */
    jx = jy = kx = ky = 0;
    /* Function Body */
    info = 0;
    if (! template_blas_lsame(uplo, "U") && ! template_blas_lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < maxMACRO(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	template_blas_erbla("SYR2  ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0 || *alpha == 0.) {
	return 0;
    }
/*     Set up the start points in X and Y if the increments are not both   
       unity. */
    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */
    if (template_blas_lsame(uplo, "U")) {
/*        Form  A  when A is stored in the upper triangle. */
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) + x[i__] * temp1 + y[
				i__] * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) + x[ix] * temp1 + y[iy] 
				* temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {
/*        Form  A  when A is stored in the lower triangle. */
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) + x[i__] * temp1 + y[
				i__] * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) + x[ix] * temp1 + y[iy] 
				* temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }
    return 0;
/*     End of DSYR2 . */
} /* dsyr2_ */
#undef a_ref

#endif
