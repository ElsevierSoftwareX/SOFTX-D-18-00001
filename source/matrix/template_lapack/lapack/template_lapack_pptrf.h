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
 

#ifndef TEMPLATE_LAPACK_PPTRF_HEADER
#define TEMPLATE_LAPACK_PPTRF_HEADER

#include "template_lapack_common.h"

template<class Treal>
int template_lapack_pptrf(const char *uplo, const integer *n, Treal *ap, integer *
	info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPPTRF computes the Cholesky factorization of a real symmetric   
    positive definite matrix A stored in packed format.   

    The factorization has the form   
       A = U**T * U,  if UPLO = 'U', or   
       A = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            On entry, the upper or lower triangle of the symmetric matrix   
            A, packed columnwise in a linear array.  The j-th column of A   
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            See below for further details.   

            On exit, if INFO = 0, the triangular factor U or L from the   
            Cholesky factorization A = U**T*U or A = L*L**T, in the same   
            storage format as A.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   

    Further Details   
    ======= =======   

    The packed storage scheme is illustrated by the following example   
    when N = 4, UPLO = 'U':   

    Two-dimensional storage of the symmetric matrix A:   

       a11 a12 a13 a14   
           a22 a23 a24   
               a33 a34     (aij = aji)   
                   a44   

    Packed storage of the upper triangle of A:   

    AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     Treal c_b16 = -1.;
    
    /* System generated locals */
    integer i__1, i__2;
    Treal d__1;
    /* Local variables */
     integer j;
     logical upper;
     integer jc, jj;
     Treal ajj;


    --ap;

    /* Function Body */
    *info = 0;
    upper = template_blas_lsame(uplo, "U");
    if (! upper && ! template_blas_lsame(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("DPPTRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

	jj = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    jc = jj + 1;
	    jj += j;

/*           Compute elements 1:J-1 of column J. */

	    if (j > 1) {
		i__2 = j - 1;
		template_blas_tpsv("Upper", "Transpose", "Non-unit", &i__2, &ap[1], &ap[
			jc], &c__1);
	    }

/*           Compute U(J,J) and test for non-positive-definiteness. */

	    i__2 = j - 1;
	    ajj = ap[jj] - template_blas_dot(&i__2, &ap[jc], &c__1, &ap[jc], &c__1);
	    if (ajj <= 0.) {
		ap[jj] = ajj;
		goto L30;
	    }
	    ap[jj] = template_blas_sqrt(ajj);
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

	jj = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

	    ajj = ap[jj];
	    if (ajj <= 0.) {
		ap[jj] = ajj;
		goto L30;
	    }
	    ajj = template_blas_sqrt(ajj);
	    ap[jj] = ajj;

/*           Compute elements J+1:N of column J and update the trailing   
             submatrix. */

	    if (j < *n) {
		i__2 = *n - j;
		d__1 = 1. / ajj;
		template_blas_scal(&i__2, &d__1, &ap[jj + 1], &c__1);
		i__2 = *n - j;
		template_blas_spr("Lower", &i__2, &c_b16, &ap[jj + 1], &c__1, &ap[jj + *n 
			- j + 1]);
		jj = jj + *n - j + 1;
	    }
/* L20: */
	}
    }
    goto L40;

L30:
    *info = j;

L40:
    return 0;

/*     End of DPPTRF */

} /* dpptrf_ */

#endif
