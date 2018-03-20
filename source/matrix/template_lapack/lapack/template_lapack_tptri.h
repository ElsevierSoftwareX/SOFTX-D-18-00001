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
 

#ifndef TEMPLATE_LAPACK_TPTRI_HEADER
#define TEMPLATE_LAPACK_TPTRI_HEADER

#include "template_lapack_common.h"

template<class Treal>
int template_lapack_tptri(const char *uplo, const char *diag, const integer *n, Treal *
	ap, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DTPTRI computes the inverse of a real upper or lower triangular   
    matrix A stored in packed format.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            On entry, the upper or lower triangular matrix A, stored   
            columnwise in a linear array.  The j-th column of A is stored   
            in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n.   
            See below for further details.   
            On exit, the (triangular) inverse of the original matrix, in   
            the same packed storage format.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular   
                  matrix is singular and its inverse can not be computed.   

    Further Details   
    ===============   

    A triangular matrix A can be transferred to packed storage using one   
    of the following program segments:   

    UPLO = 'U':                      UPLO = 'L':   

          JC = 1                           JC = 1   
          DO 2 J = 1, N                    DO 2 J = 1, N   
             DO 1 I = 1, J                    DO 1 I = J, N   
                AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J)   
        1    CONTINUE                    1    CONTINUE   
             JC = JC + J                      JC = JC + N - J + 1   
        2 CONTINUE                       2 CONTINUE   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
    
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
     integer j;
     logical upper;
     integer jc, jj;
     integer jclast;
     logical nounit;
     Treal ajj;


    --ap;

    /* Initialization added by Elias to get rid of compiler warnings. */
    jclast = 0;
    /* Function Body */
    *info = 0;
    upper = template_blas_lsame(uplo, "U");
    nounit = template_blas_lsame(diag, "N");
    if (! upper && ! template_blas_lsame(uplo, "L")) {
	*info = -1;
    } else if (! nounit && ! template_blas_lsame(diag, "U")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("TPTRI ", &i__1);
	return 0;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
	if (upper) {
	    jj = 0;
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		jj += *info;
		if (ap[jj] == 0.) {
		    return 0;
		}
/* L10: */
	    }
	} else {
	    jj = 1;
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		if (ap[jj] == 0.) {
		    return 0;
		}
		jj = jj + *n - *info + 1;
/* L20: */
	    }
	}
	*info = 0;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

	jc = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (nounit) {
		ap[jc + j - 1] = 1. / ap[jc + j - 1];
		ajj = -ap[jc + j - 1];
	    } else {
		ajj = -1.;
	    }

/*           Compute elements 1:j-1 of j-th column. */

	    i__2 = j - 1;
	    template_blas_tpmv("Upper", "No transpose", diag, &i__2, &ap[1], &ap[jc], &
		    c__1);
	    i__2 = j - 1;
	    template_blas_scal(&i__2, &ajj, &ap[jc], &c__1);
	    jc += j;
/* L30: */
	}

    } else {

/*        Compute inverse of lower triangular matrix. */

	jc = *n * (*n + 1) / 2;
	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		ap[jc] = 1. / ap[jc];
		ajj = -ap[jc];
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

		i__1 = *n - j;
		template_blas_tpmv("Lower", "No transpose", diag, &i__1, &ap[jclast], &ap[
			jc + 1], &c__1);
		i__1 = *n - j;
		template_blas_scal(&i__1, &ajj, &ap[jc + 1], &c__1);
	    }
	    jclast = jc;
	    jc = jc - *n + j - 2;
/* L40: */
	}
    }

    return 0;

/*     End of DTPTRI */

} /* dtptri_ */

#endif
