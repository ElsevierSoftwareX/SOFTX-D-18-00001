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
 

#ifndef TEMPLATE_LAPACK_TRTI2_HEADER
#define TEMPLATE_LAPACK_TRTI2_HEADER


template<class Treal>
int template_lapack_trti2(const char *uplo, const char *diag, const integer *n, Treal *
	a, const integer *lda, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DTRTI2 computes the inverse of a real upper or lower triangular   
    matrix.   

    This is the Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular.   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the triangular matrix A.  If UPLO = 'U', the   
            leading n by n upper triangular part of the array A contains   
            the upper triangular matrix, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of the array A contains   
            the lower triangular matrix, and the strictly upper   
            triangular part of A is not referenced.  If DIAG = 'U', the   
            diagonal elements of A are also not referenced and are   
            assumed to be 1.   

            On exit, the (triangular) inverse of the original matrix, in   
            the same storage format.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
     integer j;
     logical upper;
     logical nounit;
     Treal ajj;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

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
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("TRTI2 ", &i__1);
	return 0;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (nounit) {
		a_ref(j, j) = 1. / a_ref(j, j);
		ajj = -a_ref(j, j);
	    } else {
		ajj = -1.;
	    }

/*           Compute elements 1:j-1 of j-th column. */

	    i__2 = j - 1;
	    template_blas_trmv("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &
		    a_ref(1, j), &c__1);
	    i__2 = j - 1;
	    template_blas_scal(&i__2, &ajj, &a_ref(1, j), &c__1);
/* L10: */
	}
    } else {

/*        Compute inverse of lower triangular matrix. */

	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		a_ref(j, j) = 1. / a_ref(j, j);
		ajj = -a_ref(j, j);
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

		i__1 = *n - j;
		template_blas_trmv("Lower", "No transpose", diag, &i__1, &a_ref(j + 1, j 
			+ 1), lda, &a_ref(j + 1, j), &c__1);
		i__1 = *n - j;
		template_blas_scal(&i__1, &ajj, &a_ref(j + 1, j), &c__1);
	    }
/* L20: */
	}
    }

    return 0;

/*     End of DTRTI2 */

} /* dtrti2_ */

#undef a_ref


#endif
