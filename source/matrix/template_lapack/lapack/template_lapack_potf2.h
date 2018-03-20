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
 

#ifndef TEMPLATE_LAPACK_POTF2_HEADER
#define TEMPLATE_LAPACK_POTF2_HEADER


template<class Treal>
int template_lapack_potf2(const char *uplo, const integer *n, Treal *a, const integer *
	lda, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DPOTF2 computes the Cholesky factorization of a real symmetric   
    positive definite matrix A.   

    The factorization has the form   
       A = U' * U ,  if UPLO = 'U', or   
       A = L  * L',  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the unblocked version of the algorithm, calling Level 2 BLAS.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored.   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            n by n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U'*U  or A = L*L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, the leading minor of order k is not   
                 positive definite, and the factorization could not be   
                 completed.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     Treal c_b10 = -1.;
     Treal c_b12 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    Treal d__1;
    /* Local variables */
     integer j;
     logical upper;
     Treal ajj;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    upper = template_blas_lsame(uplo, "U");
    if (! upper && ! template_blas_lsame(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("POTF2 ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness. */

	    i__2 = j - 1;
	    ajj = a_ref(j, j) - template_blas_dot(&i__2, &a_ref(1, j), &c__1, &a_ref(1, j)
		    , &c__1);
	    if (ajj <= 0.) {
		a_ref(j, j) = ajj;
		goto L30;
	    }
	    ajj = template_blas_sqrt(ajj);
	    a_ref(j, j) = ajj;

/*           Compute elements J+1:N of row J. */

	    if (j < *n) {
		i__2 = j - 1;
		i__3 = *n - j;
		template_blas_gemv("Transpose", &i__2, &i__3, &c_b10, &a_ref(1, j + 1), 
			lda, &a_ref(1, j), &c__1, &c_b12, &a_ref(j, j + 1), 
			lda);
		i__2 = *n - j;
		d__1 = 1. / ajj;
		template_blas_scal(&i__2, &d__1, &a_ref(j, j + 1), lda);
	    }
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

	    i__2 = j - 1;
	    ajj = a_ref(j, j) - template_blas_dot(&i__2, &a_ref(j, 1), lda, &a_ref(j, 1), 
		    lda);
	    if (ajj <= 0.) {
		a_ref(j, j) = ajj;
		goto L30;
	    }
	    ajj = template_blas_sqrt(ajj);
	    a_ref(j, j) = ajj;

/*           Compute elements J+1:N of column J. */

	    if (j < *n) {
		i__2 = *n - j;
		i__3 = j - 1;
		template_blas_gemv("No transpose", &i__2, &i__3, &c_b10, &a_ref(j + 1, 1),
			 lda, &a_ref(j, 1), lda, &c_b12, &a_ref(j + 1, j), &
			c__1);
		i__2 = *n - j;
		d__1 = 1. / ajj;
		template_blas_scal(&i__2, &d__1, &a_ref(j + 1, j), &c__1);
	    }
/* L20: */
	}
    }
    goto L40;

L30:
    *info = j;

L40:
    return 0;

/*     End of DPOTF2 */

} /* dpotf2_ */

#undef a_ref


#endif
