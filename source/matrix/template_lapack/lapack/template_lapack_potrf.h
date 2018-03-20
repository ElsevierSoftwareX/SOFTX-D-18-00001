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
 

#ifndef TEMPLATE_LAPACK_POTRF_HEADER
#define TEMPLATE_LAPACK_POTRF_HEADER

#include "template_lapack_potf2.h"

template<class Treal>
int template_lapack_potrf(const char *uplo, const integer *n, Treal *a, const integer *
	lda, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPOTRF computes the Cholesky factorization of a real symmetric   
    positive definite matrix A.   

    The factorization has the form   
       A = U**T * U,  if UPLO = 'U', or   
       A = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c_n1 = -1;
     Treal c_b13 = -1.;
     Treal c_b14 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
     integer j;
     logical upper;
     integer jb, nb;
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
	template_blas_erbla("POTRF ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = template_lapack_ilaenv(&c__1, "DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

	template_lapack_potf2(uplo, n, &a[a_offset], lda, info);
    } else {

/*        Use blocked code. */

	if (upper) {

/*           Compute the Cholesky factorization A = U'*U. */

	    i__1 = *n;
	    i__2 = nb;
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Update and factorize the current diagonal block and test   
                for non-positive-definiteness.   

   Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = minMACRO(i__3,i__4);
		i__3 = j - 1;
		template_blas_syrk("Upper", "Transpose", &jb, &i__3, &c_b13, &a_ref(1, j),
			 lda, &c_b14, &a_ref(j, j), lda)
			;
		template_lapack_potf2("Upper", &jb, &a_ref(j, j), lda, info);
		if (*info != 0) {
		    goto L30;
		}
		if (j + jb <= *n) {

/*                 Compute the current block row. */

		    i__3 = *n - j - jb + 1;
		    i__4 = j - 1;
		    template_blas_gemm("Transpose", "No transpose", &jb, &i__3, &i__4, &
			    c_b13, &a_ref(1, j), lda, &a_ref(1, j + jb), lda, 
			    &c_b14, &a_ref(j, j + jb), lda);
		    i__3 = *n - j - jb + 1;
		    template_blas_trsm("Left", "Upper", "Transpose", "Non-unit", &jb, &
			    i__3, &c_b14, &a_ref(j, j), lda, &a_ref(j, j + jb)
			    , lda)
			    ;
		}
/* L10: */
	    }

	} else {

/*           Compute the Cholesky factorization A = L*L'. */

	    i__2 = *n;
	    i__1 = nb;
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Update and factorize the current diagonal block and test   
                for non-positive-definiteness.   

   Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = minMACRO(i__3,i__4);
		i__3 = j - 1;
		template_blas_syrk("Lower", "No transpose", &jb, &i__3, &c_b13, &a_ref(j, 
			1), lda, &c_b14, &a_ref(j, j), lda);
		template_lapack_potf2("Lower", &jb, &a_ref(j, j), lda, info);
		if (*info != 0) {
		    goto L30;
		}
		if (j + jb <= *n) {

/*                 Compute the current block column. */

		    i__3 = *n - j - jb + 1;
		    i__4 = j - 1;
		    template_blas_gemm("No transpose", "Transpose", &i__3, &jb, &i__4, &
			    c_b13, &a_ref(j + jb, 1), lda, &a_ref(j, 1), lda, 
			    &c_b14, &a_ref(j + jb, j), lda);
		    i__3 = *n - j - jb + 1;
		    template_blas_trsm("Right", "Lower", "Transpose", "Non-unit", &i__3, &
			    jb, &c_b14, &a_ref(j, j), lda, &a_ref(j + jb, j), 
			    lda);
		}
/* L20: */
	    }
	}
    }
    goto L40;

L30:
    *info = *info + j - 1;

L40:
    return 0;

/*     End of DPOTRF */

} /* dpotrf_ */

#undef a_ref


#endif
