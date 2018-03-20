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
 

#ifndef TEMPLATE_LAPACK_SYGST_HEADER
#define TEMPLATE_LAPACK_SYGST_HEADER

#include "template_lapack_sygs2.h"

template<class Treal>
int template_lapack_sygst(const integer *itype, const char *uplo, const integer *n, 
	Treal *a, const integer *lda, Treal *b, const integer *ldb, integer *
	info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYGST reduces a real symmetric-definite generalized eigenproblem   
    to standard form.   

    If ITYPE = 1, the problem is A*x = lambda*B*x,   
    and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)   

    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or   
    B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.   

    B must have been previously factorized as U**T*U or L*L**T by DPOTRF.   

    Arguments   
    =========   

    ITYPE   (input) INTEGER   
            = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);   
            = 2 or 3: compute U*A*U**T or L**T*A*L.   

    UPLO    (input) CHARACTER   
            = 'U':  Upper triangle of A is stored and B is factored as   
                    U**T*U;   
            = 'L':  Lower triangle of A is stored and B is factored as   
                    L*L**T.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the transformed matrix, stored in the   
            same format as A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input) DOUBLE PRECISION array, dimension (LDB,N)   
            The triangular factor from the Cholesky factorization of B,   
            as returned by DPOTRF.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c_n1 = -1;
     Treal c_b14 = 1.;
     Treal c_b16 = -.5;
     Treal c_b19 = -1.;
     Treal c_b52 = .5;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
     integer k;
     logical upper;
     integer kb;
     integer nb;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    upper = template_blas_lsame(uplo, "U");
    if (*itype < 1 || *itype > 3) {
	*info = -1;
    } else if (! upper && ! template_blas_lsame(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -5;
    } else if (*ldb < maxMACRO(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("SYGST ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = template_lapack_ilaenv(&c__1, "DSYGST", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

	template_lapack_sygs2(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info);
    } else {

/*        Use blocked code */

	if (*itype == 1) {
	    if (upper) {

/*              Compute inv(U')*A*inv(U) */

		i__1 = *n;
		i__2 = nb;
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = minMACRO(i__3,nb);

/*                 Update the upper triangle of A(k:n,k:n) */

		    template_lapack_sygs2(itype, uplo, &kb, &a_ref(k, k), lda, &b_ref(k, k),
			     ldb, info);
		    if (k + kb <= *n) {
			i__3 = *n - k - kb + 1;
			template_blas_trsm("Left", uplo, "Transpose", "Non-unit", &kb, &
				i__3, &c_b14, &b_ref(k, k), ldb, &a_ref(k, k 
				+ kb), lda);
			i__3 = *n - k - kb + 1;
			template_blas_symm("Left", uplo, &kb, &i__3, &c_b16, &a_ref(k, k),
				 lda, &b_ref(k, k + kb), ldb, &c_b14, &a_ref(
				k, k + kb), lda);
			i__3 = *n - k - kb + 1;
			template_blas_syr2k(uplo, "Transpose", &i__3, &kb, &c_b19, &a_ref(
				k, k + kb), lda, &b_ref(k, k + kb), ldb, &
				c_b14, &a_ref(k + kb, k + kb), lda);
			i__3 = *n - k - kb + 1;
			template_blas_symm("Left", uplo, &kb, &i__3, &c_b16, &a_ref(k, k),
				 lda, &b_ref(k, k + kb), ldb, &c_b14, &a_ref(
				k, k + kb), lda);
			i__3 = *n - k - kb + 1;
			template_blas_trsm("Right", uplo, "No transpose", "Non-unit", &kb,
				 &i__3, &c_b14, &b_ref(k + kb, k + kb), ldb, &
				a_ref(k, k + kb), lda);
		    }
/* L10: */
		}
	    } else {

/*              Compute inv(L)*A*inv(L') */

		i__2 = *n;
		i__1 = nb;
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = minMACRO(i__3,nb);

/*                 Update the lower triangle of A(k:n,k:n) */

		    template_lapack_sygs2(itype, uplo, &kb, &a_ref(k, k), lda, &b_ref(k, k),
			     ldb, info);
		    if (k + kb <= *n) {
			i__3 = *n - k - kb + 1;
			template_blas_trsm("Right", uplo, "Transpose", "Non-unit", &i__3, 
				&kb, &c_b14, &b_ref(k, k), ldb, &a_ref(k + kb,
				 k), lda);
			i__3 = *n - k - kb + 1;
			template_blas_symm("Right", uplo, &i__3, &kb, &c_b16, &a_ref(k, k)
				, lda, &b_ref(k + kb, k), ldb, &c_b14, &a_ref(
				k + kb, k), lda);
			i__3 = *n - k - kb + 1;
			template_blas_syr2k(uplo, "No transpose", &i__3, &kb, &c_b19, &
				a_ref(k + kb, k), lda, &b_ref(k + kb, k), ldb,
				 &c_b14, &a_ref(k + kb, k + kb), lda);
			i__3 = *n - k - kb + 1;
			template_blas_symm("Right", uplo, &i__3, &kb, &c_b16, &a_ref(k, k)
				, lda, &b_ref(k + kb, k), ldb, &c_b14, &a_ref(
				k + kb, k), lda);
			i__3 = *n - k - kb + 1;
			template_blas_trsm("Left", uplo, "No transpose", "Non-unit", &
				i__3, &kb, &c_b14, &b_ref(k + kb, k + kb), 
				ldb, &a_ref(k + kb, k), lda);
		    }
/* L20: */
		}
	    }
	} else {
	    if (upper) {

/*              Compute U*A*U' */

		i__1 = *n;
		i__2 = nb;
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = minMACRO(i__3,nb);

/*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */

		    i__3 = k - 1;
		    template_blas_trmm("Left", uplo, "No transpose", "Non-unit", &i__3, &
			    kb, &c_b14, &b[b_offset], ldb, &a_ref(1, k), lda);
		    i__3 = k - 1;
		    template_blas_symm("Right", uplo, &i__3, &kb, &c_b52, &a_ref(k, k), 
			    lda, &b_ref(1, k), ldb, &c_b14, &a_ref(1, k), lda);
		    i__3 = k - 1;
		    template_blas_syr2k(uplo, "No transpose", &i__3, &kb, &c_b14, &a_ref(
			    1, k), lda, &b_ref(1, k), ldb, &c_b14, &a[
			    a_offset], lda);
		    i__3 = k - 1;
		    template_blas_symm("Right", uplo, &i__3, &kb, &c_b52, &a_ref(k, k), 
			    lda, &b_ref(1, k), ldb, &c_b14, &a_ref(1, k), lda);
		    i__3 = k - 1;
		    template_blas_trmm("Right", uplo, "Transpose", "Non-unit", &i__3, &kb,
			     &c_b14, &b_ref(k, k), ldb, &a_ref(1, k), lda);
		    template_lapack_sygs2(itype, uplo, &kb, &a_ref(k, k), lda, &b_ref(k, k),
			     ldb, info);
/* L30: */
		}
	    } else {

/*              Compute L'*A*L */

		i__2 = *n;
		i__1 = nb;
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = minMACRO(i__3,nb);

/*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */

		    i__3 = k - 1;
		    template_blas_trmm("Right", uplo, "No transpose", "Non-unit", &kb, &
			    i__3, &c_b14, &b[b_offset], ldb, &a_ref(k, 1), 
			    lda);
		    i__3 = k - 1;
		    template_blas_symm("Left", uplo, &kb, &i__3, &c_b52, &a_ref(k, k), 
			    lda, &b_ref(k, 1), ldb, &c_b14, &a_ref(k, 1), lda);
		    i__3 = k - 1;
		    template_blas_syr2k(uplo, "Transpose", &i__3, &kb, &c_b14, &a_ref(k, 
			    1), lda, &b_ref(k, 1), ldb, &c_b14, &a[a_offset], 
			    lda);
		    i__3 = k - 1;
		    template_blas_symm("Left", uplo, &kb, &i__3, &c_b52, &a_ref(k, k), 
			    lda, &b_ref(k, 1), ldb, &c_b14, &a_ref(k, 1), lda);
		    i__3 = k - 1;
		    template_blas_trmm("Left", uplo, "Transpose", "Non-unit", &kb, &i__3, 
			    &c_b14, &b_ref(k, k), ldb, &a_ref(k, 1), lda);
		    template_lapack_sygs2(itype, uplo, &kb, &a_ref(k, k), lda, &b_ref(k, k),
			     ldb, info);
/* L40: */
		}
	    }
	}
    }
    return 0;

/*     End of DSYGST */

} /* dsygst_ */

#undef b_ref
#undef a_ref


#endif
