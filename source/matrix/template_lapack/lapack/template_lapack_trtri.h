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
 

#ifndef TEMPLATE_LAPACK_TRTRI_HEADER
#define TEMPLATE_LAPACK_TRTRI_HEADER


template<class Treal>
int template_lapack_trtri(const char *uplo, const char *diag,
			  const integer *n, Treal *a,
			  const integer *lda, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DTRTRI computes the inverse of a real upper or lower triangular   
    matrix A.   

    This is the Level 3 BLAS version of the algorithm.   

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

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the triangular matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of the array A contains   
            the upper triangular matrix, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of the array A contains   
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
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, A(i,i) is exactly zero.  The triangular   
                 matrix is singular and its inverse can not be computed.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c_n1 = -1;
     integer c__2 = 2;
     Treal c_b18 = 1.;
     Treal c_b22 = -1.;
    
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, i__1, i__2[2], i__3, i__4, i__5;
    char ch__1[2];
    /* Local variables */
     integer j;
     logical upper;
     integer jb, nb, nn;
     logical nounit;
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
	template_blas_erbla("TRTRI ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (a_ref(*info, *info) == 0.) {
		return 0;
	    }
/* L10: */
	}
	*info = 0;
    }

/*     Determine the block size for this environment.   

   Writing concatenation */
    i__2[0] = 1, a__1[0] = (char*)uplo;
    i__2[1] = 1, a__1[1] = (char*)diag;
    template_blas_s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
    nb = template_lapack_ilaenv(&c__1, "DTRTRI", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)2);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

	template_lapack_trti2(uplo, diag, n, &a[a_offset], lda, info);
    } else {

/*        Use blocked code */

	if (upper) {

/*           Compute inverse of upper triangular matrix */

	    i__1 = *n;
	    i__3 = nb;
	    for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
/* Computing MIN */
		i__4 = nb, i__5 = *n - j + 1;
		jb = minMACRO(i__4,i__5);

/*              Compute rows 1:j-1 of current block column */

		i__4 = j - 1;
		template_blas_trmm("Left", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b18, &a[a_offset], lda, &a_ref(1, j), lda);
		i__4 = j - 1;
		template_blas_trsm("Right", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b22, &a_ref(j, j), lda, &a_ref(1, j), lda);

/*              Compute inverse of current diagonal block */

		template_lapack_trti2("Upper", diag, &jb, &a_ref(j, j), lda, info);
/* L20: */
	    }
	} else {

/*           Compute inverse of lower triangular matrix */

	    nn = (*n - 1) / nb * nb + 1;
	    i__3 = -nb;
	    for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3) {
/* Computing MIN */
		i__1 = nb, i__4 = *n - j + 1;
		jb = minMACRO(i__1,i__4);
		if (j + jb <= *n) {

/*                 Compute rows j+jb:n of current block column */

		    i__1 = *n - j - jb + 1;
		    template_blas_trmm("Left", "Lower", "No transpose", diag, &i__1, &jb, 
			    &c_b18, &a_ref(j + jb, j + jb), lda, &a_ref(j + 
			    jb, j), lda);
		    i__1 = *n - j - jb + 1;
		    template_blas_trsm("Right", "Lower", "No transpose", diag, &i__1, &jb,
			     &c_b22, &a_ref(j, j), lda, &a_ref(j + jb, j), 
			    lda);
		}

/*              Compute inverse of current diagonal block */

		template_lapack_trti2("Lower", diag, &jb, &a_ref(j, j), lda, info);
/* L30: */
	    }
	}
    }

    return 0;

/*     End of DTRTRI */

} /* dtrtri_ */

#undef a_ref


#endif
