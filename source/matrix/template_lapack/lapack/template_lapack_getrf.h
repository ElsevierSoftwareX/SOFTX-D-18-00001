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
 

#ifndef TEMPLATE_LAPACK_GETRF_HEADER
#define TEMPLATE_LAPACK_GETRF_HEADER


template<class Treal>
int template_lapack_getrf(const integer *m, const integer *n, Treal *a, const integer *
	lda, integer *ipiv, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGETRF computes an LU factorization of a general M-by-N matrix A   
    using partial pivoting with row interchanges.   

    The factorization has the form   
       A = P * L * U   
    where P is a permutation matrix, L is lower triangular with unit   
    diagonal elements (lower trapezoidal if m > n), and U is upper   
    triangular (upper trapezoidal if m < n).   

    This is the right-looking Level 3 BLAS version of the algorithm.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix to be factored.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    IPIV    (output) INTEGER array, dimension (min(M,N))   
            The pivot indices; for 1 <= i <= min(M,N), row i of the   
            matrix was interchanged with row IPIV(i).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization   
                  has been completed, but the factor U is exactly   
                  singular, and division by zero will occur if it is used   
                  to solve a system of equations.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c_n1 = -1;
     Treal c_b16 = 1.;
     Treal c_b19 = -1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
     integer i__, j;
     integer iinfo;
     integer jb, nb;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < maxMACRO(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("GETRF ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = template_lapack_ilaenv(&c__1, "DGETRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
				1);
    if (nb <= 1 || nb >= minMACRO(*m,*n)) {

/*        Use unblocked code. */

	template_lapack_getf2(m, n, &a[a_offset], lda, &ipiv[1], info);
    } else {

/*        Use blocked code. */

	i__1 = minMACRO(*m,*n);
	i__2 = nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = minMACRO(*m,*n) - j + 1;
	    jb = minMACRO(i__3,nb);

/*           Factor diagonal and subdiagonal blocks and test for exact   
             singularity. */

	    i__3 = *m - j + 1;
	    template_lapack_getf2(&i__3, &jb, &a_ref(j, j), lda, &ipiv[j], &iinfo);

/*           Adjust INFO and the pivot indices. */

	    if (*info == 0 && iinfo > 0) {
		*info = iinfo + j - 1;
	    }
/* Computing MIN */
	    i__4 = *m, i__5 = j + jb - 1;
	    i__3 = minMACRO(i__4,i__5);
	    for (i__ = j; i__ <= i__3; ++i__) {
		ipiv[i__] = j - 1 + ipiv[i__];
/* L10: */
	    }

/*           Apply interchanges to columns 1:J-1. */

	    i__3 = j - 1;
	    i__4 = j + jb - 1;
	    template_lapack_laswp(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

	    if (j + jb <= *n) {

/*              Apply interchanges to columns J+JB:N. */

		i__3 = *n - j - jb + 1;
		i__4 = j + jb - 1;
		template_lapack_laswp(&i__3, &a_ref(1, j + jb), lda, &j, &i__4, &ipiv[1], &
			c__1);

/*              Compute block row of U. */

		i__3 = *n - j - jb + 1;
		template_blas_trsm("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &
			c_b16, &a_ref(j, j), lda, &a_ref(j, j + jb), lda);
		if (j + jb <= *m) {

/*                 Update trailing submatrix. */

		    i__3 = *m - j - jb + 1;
		    i__4 = *n - j - jb + 1;
		    template_blas_gemm("No transpose", "No transpose", &i__3, &i__4, &jb, 
				       &c_b19, &a_ref(j + jb, j), lda, &a_ref(j, j + jb),
				       lda, &c_b16, &a_ref(j + jb, j + jb), lda);
		}
	    }
/* L20: */
	}
    }
    return 0;

/*     End of DGETRF */

} /* dgetrf_ */

#undef a_ref


#endif
