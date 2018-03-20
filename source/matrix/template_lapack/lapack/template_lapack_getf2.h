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
 

#ifndef TEMPLATE_LAPACK_GETF2_HEADER
#define TEMPLATE_LAPACK_GETF2_HEADER


template<class Treal>
int template_lapack_getf2(const integer *m, const integer *n, Treal *a, const integer *
	lda, integer *ipiv, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1992   


    Purpose   
    =======   

    DGETF2 computes an LU factorization of a general m-by-n matrix A   
    using partial pivoting with row interchanges.   

    The factorization has the form   
       A = P * L * U   
    where P is a permutation matrix, L is lower triangular with unit   
    diagonal elements (lower trapezoidal if m > n), and U is upper   
    triangular (upper trapezoidal if m < n).   

    This is the right-looking Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix to be factored.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    IPIV    (output) INTEGER array, dimension (min(M,N))   
            The pivot indices; for 1 <= i <= min(M,N), row i of the   
            matrix was interchanged with row IPIV(i).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, U(k,k) is exactly zero. The factorization   
                 has been completed, but the factor U is exactly   
                 singular, and division by zero will occur if it is used   
                 to solve a system of equations.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     Treal c_b6 = -1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    Treal d__1;
    /* Local variables */
     integer j;
     integer jp;
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
	template_blas_erbla("GETF2 ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = minMACRO(*m,*n);
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

	i__2 = *m - j + 1;
	jp = j - 1 + template_blas_idamax(&i__2, &a_ref(j, j), &c__1);
	ipiv[j] = jp;
	if (a_ref(jp, j) != 0.) {

/*           Apply the interchange to columns 1:N. */

	    if (jp != j) {
	      template_blas_swap(n, &a_ref(j, 1), lda, &a_ref(jp, 1), lda);
	    }

/*           Compute elements J+1:M of J-th column. */

	    if (j < *m) {
		i__2 = *m - j;
		d__1 = 1. / a_ref(j, j);
		template_blas_scal(&i__2, &d__1, &a_ref(j + 1, j), &c__1);
	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < minMACRO(*m,*n)) {

/*           Update trailing submatrix. */

	    i__2 = *m - j;
	    i__3 = *n - j;
	    template_blas_ger(&i__2, &i__3, &c_b6, &a_ref(j + 1, j), &c__1, &a_ref(j, j + 
		    1), lda, &a_ref(j + 1, j + 1), lda);
	}
/* L10: */
    }
    return 0;

/*     End of DGETF2 */

} /* dgetf2_ */

#undef a_ref


#endif
