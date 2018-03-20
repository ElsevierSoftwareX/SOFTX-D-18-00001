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
 

#ifndef TEMPLATE_LAPACK_GESV_HEADER
#define TEMPLATE_LAPACK_GESV_HEADER


template<class Treal>
int template_lapack_gesv(const integer *n, const integer *nrhs, Treal *a, const integer 
	*lda, integer *ipiv, Treal *b, const integer *ldb, integer *info)
{
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGESV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   

    The LU decomposition with partial pivoting and row interchanges is   
    used to factor A as   
       A = P * L * U,   
    where P is a permutation matrix, L is unit lower triangular, and U is   
    upper triangular.  The factored form of A is then used to solve the   
    system of equations A * X = B.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N coefficient matrix A.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (output) INTEGER array, dimension (N)   
            The pivot indices that define the permutation matrix P;   
            row i of the matrix was interchanged with row IPIV(i).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS matrix of right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization   
                  has been completed, but the factor U is exactly   
                  singular, so the solution could not be computed.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -4;
    } else if (*ldb < maxMACRO(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("GESV  ", &i__1);
	return 0;
    }

/*     Compute the LU factorization of A. */

    template_lapack_getrf(n, n, &a[a_offset], lda, &ipiv[1], info);
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

      template_lapack_getrs("No transpose", n, nrhs, &a[a_offset], lda, &ipiv[1], &b[
						  b_offset], ldb, info);
    }
    return 0;

/*     End of DGESV */

} /* dgesv_ */

#endif
