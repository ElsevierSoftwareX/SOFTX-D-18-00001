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
 

#ifndef TEMPLATE_LAPACK_GETRS_HEADER
#define TEMPLATE_LAPACK_GETRS_HEADER


template<class Treal>
int template_lapack_getrs(const char *trans, const integer *n, const integer *nrhs, 
	const Treal *a, const integer *lda, const integer *ipiv, Treal *b, const integer *
	ldb, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGETRS solves a system of linear equations   
       A * X = B  or  A' * X = B   
    with a general N-by-N matrix A using the LU factorization computed   
    by DGETRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A'* X = B  (Transpose)   
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The factors L and U from the factorization A = P*L*U   
            as computed by DGETRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices from DGETRF; for 1<=i<=N, row i of the   
            matrix was interchanged with row IPIV(i).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

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
     Treal c_b12 = 1.;
     integer c_n1 = -1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
     logical notran;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    notran = template_blas_lsame(trans, "N");
    if (! notran && ! template_blas_lsame(trans, "T") && ! template_blas_lsame(
	    trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -5;
    } else if (*ldb < maxMACRO(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("GETRS ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (notran) {

/*        Solve A * X = B.   

          Apply row interchanges to the right hand sides. */

	template_lapack_laswp(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);

/*        Solve L*X = B, overwriting B with X. */

	template_blas_trsm("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b12, &a[
		a_offset], lda, &b[b_offset], ldb);

/*        Solve U*X = B, overwriting B with X. */

	template_blas_trsm("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b12, &
		a[a_offset], lda, &b[b_offset], ldb);
    } else {

/*        Solve A' * X = B.   

          Solve U'*X = B, overwriting B with X. */

	template_blas_trsm("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b12, &a[
		a_offset], lda, &b[b_offset], ldb);

/*        Solve L'*X = B, overwriting B with X. */

	template_blas_trsm("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b12, &a[
		a_offset], lda, &b[b_offset], ldb);

/*        Apply row interchanges to the solution vectors. */

	template_lapack_laswp(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }

    return 0;

/*     End of DGETRS */

} /* dgetrs_ */

#endif
