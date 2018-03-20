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
 

#ifndef TEMPLATE_LAPACK_SYGV_HEADER
#define TEMPLATE_LAPACK_SYGV_HEADER


template<class Treal>
int template_lapack_sygv(const integer *itype, const char *jobz, const char *uplo, const integer *
	n, Treal *a, const integer *lda, Treal *b, const integer *ldb, 
	Treal *w, Treal *work, const integer *lwork, integer *info)
{

  //printf("entering template_lapack_sygv\n");

/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DSYGV computes all the eigenvalues, and optionally, the eigenvectors   
    of a real generalized symmetric-definite eigenproblem, of the form   
    A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.   
    Here A and B are assumed to be symmetric and B is also   
    positive definite.   

    Arguments   
    =========   

    ITYPE   (input) INTEGER   
            Specifies the problem type to be solved:   
            = 1:  A*x = (lambda)*B*x   
            = 2:  A*B*x = (lambda)*x   
            = 3:  B*A*x = (lambda)*x   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangles of A and B are stored;   
            = 'L':  Lower triangles of A and B are stored.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   

            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            matrix Z of eigenvectors.  The eigenvectors are normalized   
            as follows:   
            if ITYPE = 1 or 2, Z**T*B*Z = I;   
            if ITYPE = 3, Z**T*inv(B)*Z = I.   
            If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')   
            or the lower triangle (if UPLO='L') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the symmetric positive definite matrix B.   
            If UPLO = 'U', the leading N-by-N upper triangular part of B   
            contains the upper triangular part of the matrix B.   
            If UPLO = 'L', the leading N-by-N lower triangular part of B   
            contains the lower triangular part of the matrix B.   

            On exit, if INFO <= N, the part of B containing the matrix is   
            overwritten by the triangular factor U or L from the Cholesky   
            factorization B = U**T*U or B = L*L**T.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= max(1,3*N-1).   
            For optimal efficiency, LWORK >= (NB+2)*N,   
            where NB is the blocksize for DSYTRD returned by ILAENV.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  DPOTRF or DSYEV returned an error code:   
               <= N:  if INFO = i, DSYEV failed to converge;   
                      i off-diagonal elements of an intermediate   
                      tridiagonal form did not converge to zero;   
               > N:   if INFO = N + i, for 1 <= i <= N, then the leading   
                      minor of order i of B is not positive definite.   
                      The factorization of B could not be completed and   
                      no eigenvalues or eigenvectors were computed.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c_n1 = -1;
     Treal c_b16 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    /* Local variables */
     integer neig;
     char trans[1];
     logical upper;
     logical wantz;
     integer nb;
     integer lwkopt;
     logical lquery;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    --w;
    --work;

    /* Initialization added by Elias to get rid of compiler warnings. */
    lwkopt = 0;
    /* Function Body */
    wantz = template_blas_lsame(jobz, "V");
    upper = template_blas_lsame(uplo, "U");
    lquery = *lwork == -1;

    *info = 0;
    if (*itype < 1 || *itype > 3) {
	*info = -1;
    } else if (! (wantz || template_blas_lsame(jobz, "N"))) {
	*info = -2;
    } else if (! (upper || template_blas_lsame(uplo, "L"))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -6;
    } else if (*ldb < maxMACRO(1,*n)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3 - 1;
	if (*lwork < maxMACRO(i__1,i__2) && ! lquery) {
	    *info = -11;
	}
    }

    if (*info == 0) {
	nb = template_lapack_ilaenv(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
	lwkopt = (nb + 2) * *n;
	work[1] = (Treal) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("SYGV  ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Form a Cholesky factorization of B. */

    //printf("calling template_lapack_potrf\n");
    template_lapack_potrf(uplo, n, &b[b_offset], ldb, info);
    //printf("template_lapack_potrf returned\n");


    if (*info != 0) {
	*info = *n + *info;
	return 0;
    }

/*     Transform problem to standard eigenvalue problem and solve. */

    //printf("calling template_lapack_sygst\n");
    template_lapack_sygst(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info);
    //printf("template_lapack_sygst returned\n");

    //printf("calling template_lapack_syev\n");
    template_lapack_syev(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, info);
    //printf("template_lapack_syev returned\n");

    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

	neig = *n;
	if (*info > 0) {
	    neig = *info - 1;
	}
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;   
             backtransform eigenvectors: x = inv(L)'*y or inv(U)*y */

	    if (upper) {
		*(unsigned char *)trans = 'N';
	    } else {
		*(unsigned char *)trans = 'T';
	    }

	    template_blas_trsm("Left", uplo, trans, "Non-unit", n, &neig, &c_b16, &b[
		    b_offset], ldb, &a[a_offset], lda);

	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x;   
             backtransform eigenvectors: x = L*y or U'*y */

	    if (upper) {
		*(unsigned char *)trans = 'T';
	    } else {
		*(unsigned char *)trans = 'N';
	    }

	    template_blas_trmm("Left", uplo, trans, "Non-unit", n, &neig, &c_b16, &b[
		    b_offset], ldb, &a[a_offset], lda);
	}
    }

    work[1] = (Treal) lwkopt;
    return 0;

/*     End of DSYGV */

} /* dsygv_ */

#endif
