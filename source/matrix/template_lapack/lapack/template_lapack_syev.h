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
 

#ifndef TEMPLATE_LAPACK_SYEV_HEADER
#define TEMPLATE_LAPACK_SYEV_HEADER


template<class Treal>
int template_lapack_syev(const char *jobz, const char *uplo, const integer *n, Treal *a,
	 const integer *lda, Treal *w, Treal *work, const integer *lwork, 
	integer *info)
{

  //printf("entering template_lapack_syev\n");

/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DSYEV computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric matrix A.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            orthonormal eigenvectors of the matrix A.   
            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')   
            or the upper triangle (if UPLO='U') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

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
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c_n1 = -1;
     integer c__0 = 0;
     Treal c_b17 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    Treal d__1;
    /* Local variables */
     integer inde;
     Treal anrm;
     integer imax;
     Treal rmin, rmax;
     Treal sigma;
     integer iinfo;
     logical lower, wantz;
     integer nb;
     integer iscale;
     Treal safmin;
     Treal bignum;
     integer indtau;
     integer indwrk;
     integer llwork;
     Treal smlnum;
     integer lwkopt;
     logical lquery;
     Treal eps;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --w;
    --work;

    /* Initialization added by Elias to get rid of compiler warnings. */
    lwkopt = 0;
    /* Function Body */
    wantz = template_blas_lsame(jobz, "V");
    lower = template_blas_lsame(uplo, "L");
    lquery = *lwork == -1;

    *info = 0;
    if (! (wantz || template_blas_lsame(jobz, "N"))) {
	*info = -1;
    } else if (! (lower || template_blas_lsame(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3 - 1;
	if (*lwork < maxMACRO(i__1,i__2) && ! lquery) {
	    *info = -8;
	}
    }

    if (*info == 0) {
        nb = template_lapack_ilaenv(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
/* Computing MAX */
	i__1 = 1, i__2 = (nb + 2) * *n;
	lwkopt = maxMACRO(i__1,i__2);
	work[1] = (Treal) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("SYEV  ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	work[1] = 1.;
	return 0;
    }

    if (*n == 1) {
	w[1] = a_ref(1, 1);
	work[1] = 3.;
	if (wantz) {
	    a_ref(1, 1) = 1.;
	}
	return 0;
    }

/*     Get machine constants. */

    //printf("before getting machine constants.\n");

    //printf("calling template_lapack_lamch for Safe minimum\n");
    safmin = template_lapack_lamch("Safe minimum", (Treal)0);
    //printf("template_lapack_lamch for Safe minimum returned\n");

    eps = template_lapack_lamch("Precision", (Treal)0);

    //printf("after getting machine constants.\n");

    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = template_blas_sqrt(smlnum);
    rmax = template_blas_sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

    anrm = template_lapack_lansy("M", uplo, n, &a[a_offset], lda, &work[1]);
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	template_lapack_lascl(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, 
		info);
    }

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    template_lapack_sytrd(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &
	    work[indwrk], &llwork, &iinfo);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call   
       DORGTR to generate the orthogonal matrix, then call DSTEQR. */

    if (! wantz) {
	template_lapack_sterf(n, &w[1], &work[inde], info);
    } else {
	template_lapack_orgtr(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
		llwork, &iinfo);
	template_lapack_steqr(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau],
		 info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *n;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	template_blas_scal(&imax, &d__1, &w[1], &c__1);
    }

/*     Set WORK(1) to optimal workspace size. */

    work[1] = (Treal) lwkopt;

    return 0;

/*     End of DSYEV */

} /* dsyev_ */

#undef a_ref


#endif
