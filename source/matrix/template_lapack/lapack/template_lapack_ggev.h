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
 

#ifndef TEMPLATE_LAPACK_GGEV_HEADER
#define TEMPLATE_LAPACK_GGEV_HEADER


template<class Treal>
int template_lapack_ggev(const char *jobvl, const char *jobvr, const integer *n, Treal *
	a, const integer *lda, Treal *b, const integer *ldb, Treal *alphar, 
	Treal *alphai, Treal *beta, Treal *vl, const integer *ldvl, 
	Treal *vr, const integer *ldvr, Treal *work, const integer *lwork, 
	integer *info)
{
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)   
    the generalized eigenvalues, and optionally, the left and/or right   
    generalized eigenvectors.   

    A generalized eigenvalue for a pair of matrices (A,B) is a scalar   
    lambda or a ratio alpha/beta = lambda, such that A - lambda*B is   
    singular. It is usually represented as the pair (alpha,beta), as   
    there is a reasonable interpretation for beta=0, and even for both   
    being zero.   

    The right eigenvector v(j) corresponding to the eigenvalue lambda(j)   
    of (A,B) satisfies   

                     A * v(j) = lambda(j) * B * v(j).   

    The left eigenvector u(j) corresponding to the eigenvalue lambda(j)   
    of (A,B) satisfies   

                     u(j)**H * A  = lambda(j) * u(j)**H * B .   

    where u(j)**H is the conjugate-transpose of u(j).   


    Arguments   
    =========   

    JOBVL   (input) CHARACTER*1   
            = 'N':  do not compute the left generalized eigenvectors;   
            = 'V':  compute the left generalized eigenvectors.   

    JOBVR   (input) CHARACTER*1   
            = 'N':  do not compute the right generalized eigenvectors;   
            = 'V':  compute the right generalized eigenvectors.   

    N       (input) INTEGER   
            The order of the matrices A, B, VL, and VR.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the matrix A in the pair (A,B).   
            On exit, A has been overwritten.   

    LDA     (input) INTEGER   
            The leading dimension of A.  LDA >= max(1,N).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the matrix B in the pair (A,B).   
            On exit, B has been overwritten.   

    LDB     (input) INTEGER   
            The leading dimension of B.  LDB >= max(1,N).   

    ALPHAR  (output) DOUBLE PRECISION array, dimension (N)   
    ALPHAI  (output) DOUBLE PRECISION array, dimension (N)   
    BETA    (output) DOUBLE PRECISION array, dimension (N)   
            On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will   
            be the generalized eigenvalues.  If ALPHAI(j) is zero, then   
            the j-th eigenvalue is real; if positive, then the j-th and   
            (j+1)-st eigenvalues are a complex conjugate pair, with   
            ALPHAI(j+1) negative.   

            Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)   
            may easily over- or underflow, and BETA(j) may even be zero.   
            Thus, the user should avoid naively computing the ratio   
            alpha/beta.  However, ALPHAR and ALPHAI will be always less   
            than and usually comparable with norm(A) in magnitude, and   
            BETA always less than and usually comparable with norm(B).   

    VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)   
            If JOBVL = 'V', the left eigenvectors u(j) are stored one   
            after another in the columns of VL, in the same order as   
            their eigenvalues. If the j-th eigenvalue is real, then   
            u(j) = VL(:,j), the j-th column of VL. If the j-th and   
            (j+1)-th eigenvalues form a complex conjugate pair, then   
            u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).   
            Each eigenvector will be scaled so the largest component have   
            abs(real part)+abs(imag. part)=1.   
            Not referenced if JOBVL = 'N'.   

    LDVL    (input) INTEGER   
            The leading dimension of the matrix VL. LDVL >= 1, and   
            if JOBVL = 'V', LDVL >= N.   

    VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)   
            If JOBVR = 'V', the right eigenvectors v(j) are stored one   
            after another in the columns of VR, in the same order as   
            their eigenvalues. If the j-th eigenvalue is real, then   
            v(j) = VR(:,j), the j-th column of VR. If the j-th and   
            (j+1)-th eigenvalues form a complex conjugate pair, then   
            v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).   
            Each eigenvector will be scaled so the largest component have   
            abs(real part)+abs(imag. part)=1.   
            Not referenced if JOBVR = 'N'.   

    LDVR    (input) INTEGER   
            The leading dimension of the matrix VR. LDVR >= 1, and   
            if JOBVR = 'V', LDVR >= N.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,8*N).   
            For good performance, LWORK must generally be larger.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            = 1,...,N:   
                  The QZ iteration failed.  No eigenvectors have been   
                  calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)   
                  should be correct for j=INFO+1,...,N.   
            > N:  =N+1: other than QZ iteration failed in DHGEQZ.   
                  =N+2: error return from DTGEVC.   

    =====================================================================   


       Decode the input arguments   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c__0 = 0;
     Treal c_b26 = 0.;
     Treal c_b27 = 1.;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2;
    Treal d__1, d__2, d__3, d__4;
    /* Local variables */
     Treal anrm, bnrm;
     integer ierr, itau;
     Treal temp;
     logical ilvl, ilvr;
     integer iwrk;
     integer ileft, icols, irows;
     integer jc;
     integer in;
     integer jr;
     logical ilascl, ilbscl;
     logical ldumma[1];
     char chtemp[1];
     Treal bignum;
     integer ijobvl, iright, ijobvr;
     Treal anrmto, bnrmto;
     integer minwrk, maxwrk;
     Treal smlnum;
     logical lquery;
     integer ihi, ilo;
     Treal eps;
     logical ilv;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define vl_ref(a_1,a_2) vl[(a_2)*vl_dim1 + a_1]
#define vr_ref(a_1,a_2) vr[(a_2)*vr_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    --alphar;
    --alphai;
    --beta;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1 * 1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1 * 1;
    vr -= vr_offset;
    --work;

    /* Initialization added by Elias to get rid of compiler warnings. */
    maxwrk = 0;
    /* Function Body */
    if (template_blas_lsame(jobvl, "N")) {
	ijobvl = 1;
	ilvl = FALSE_;
    } else if (template_blas_lsame(jobvl, "V")) {
	ijobvl = 2;
	ilvl = TRUE_;
    } else {
	ijobvl = -1;
	ilvl = FALSE_;
    }

    if (template_blas_lsame(jobvr, "N")) {
	ijobvr = 1;
	ilvr = FALSE_;
    } else if (template_blas_lsame(jobvr, "V")) {
	ijobvr = 2;
	ilvr = TRUE_;
    } else {
	ijobvr = -1;
	ilvr = FALSE_;
    }
    ilv = ilvl || ilvr;

/*     Test the input arguments */

    *info = 0;
    lquery = *lwork == -1;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -5;
    } else if (*ldb < maxMACRO(1,*n)) {
	*info = -7;
    } else if (*ldvl < 1 || ( ilvl && *ldvl < *n ) ) {
	*info = -12;
    } else if (*ldvr < 1 || ( ilvr && *ldvr < *n ) ) {
	*info = -14;
    }

/*     Compute workspace   
        (Note: Comments in the code beginning "Workspace:" describe the   
         minimal amount of workspace needed at that point in the code,   
         as well as the preferred amount for good performance.   
         NB refers to the optimal block size for the immediately   
         following subroutine, as returned by ILAENV. The workspace is   
         computed assuming ILO = 1 and IHI = N, the worst case.) */

    minwrk = 1;
    if (*info == 0 && (*lwork >= 1 || lquery)) {
      maxwrk = *n * 7 + *n * template_lapack_ilaenv(&c__1, "DGEQRF", " ", n, &c__1, n, &
						    c__0, (ftnlen)6, (ftnlen)1);
      /* Computing MAX */
      i__1 = 1, i__2 = *n << 3;
      minwrk = maxMACRO(i__1,i__2);
      work[1] = (Treal) maxwrk;
    }

    if (*lwork < minwrk && ! lquery) {
	*info = -16;
    }

    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("GGEV  ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = template_lapack_lamch("P", (Treal)0);
    smlnum = template_lapack_lamch("S", (Treal)0);
    bignum = 1. / smlnum;
    template_lapack_labad(&smlnum, &bignum);
    smlnum = template_blas_sqrt(smlnum) / eps;
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = template_lapack_lange("M", n, n, &a[a_offset], lda, &work[1]);
    ilascl = FALSE_;
    if (anrm > 0. && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = TRUE_;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = TRUE_;
    }
    if (ilascl) {
	template_lapack_lascl("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &
		ierr);
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

    bnrm = template_lapack_lange("M", n, n, &b[b_offset], ldb, &work[1]);
    ilbscl = FALSE_;
    if (bnrm > 0. && bnrm < smlnum) {
	bnrmto = smlnum;
	ilbscl = TRUE_;
    } else if (bnrm > bignum) {
	bnrmto = bignum;
	ilbscl = TRUE_;
    }
    if (ilbscl) {
	template_lapack_lascl("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &
		ierr);
    }

/*     Permute the matrices A, B to isolate eigenvalues if possible   
       (Workspace: need 6*N) */

    ileft = 1;
    iright = *n + 1;
    iwrk = iright + *n;
    template_lapack_ggbal("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
			    ileft], &work[iright], &work[iwrk], &ierr);

/*     Reduce B to triangular form (QR decomposition of B)   
       (Workspace: need N, prefer N*NB) */

    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = *n + 1 - ilo;
    } else {
	icols = irows;
    }
    itau = iwrk;
    iwrk = itau + irows;
    i__1 = *lwork + 1 - iwrk;
    template_lapack_geqrf(&irows, &icols, &b_ref(ilo, ilo), ldb, &work[itau], &work[iwrk], &
	    i__1, &ierr);

/*     Apply the orthogonal transformation to matrix A   
       (Workspace: need N, prefer N*NB) */

    i__1 = *lwork + 1 - iwrk;
    /* Local char arrays added by Elias to get rid of compiler warnings. */
    char str_L[] = {'L', 0};
    char str_T[] = {'T', 0};
    template_lapack_ormqr(str_L, str_T, &irows, &icols, &irows, &b_ref(ilo, ilo), ldb, &work[
	    itau], &a_ref(ilo, ilo), lda, &work[iwrk], &i__1, &ierr);

/*     Initialize VL   
       (Workspace: need N, prefer N*NB) */

    if (ilvl) {
	template_lapack_laset("Full", n, n, &c_b26, &c_b27, &vl[vl_offset], ldvl)
		;
	i__1 = irows - 1;
	i__2 = irows - 1;
	template_lapack_lacpy("L", &i__1, &i__2, &b_ref(ilo + 1, ilo), ldb, &vl_ref(ilo + 1,
		 ilo), ldvl);
	i__1 = *lwork + 1 - iwrk;
	template_lapack_orgqr(&irows, &irows, &irows, &vl_ref(ilo, ilo), ldvl, &work[itau], 
		&work[iwrk], &i__1, &ierr);
    }

/*     Initialize VR */

    if (ilvr) {
	template_lapack_laset("Full", n, n, &c_b26, &c_b27, &vr[vr_offset], ldvr)
		;
    }

/*     Reduce to generalized Hessenberg form   
       (Workspace: none needed) */

    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

	template_lapack_gghrd(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &ierr);
    } else {
	template_lapack_gghrd("N", "N", &irows, &c__1, &irows, &a_ref(ilo, ilo), lda, &
		b_ref(ilo, ilo), ldb, &vl[vl_offset], ldvl, &vr[vr_offset], 
		ldvr, &ierr);
    }

/*     Perform QZ algorithm (Compute eigenvalues, and optionally, the   
       Schur forms and Schur vectors)   
       (Workspace: need N) */

    iwrk = itau;
    if (ilv) {
	*(unsigned char *)chtemp = 'S';
    } else {
	*(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwrk;
    template_lapack_hgeqz(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], 
	    ldvl, &vr[vr_offset], ldvr, &work[iwrk], &i__1, &ierr);
    if (ierr != 0) {
	if (ierr > 0 && ierr <= *n) {
	    *info = ierr;
	} else if (ierr > *n && ierr <= *n << 1) {
	    *info = ierr - *n;
	} else {
	    *info = *n + 1;
	}
	goto L110;
    }

/*     Compute Eigenvectors   
       (Workspace: need 6*N) */

    if (ilv) {
	if (ilvl) {
	    if (ilvr) {
		*(unsigned char *)chtemp = 'B';
	    } else {
		*(unsigned char *)chtemp = 'L';
	    }
	} else {
	    *(unsigned char *)chtemp = 'R';
	}
	template_lapack_tgevc(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwrk], &ierr);
	if (ierr != 0) {
	    *info = *n + 2;
	    goto L110;
	}

/*        Undo balancing on VL and VR and normalization   
          (Workspace: none needed) */

	if (ilvl) {
	    template_lapack_ggbak("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vl[vl_offset], ldvl, &ierr);
	    i__1 = *n;
	    for (jc = 1; jc <= i__1; ++jc) {
		if (alphai[jc] < 0.) {
		    goto L50;
		}
		temp = 0.;
		if (alphai[jc] == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
			d__2 = temp, d__3 = (d__1 = vl_ref(jr, jc), absMACRO(d__1))
				;
			temp = maxMACRO(d__2,d__3);
/* L10: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
			d__3 = temp, d__4 = (d__1 = vl_ref(jr, jc), absMACRO(d__1))
				 + (d__2 = vl_ref(jr, jc + 1), absMACRO(d__2));
			temp = maxMACRO(d__3,d__4);
/* L20: */
		    }
		}
		if (temp < smlnum) {
		    goto L50;
		}
		temp = 1. / temp;
		if (alphai[jc] == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vl_ref(jr, jc) = vl_ref(jr, jc) * temp;
/* L30: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vl_ref(jr, jc) = vl_ref(jr, jc) * temp;
			vl_ref(jr, jc + 1) = vl_ref(jr, jc + 1) * temp;
/* L40: */
		    }
		}
L50:
		;
	    }
	}
	if (ilvr) {
	    template_lapack_ggbak("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vr[vr_offset], ldvr, &ierr);
	    i__1 = *n;
	    for (jc = 1; jc <= i__1; ++jc) {
		if (alphai[jc] < 0.) {
		    goto L100;
		}
		temp = 0.;
		if (alphai[jc] == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
			d__2 = temp, d__3 = (d__1 = vr_ref(jr, jc), absMACRO(d__1))
				;
			temp = maxMACRO(d__2,d__3);
/* L60: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
			d__3 = temp, d__4 = (d__1 = vr_ref(jr, jc), absMACRO(d__1))
				 + (d__2 = vr_ref(jr, jc + 1), absMACRO(d__2));
			temp = maxMACRO(d__3,d__4);
/* L70: */
		    }
		}
		if (temp < smlnum) {
		    goto L100;
		}
		temp = 1. / temp;
		if (alphai[jc] == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vr_ref(jr, jc) = vr_ref(jr, jc) * temp;
/* L80: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vr_ref(jr, jc) = vr_ref(jr, jc) * temp;
			vr_ref(jr, jc + 1) = vr_ref(jr, jc + 1) * temp;
/* L90: */
		    }
		}
L100:
		;
	    }
	}

/*        End of eigenvector calculation */

    }

/*     Undo scaling if necessary */

    if (ilascl) {
	template_lapack_lascl("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &
		ierr);
	template_lapack_lascl("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &
		ierr);
    }

    if (ilbscl) {
	template_lapack_lascl("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &
		ierr);
    }

L110:

    work[1] = (Treal) maxwrk;

    return 0;

/*     End of DGGEV */

} /* dggev_ */

#undef vr_ref
#undef vl_ref
#undef b_ref
#undef a_ref


#endif
