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
 

#ifndef TEMPLATE_LAPACK_GEQRF_HEADER
#define TEMPLATE_LAPACK_GEQRF_HEADER


template<class Treal>
int template_lapack_geqrf(const integer *m, const integer *n, Treal *a, const integer *
	lda, Treal *tau, Treal *work, const integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DGEQRF computes a QR factorization of a real M-by-N matrix A:   
    A = Q * R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the elements on and above the diagonal of the array   
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is   
            upper triangular if m >= n); the elements below the diagonal,   
            with the array TAU, represent the orthogonal matrix Q as a   
            product of min(m,n) elementary reflectors (see Further   
            Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors (see Further   
            Details).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is   
            the optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),   
    and tau in TAU(i).   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c_n1 = -1;
     integer c__3 = 3;
     integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
     integer i__, k, nbmin, iinfo;
     integer ib, nb;
     integer nx;
     integer ldwork, lwkopt;
     logical lquery;
     integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = template_lapack_ilaenv(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
    lwkopt = *n * nb;
    work[1] = (Treal) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < maxMACRO(1,*m)) {
	*info = -4;
    } else if (*lwork < maxMACRO(1,*n) && ! lquery) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("GEQRF ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    k = minMACRO(*m,*n);
    if (k == 0) {
	work[1] = 1.;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code.   

   Computing MAX */
	i__1 = 0, i__2 = template_lapack_ilaenv(&c__3, "DGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = maxMACRO(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = template_lapack_ilaenv(&c__2, "DGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
		nbmin = maxMACRO(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = k - i__ + 1;
	    ib = minMACRO(i__3,nb);

/*           Compute the QR factorization of the current block   
             A(i:m,i:i+ib-1) */

	    i__3 = *m - i__ + 1;
	    template_lapack_geqr2(&i__3, &ib, &a_ref(i__, i__), lda, &tau[i__], &work[1], &
		    iinfo);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *m - i__ + 1;
		template_lapack_larft("Forward", "Columnwise", &i__3, &ib, &a_ref(i__, i__),
			 lda, &tau[i__], &work[1], &ldwork);

/*              Apply H' to A(i:m,i+ib:n) from the left */

		i__3 = *m - i__ + 1;
		i__4 = *n - i__ - ib + 1;
		template_lapack_larfb("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &a_ref(i__, i__), lda, &work[1], &ldwork, &
			a_ref(i__, i__ + ib), lda, &work[ib + 1], &ldwork);
	    }
/* L10: */
	}
    } else {
	i__ = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	template_lapack_geqr2(&i__2, &i__1, &a_ref(i__, i__), lda, &tau[i__], &work[1], &
		iinfo);
    }

    work[1] = (Treal) iws;
    return 0;

/*     End of DGEQRF */

} /* dgeqrf_ */

#undef a_ref


#endif
