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
 

#ifndef TEMPLATE_LAPACK_ORGQL_HEADER
#define TEMPLATE_LAPACK_ORGQL_HEADER


template<class Treal>
int template_lapack_orgql(const integer *m, const integer *n, const integer *k, Treal *
	a, const integer *lda, const Treal *tau, Treal *work, const integer *lwork, 
	integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DORGQL generates an M-by-N real matrix Q with orthonormal columns,   
    which is defined as the last N columns of a product of K elementary   
    reflectors of order M   

          Q  =  H(k) . . . H(2) H(1)   

    as returned by DGEQLF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the   
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the (n-k+i)-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as   
            returned by DGEQLF in the last k columns of its array   
            argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQLF.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

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
     integer i__, j, l, nbmin, iinfo;
     integer ib, nb, kk;
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
    nb = template_lapack_ilaenv(&c__1, "DORGQL", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = maxMACRO(1,*n) * nb;
    work[1] = (Treal) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < maxMACRO(1,*m)) {
	*info = -5;
    } else if (*lwork < maxMACRO(1,*n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("ORGQL ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	work[1] = 1.;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.   

   Computing MAX */
	i__1 = 0, i__2 = template_lapack_ilaenv(&c__3, "DORGQL", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = maxMACRO(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = template_lapack_ilaenv(&c__2, "DORGQL", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
		nbmin = maxMACRO(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block.   
          The last kk columns are handled by the block method.   

   Computing MIN */
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
	kk = minMACRO(i__1,i__2);

/*        Set A(m-kk+1:m,1:n-kk) to zero. */

	i__1 = *n - kk;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = *m - kk + 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = 0.;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the first or only block. */

    i__1 = *m - kk;
    i__2 = *n - kk;
    i__3 = *k - kk;
    template_lapack_org2l(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
	    ;

    if (kk > 0) {

/*        Use blocked code */

	i__1 = *k;
	i__2 = nb;
	for (i__ = *k - kk + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {
/* Computing MIN */
	    i__3 = nb, i__4 = *k - i__ + 1;
	    ib = minMACRO(i__3,i__4);
	    if (*n - *k + i__ > 1) {

/*              Form the triangular factor of the block reflector   
                H = H(i+ib-1) . . . H(i+1) H(i) */

		i__3 = *m - *k + i__ + ib - 1;
		template_lapack_larft("Backward", "Columnwise", &i__3, &ib, &a_ref(1, *n - *
			k + i__), lda, &tau[i__], &work[1], &ldwork);

/*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

		i__3 = *m - *k + i__ + ib - 1;
		i__4 = *n - *k + i__ - 1;
		template_lapack_larfb("Left", "No transpose", "Backward", "Columnwise", &
			i__3, &i__4, &ib, &a_ref(1, *n - *k + i__), lda, &
			work[1], &ldwork, &a[a_offset], lda, &work[ib + 1], &
			ldwork);
	    }

/*           Apply H to rows 1:m-k+i+ib-1 of current block */

	    i__3 = *m - *k + i__ + ib - 1;
	    template_lapack_org2l(&i__3, &ib, &ib, &a_ref(1, *n - *k + i__), lda, &tau[i__],
		     &work[1], &iinfo);

/*           Set rows m-k+i+ib:m of current block to zero */

	    i__3 = *n - *k + i__ + ib - 1;
	    for (j = *n - *k + i__; j <= i__3; ++j) {
		i__4 = *m;
		for (l = *m - *k + i__ + ib; l <= i__4; ++l) {
		    a_ref(l, j) = 0.;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    work[1] = (Treal) iws;
    return 0;

/*     End of DORGQL */

} /* dorgql_ */

#undef a_ref


#endif
