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
 

#ifndef TEMPLATE_LAPACK_ORGTR_HEADER
#define TEMPLATE_LAPACK_ORGTR_HEADER


template<class Treal>
int template_lapack_orgtr(const char *uplo, const integer *n, Treal *a, const integer *
	lda, const Treal *tau, Treal *work, const integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DORGTR generates a real orthogonal matrix Q which is defined as the   
    product of n-1 elementary reflectors of order N, as returned by   
    DSYTRD:   

    if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),   

    if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangle of A contains elementary reflectors   
                   from DSYTRD;   
            = 'L': Lower triangle of A contains elementary reflectors   
                   from DSYTRD.   

    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors,   
            as returned by DSYTRD.   
            On exit, the N-by-N orthogonal matrix Q.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,N).   

    TAU     (input) DOUBLE PRECISION array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DSYTRD.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,N-1).   
            For optimum performance LWORK >= (N-1)*NB, where NB is   
            the optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     integer c_n1 = -1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
     integer i__, j;
     integer iinfo;
     logical upper;
     integer nb;
     integer lwkopt;
     logical lquery;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Initialization added by Elias to get rid of compiler warnings. */
    lwkopt = 0;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    upper = template_blas_lsame(uplo, "U");
    if (! upper && ! template_blas_lsame(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	if (*lwork < maxMACRO(i__1,i__2) && ! lquery) {
	    *info = -7;
	}
    }

    if (*info == 0) {
	if (upper) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    nb = template_lapack_ilaenv(&c__1, "DORGQL", " ", &i__1, &i__2, &i__3, &c_n1, (
		    ftnlen)6, (ftnlen)1);
	} else {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    nb = template_lapack_ilaenv(&c__1, "DORGQR", " ", &i__1, &i__2, &i__3, &c_n1, (
		    ftnlen)6, (ftnlen)1);
	}
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	lwkopt = maxMACRO(i__1,i__2) * nb;
	work[1] = (Treal) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("ORGTR ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	work[1] = 1.;
	return 0;
    }

    if (upper) {

/*        Q was determined by a call to DSYTRD with UPLO = 'U'   

          Shift the vectors which define the elementary reflectors one   
          column to the left, and set the last row and column of Q to   
          those of the unit matrix */

	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j + 1);
/* L10: */
	    }
	    a_ref(*n, j) = 0.;
/* L20: */
	}
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    a_ref(i__, *n) = 0.;
/* L30: */
	}
	a_ref(*n, *n) = 1.;

/*        Generate Q(1:n-1,1:n-1) */

	i__1 = *n - 1;
	i__2 = *n - 1;
	i__3 = *n - 1;
	template_lapack_orgql(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], 
		lwork, &iinfo);

    } else {

/*        Q was determined by a call to DSYTRD with UPLO = 'L'.   

          Shift the vectors which define the elementary reflectors one   
          column to the right, and set the first row and column of Q to   
          those of the unit matrix */

	for (j = *n; j >= 2; --j) {
	    a_ref(1, j) = 0.;
	    i__1 = *n;
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
		a_ref(i__, j) = a_ref(i__, j - 1);
/* L40: */
	    }
/* L50: */
	}
	a_ref(1, 1) = 1.;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    a_ref(i__, 1) = 0.;
/* L60: */
	}
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    template_lapack_orgqr(
				  &i__1, 
				  &i__2, 
				  &i__3, 
				  &a_ref(2, 2), 
				  lda, 
				  &tau[1], 
				  &work[1],
				  lwork, 
				  &iinfo
				  );
	}
    }
    work[1] = (Treal) lwkopt;
    return 0;

/*     End of DORGTR */

} /* dorgtr_ */

#undef a_ref


#endif
