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
 

#ifndef TEMPLATE_LAPACK_LAG2_HEADER
#define TEMPLATE_LAPACK_LAG2_HEADER


template<class Treal>
int template_lapack_lag2(const Treal *a, const integer *lda, const Treal *b, 
	const integer *ldb, const Treal *safmin, Treal *scale1, Treal *
	scale2, Treal *wr1, Treal *wr2, Treal *wi)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue   
    problem  A - w B, with scaling as necessary to avoid over-/underflow.   

    The scaling factor "s" results in a modified eigenvalue equation   

        s A - w B   

    where  s  is a non-negative scaling factor chosen so that  w,  w B,   
    and  s A  do not overflow and, if possible, do not underflow, either.   

    Arguments   
    =========   

    A       (input) DOUBLE PRECISION array, dimension (LDA, 2)   
            On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm   
            is less than 1/SAFMIN.  Entries less than   
            sqrt(SAFMIN)*norm(A) are subject to being treated as zero.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= 2.   

    B       (input) DOUBLE PRECISION array, dimension (LDB, 2)   
            On entry, the 2 x 2 upper triangular matrix B.  It is   
            assumed that the one-norm of B is less than 1/SAFMIN.  The   
            diagonals should be at least sqrt(SAFMIN) times the largest   
            element of B (in absolute value); if a diagonal is smaller   
            than that, then  +/- sqrt(SAFMIN) will be used instead of   
            that diagonal.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= 2.   

    SAFMIN  (input) DOUBLE PRECISION   
            The smallest positive number s.t. 1/SAFMIN does not   
            overflow.  (This should always be DLAMCH('S') -- it is an   
            argument in order to avoid having to call DLAMCH frequently.)   

    SCALE1  (output) DOUBLE PRECISION   
            A scaling factor used to avoid over-/underflow in the   
            eigenvalue equation which defines the first eigenvalue.  If   
            the eigenvalues are complex, then the eigenvalues are   
            ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the   
            exponent range of the machine), SCALE1=SCALE2, and SCALE1   
            will always be positive.  If the eigenvalues are real, then   
            the first (real) eigenvalue is  WR1 / SCALE1 , but this may   
            overflow or underflow, and in fact, SCALE1 may be zero or   
            less than the underflow threshhold if the exact eigenvalue   
            is sufficiently large.   

    SCALE2  (output) DOUBLE PRECISION   
            A scaling factor used to avoid over-/underflow in the   
            eigenvalue equation which defines the second eigenvalue.  If   
            the eigenvalues are complex, then SCALE2=SCALE1.  If the   
            eigenvalues are real, then the second (real) eigenvalue is   
            WR2 / SCALE2 , but this may overflow or underflow, and in   
            fact, SCALE2 may be zero or less than the underflow   
            threshhold if the exact eigenvalue is sufficiently large.   

    WR1     (output) DOUBLE PRECISION   
            If the eigenvalue is real, then WR1 is SCALE1 times the   
            eigenvalue closest to the (2,2) element of A B**(-1).  If the   
            eigenvalue is complex, then WR1=WR2 is SCALE1 times the real   
            part of the eigenvalues.   

    WR2     (output) DOUBLE PRECISION   
            If the eigenvalue is real, then WR2 is SCALE2 times the   
            other eigenvalue.  If the eigenvalue is complex, then   
            WR1=WR2 is SCALE1 times the real part of the eigenvalues.   

    WI      (output) DOUBLE PRECISION   
            If the eigenvalue is real, then WI is zero.  If the   
            eigenvalue is complex, then WI is SCALE1 times the imaginary   
            part of the eigenvalues.  WI will always be non-negative.   

    =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;
    Treal d__1, d__2, d__3, d__4, d__5, d__6;
    /* Local variables */
     Treal diff, bmin, wbig, wabs, wdet, r__, binv11, binv22, 
	    discr, anorm, bnorm, bsize, shift, c1, c2, c3, c4, c5, rtmin, 
	    rtmax, wsize, s1, s2, a11, a12, a21, a22, b11, b12, b22, ascale, 
	    bscale, pp, qq, ss, wscale, safmax, wsmall, as11, as12, as22, sum,
	     abi22;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    rtmin = template_blas_sqrt(*safmin);
    rtmax = 1. / rtmin;
    safmax = 1. / *safmin;

/*     Scale A   

   Computing MAX */
    d__5 = (d__1 = a_ref(1, 1), absMACRO(d__1)) + (d__2 = a_ref(2, 1), absMACRO(d__2)), 
	    d__6 = (d__3 = a_ref(1, 2), absMACRO(d__3)) + (d__4 = a_ref(2, 2), absMACRO(
	    d__4)), d__5 = maxMACRO(d__5,d__6);
    anorm = maxMACRO(d__5,*safmin);
    ascale = 1. / anorm;
    a11 = ascale * a_ref(1, 1);
    a21 = ascale * a_ref(2, 1);
    a12 = ascale * a_ref(1, 2);
    a22 = ascale * a_ref(2, 2);

/*     Perturb B if necessary to insure non-singularity */

    b11 = b_ref(1, 1);
    b12 = b_ref(1, 2);
    b22 = b_ref(2, 2);
/* Computing MAX */
    d__1 = absMACRO(b11), d__2 = absMACRO(b12), d__1 = maxMACRO(d__1,d__2), d__2 = absMACRO(b22), 
	    d__1 = maxMACRO(d__1,d__2);
    bmin = rtmin * maxMACRO(d__1,rtmin);
    if (absMACRO(b11) < bmin) {
	b11 = template_lapack_d_sign(&bmin, &b11);
    }
    if (absMACRO(b22) < bmin) {
	b22 = template_lapack_d_sign(&bmin, &b22);
    }

/*     Scale B   

   Computing MAX */
    d__1 = absMACRO(b11), d__2 = absMACRO(b12) + absMACRO(b22), d__1 = maxMACRO(d__1,d__2);
    bnorm = maxMACRO(d__1,*safmin);
/* Computing MAX */
    d__1 = absMACRO(b11), d__2 = absMACRO(b22);
    bsize = maxMACRO(d__1,d__2);
    bscale = 1. / bsize;
    b11 *= bscale;
    b12 *= bscale;
    b22 *= bscale;

/*     Compute larger eigenvalue by method described by C. van Loan   

       ( AS is A shifted by -SHIFT*B ) */

    binv11 = 1. / b11;
    binv22 = 1. / b22;
    s1 = a11 * binv11;
    s2 = a22 * binv22;
    if (absMACRO(s1) <= absMACRO(s2)) {
	as12 = a12 - s1 * b12;
	as22 = a22 - s1 * b22;
	ss = a21 * (binv11 * binv22);
	abi22 = as22 * binv22 - ss * b12;
	pp = abi22 * .5;
	shift = s1;
    } else {
	as12 = a12 - s2 * b12;
	as11 = a11 - s2 * b11;
	ss = a21 * (binv11 * binv22);
	abi22 = -ss * b12;
	pp = (as11 * binv11 + abi22) * .5;
	shift = s2;
    }
    qq = ss * as12;
    if ((d__1 = pp * rtmin, absMACRO(d__1)) >= 1.) {
/* Computing 2nd power */
	d__1 = rtmin * pp;
	discr = d__1 * d__1 + qq * *safmin;
	r__ = template_blas_sqrt((absMACRO(discr))) * rtmax;
    } else {
/* Computing 2nd power */
	d__1 = pp;
	if (d__1 * d__1 + absMACRO(qq) <= *safmin) {
/* Computing 2nd power */
	    d__1 = rtmax * pp;
	    discr = d__1 * d__1 + qq * safmax;
	    r__ = template_blas_sqrt((absMACRO(discr))) * rtmin;
	} else {
/* Computing 2nd power */
	    d__1 = pp;
	    discr = d__1 * d__1 + qq;
	    r__ = template_blas_sqrt((absMACRO(discr)));
	}
    }

/*     Note: the test of R in the following IF is to cover the case when   
             DISCR is small and negative and is flushed to zero during   
             the calculation of R.  On machines which have a consistent   
             flush-to-zero threshhold and handle numbers above that   
             threshhold correctly, it would not be necessary. */

    if (discr >= 0. || r__ == 0.) {
	sum = pp + template_lapack_d_sign(&r__, &pp);
	diff = pp - template_lapack_d_sign(&r__, &pp);
	wbig = shift + sum;

/*        Compute smaller eigenvalue */

	wsmall = shift + diff;
/* Computing MAX */
	d__1 = absMACRO(wsmall);
	if (absMACRO(wbig) * .5 > maxMACRO(d__1,*safmin)) {
	    wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22);
	    wsmall = wdet / wbig;
	}

/*        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)   
          for WR1. */

	if (pp > abi22) {
	    *wr1 = minMACRO(wbig,wsmall);
	    *wr2 = maxMACRO(wbig,wsmall);
	} else {
	    *wr1 = maxMACRO(wbig,wsmall);
	    *wr2 = minMACRO(wbig,wsmall);
	}
	*wi = 0.;
    } else {

/*        Complex eigenvalues */

	*wr1 = shift + pp;
	*wr2 = *wr1;
	*wi = r__;
    }

/*     Further scaling to avoid underflow and overflow in computing   
       SCALE1 and overflow in computing w*B.   

       This scale factor (WSCALE) is bounded from above using C1 and C2,   
       and from below using C3 and C4.   
          C1 implements the condition  s A  must never overflow.   
          C2 implements the condition  w B  must never overflow.   
          C3, with C2,   
             implement the condition that s A - w B must never overflow.   
          C4 implements the condition  s    should not underflow.   
          C5 implements the condition  max(s,|w|) should be at least 2. */

    c1 = bsize * (*safmin * maxMACRO(1.,ascale));
    c2 = *safmin * maxMACRO(1.,bnorm);
    c3 = bsize * *safmin;
    if (ascale <= 1. && bsize <= 1.) {
/* Computing MIN */
	d__1 = 1., d__2 = ascale / *safmin * bsize;
	c4 = minMACRO(d__1,d__2);
    } else {
	c4 = 1.;
    }
    if (ascale <= 1. || bsize <= 1.) {
/* Computing MIN */
	d__1 = 1., d__2 = ascale * bsize;
	c5 = minMACRO(d__1,d__2);
    } else {
	c5 = 1.;
    }

/*     Scale first eigenvalue */

    wabs = absMACRO(*wr1) + absMACRO(*wi);
/* Computing MAX   
   Computing MIN */
    d__3 = c4, d__4 = maxMACRO(wabs,c5) * .5;
    d__1 = maxMACRO(*safmin,c1), d__2 = (wabs * c2 + c3) * 1.0000100000000001, 
	    d__1 = maxMACRO(d__1,d__2), d__2 = minMACRO(d__3,d__4);
    wsize = maxMACRO(d__1,d__2);
    if (wsize != 1.) {
	wscale = 1. / wsize;
	if (wsize > 1.) {
	    *scale1 = maxMACRO(ascale,bsize) * wscale * minMACRO(ascale,bsize);
	} else {
	    *scale1 = minMACRO(ascale,bsize) * wscale * maxMACRO(ascale,bsize);
	}
	*wr1 *= wscale;
	if (*wi != 0.) {
	    *wi *= wscale;
	    *wr2 = *wr1;
	    *scale2 = *scale1;
	}
    } else {
	*scale1 = ascale * bsize;
	*scale2 = *scale1;
    }

/*     Scale second eigenvalue (if real) */

    if (*wi == 0.) {
/* Computing MAX   
   Computing MIN   
   Computing MAX */
	d__5 = absMACRO(*wr2);
	d__3 = c4, d__4 = maxMACRO(d__5,c5) * .5;
	d__1 = maxMACRO(*safmin,c1), d__2 = (absMACRO(*wr2) * c2 + c3) * 
		1.0000100000000001, d__1 = maxMACRO(d__1,d__2), d__2 = minMACRO(d__3,
		d__4);
	wsize = maxMACRO(d__1,d__2);
	if (wsize != 1.) {
	    wscale = 1. / wsize;
	    if (wsize > 1.) {
		*scale2 = maxMACRO(ascale,bsize) * wscale * minMACRO(ascale,bsize);
	    } else {
		*scale2 = minMACRO(ascale,bsize) * wscale * maxMACRO(ascale,bsize);
	    }
	    *wr2 *= wscale;
	} else {
	    *scale2 = ascale * bsize;
	}
    }

/*     End of DLAG2 */

    return 0;
} /* dlag2_ */

#undef b_ref
#undef a_ref


#endif
