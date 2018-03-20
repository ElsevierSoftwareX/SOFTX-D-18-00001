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
 

#ifndef TEMPLATE_LAPACK_LAR1V_HEADER
#define TEMPLATE_LAPACK_LAR1V_HEADER

template<class Treal>
int template_lapack_lar1v(integer *n, integer *b1, integer *bn, Treal 
	*lambda, Treal *d__, Treal *l, Treal *ld, Treal *
	lld, Treal *pivmin, Treal *gaptol, Treal *z__, logical 
	*wantnc, integer *negcnt, Treal *ztz, Treal *mingma, 
	integer *r__, integer *isuppz, Treal *nrminv, Treal *resid, 
	Treal *rqcorr, Treal *work)
{
    /* System generated locals */
    integer i__1;
    Treal d__1, d__2, d__3;


    /* Local variables */
    integer i__;
    Treal s;
    integer r1, r2;
    Treal eps, tmp;
    integer neg1, neg2, indp, inds;
    Treal dplus;
    integer indlpl, indumn;
    Treal dminus;
    logical sawnan1, sawnan2;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLAR1V computes the (scaled) r-th column of the inverse of */
/*  the sumbmatrix in rows B1 through BN of the tridiagonal matrix */
/*  L D L^T - sigma I. When sigma is close to an eigenvalue, the */
/*  computed vector is an accurate eigenvector. Usually, r corresponds */
/*  to the index where the eigenvector is largest in magnitude. */
/*  The following steps accomplish this computation : */
/*  (a) Stationary qd transform,  L D L^T - sigma I = L(+) D(+) L(+)^T, */
/*  (b) Progressive qd transform, L D L^T - sigma I = U(-) D(-) U(-)^T, */
/*  (c) Computation of the diagonal elements of the inverse of */
/*      L D L^T - sigma I by combining the above transforms, and choosing */
/*      r as the index where the diagonal of the inverse is (one of the) */
/*      largest in magnitude. */
/*  (d) Computation of the (scaled) r-th column of the inverse using the */
/*      twisted factorization obtained by combining the top part of the */
/*      the stationary and the bottom part of the progressive transform. */

/*  Arguments */
/*  ========= */

/*  N        (input) INTEGER */
/*           The order of the matrix L D L^T. */

/*  B1       (input) INTEGER */
/*           First index of the submatrix of L D L^T. */

/*  BN       (input) INTEGER */
/*           Last index of the submatrix of L D L^T. */

/*  LAMBDA    (input) DOUBLE PRECISION */
/*           The shift. In order to compute an accurate eigenvector, */
/*           LAMBDA should be a good approximation to an eigenvalue */
/*           of L D L^T. */

/*  L        (input) DOUBLE PRECISION array, dimension (N-1) */
/*           The (n-1) subdiagonal elements of the unit bidiagonal matrix */
/*           L, in elements 1 to N-1. */

/*  D        (input) DOUBLE PRECISION array, dimension (N) */
/*           The n diagonal elements of the diagonal matrix D. */

/*  LD       (input) DOUBLE PRECISION array, dimension (N-1) */
/*           The n-1 elements L(i)*D(i). */

/*  LLD      (input) DOUBLE PRECISION array, dimension (N-1) */
/*           The n-1 elements L(i)*L(i)*D(i). */

/*  PIVMIN   (input) DOUBLE PRECISION */
/*           The minimum pivot in the Sturm sequence. */

/*  GAPTOL   (input) DOUBLE PRECISION */
/*           Tolerance that indicates when eigenvector entries are negligible */
/*           w.r.t. their contribution to the residual. */

/*  Z        (input/output) DOUBLE PRECISION array, dimension (N) */
/*           On input, all entries of Z must be set to 0. */
/*           On output, Z contains the (scaled) r-th column of the */
/*           inverse. The scaling is such that Z(R) equals 1. */

/*  WANTNC   (input) LOGICAL */
/*           Specifies whether NEGCNT has to be computed. */

/*  NEGCNT   (output) INTEGER */
/*           If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin */
/*           in the  matrix factorization L D L^T, and NEGCNT = -1 otherwise. */

/*  ZTZ      (output) DOUBLE PRECISION */
/*           The square of the 2-norm of Z. */

/*  MINGMA   (output) DOUBLE PRECISION */
/*           The reciprocal of the largest (in magnitude) diagonal */
/*           element of the inverse of L D L^T - sigma I. */

/*  R        (input/output) INTEGER */
/*           The twist index for the twisted factorization used to */
/*           compute Z. */
/*           On input, 0 <= R <= N. If R is input as 0, R is set to */
/*           the index where (L D L^T - sigma I)^{-1} is largest */
/*           in magnitude. If 1 <= R <= N, R is unchanged. */
/*           On output, R contains the twist index used to compute Z. */
/*           Ideally, R designates the position of the maximum entry in the */
/*           eigenvector. */

/*  ISUPPZ   (output) INTEGER array, dimension (2) */
/*           The support of the vector in Z, i.e., the vector Z is */
/*           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ). */

/*  NRMINV   (output) DOUBLE PRECISION */
/*           NRMINV = 1/SQRT( ZTZ ) */

/*  RESID    (output) DOUBLE PRECISION */
/*           The residual of the FP vector. */
/*           RESID = ABS( MINGMA )/SQRT( ZTZ ) */

/*  RQCORR   (output) DOUBLE PRECISION */
/*           The Rayleigh Quotient correction to LAMBDA. */
/*           RQCORR = MINGMA*TMP */

/*  WORK     (workspace) DOUBLE PRECISION array, dimension (4*N) */

/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Beresford Parlett, University of California, Berkeley, USA */
/*     Jim Demmel, University of California, Berkeley, USA */
/*     Inderjit Dhillon, University of Texas, Austin, USA */
/*     Osni Marques, LBNL/NERSC, USA */
/*     Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --work;
    --isuppz;
    --z__;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    eps = template_lapack_lamch("Precision", (Treal)0);
    if (*r__ == 0) {
	r1 = *b1;
	r2 = *bn;
    } else {
	r1 = *r__;
	r2 = *r__;
    }
/*     Storage for LPLUS */
    indlpl = 0;
/*     Storage for UMINUS */
    indumn = *n;
    inds = (*n << 1) + 1;
    indp = *n * 3 + 1;
    if (*b1 == 1) {
	work[inds] = 0.;
    } else {
	work[inds + *b1 - 1] = lld[*b1 - 1];
    }

/*     Compute the stationary transform (using the differential form) */
/*     until the index R2. */

    sawnan1 = FALSE_;
    neg1 = 0;
    s = work[inds + *b1 - 1] - *lambda;
    i__1 = r1 - 1;
    for (i__ = *b1; i__ <= i__1; ++i__) {
	dplus = d__[i__] + s;
	work[indlpl + i__] = ld[i__] / dplus;
	if (dplus < 0.) {
	    ++neg1;
	}
	work[inds + i__] = s * work[indlpl + i__] * l[i__];
	s = work[inds + i__] - *lambda;
/* L50: */
    }
    sawnan1 = template_lapack_isnan(&s);
    if (sawnan1) {
	goto L60;
    }
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
	dplus = d__[i__] + s;
	work[indlpl + i__] = ld[i__] / dplus;
	work[inds + i__] = s * work[indlpl + i__] * l[i__];
	s = work[inds + i__] - *lambda;
/* L51: */
    }
    sawnan1 = template_lapack_isnan(&s);

L60:
    if (sawnan1) {
/*        Runs a slower version of the above loop if a NaN is detected */
	neg1 = 0;
	s = work[inds + *b1 - 1] - *lambda;
	i__1 = r1 - 1;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    dplus = d__[i__] + s;
	    if (absMACRO(dplus) < *pivmin) {
		dplus = -(*pivmin);
	    }
	    work[indlpl + i__] = ld[i__] / dplus;
	    if (dplus < 0.) {
		++neg1;
	    }
	    work[inds + i__] = s * work[indlpl + i__] * l[i__];
	    if (work[indlpl + i__] == 0.) {
		work[inds + i__] = lld[i__];
	    }
	    s = work[inds + i__] - *lambda;
/* L70: */
	}
	i__1 = r2 - 1;
	for (i__ = r1; i__ <= i__1; ++i__) {
	    dplus = d__[i__] + s;
	    if (absMACRO(dplus) < *pivmin) {
		dplus = -(*pivmin);
	    }
	    work[indlpl + i__] = ld[i__] / dplus;
	    work[inds + i__] = s * work[indlpl + i__] * l[i__];
	    if (work[indlpl + i__] == 0.) {
		work[inds + i__] = lld[i__];
	    }
	    s = work[inds + i__] - *lambda;
/* L71: */
	}
    }

/*     Compute the progressive transform (using the differential form) */
/*     until the index R1 */

    sawnan2 = FALSE_;
    neg2 = 0;
    work[indp + *bn - 1] = d__[*bn] - *lambda;
    i__1 = r1;
    for (i__ = *bn - 1; i__ >= i__1; --i__) {
	dminus = lld[i__] + work[indp + i__];
	tmp = d__[i__] / dminus;
	if (dminus < 0.) {
	    ++neg2;
	}
	work[indumn + i__] = l[i__] * tmp;
	work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
/* L80: */
    }
    tmp = work[indp + r1 - 1];
    sawnan2 = template_lapack_isnan(&tmp);
    if (sawnan2) {
/*        Runs a slower version of the above loop if a NaN is detected */
	neg2 = 0;
	i__1 = r1;
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
	    dminus = lld[i__] + work[indp + i__];
	    if (absMACRO(dminus) < *pivmin) {
		dminus = -(*pivmin);
	    }
	    tmp = d__[i__] / dminus;
	    if (dminus < 0.) {
		++neg2;
	    }
	    work[indumn + i__] = l[i__] * tmp;
	    work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
	    if (tmp == 0.) {
		work[indp + i__ - 1] = d__[i__] - *lambda;
	    }
/* L100: */
	}
    }

/*     Find the index (from R1 to R2) of the largest (in magnitude) */
/*     diagonal element of the inverse */

    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
    if (*mingma < 0.) {
	++neg1;
    }
    if (*wantnc) {
	*negcnt = neg1 + neg2;
    } else {
	*negcnt = -1;
    }
    if (absMACRO(*mingma) == 0.) {
	*mingma = eps * work[inds + r1 - 1];
    }
    *r__ = r1;
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
	tmp = work[inds + i__] + work[indp + i__];
	if (tmp == 0.) {
	    tmp = eps * work[inds + i__];
	}
	if (absMACRO(tmp) <= absMACRO(*mingma)) {
	    *mingma = tmp;
	    *r__ = i__ + 1;
	}
/* L110: */
    }

/*     Compute the FP vector: solve N^T v = e_r */

    isuppz[1] = *b1;
    isuppz[2] = *bn;
    z__[*r__] = 1.;
    *ztz = 1.;

/*     Compute the FP vector upwards from R */

    if (! sawnan1 && ! sawnan2) {
	i__1 = *b1;
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	    z__[i__] = -(work[indlpl + i__] * z__[i__ + 1]);
	    if (((d__1 = z__[i__], absMACRO(d__1)) + (d__2 = z__[i__ + 1], absMACRO(
		    d__2))) * (d__3 = ld[i__], absMACRO(d__3)) < *gaptol) {
		z__[i__] = 0.;
		isuppz[1] = i__ + 1;
		goto L220;
	    }
	    *ztz += z__[i__] * z__[i__];
/* L210: */
	}
L220:
	;
    } else {
/*        Run slower loop if NaN occurred. */
	i__1 = *b1;
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	    if (z__[i__ + 1] == 0.) {
		z__[i__] = -(ld[i__ + 1] / ld[i__]) * z__[i__ + 2];
	    } else {
		z__[i__] = -(work[indlpl + i__] * z__[i__ + 1]);
	    }
	    if (((d__1 = z__[i__], absMACRO(d__1)) + (d__2 = z__[i__ + 1], absMACRO(
		    d__2))) * (d__3 = ld[i__], absMACRO(d__3)) < *gaptol) {
		z__[i__] = 0.;
		isuppz[1] = i__ + 1;
		goto L240;
	    }
	    *ztz += z__[i__] * z__[i__];
/* L230: */
	}
L240:
	;
    }
/*     Compute the FP vector downwards from R in blocks of size BLKSIZ */
    if (! sawnan1 && ! sawnan2) {
	i__1 = *bn - 1;
	for (i__ = *r__; i__ <= i__1; ++i__) {
	    z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
	    if (((d__1 = z__[i__], absMACRO(d__1)) + (d__2 = z__[i__ + 1], absMACRO(
		    d__2))) * (d__3 = ld[i__], absMACRO(d__3)) < *gaptol) {
		z__[i__ + 1] = 0.;
		isuppz[2] = i__;
		goto L260;
	    }
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
/* L250: */
	}
L260:
	;
    } else {
/*        Run slower loop if NaN occurred. */
	i__1 = *bn - 1;
	for (i__ = *r__; i__ <= i__1; ++i__) {
	    if (z__[i__] == 0.) {
		z__[i__ + 1] = -(ld[i__ - 1] / ld[i__]) * z__[i__ - 1];
	    } else {
		z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
	    }
	    if (((d__1 = z__[i__], absMACRO(d__1)) + (d__2 = z__[i__ + 1], absMACRO(
		    d__2))) * (d__3 = ld[i__], absMACRO(d__3)) < *gaptol) {
		z__[i__ + 1] = 0.;
		isuppz[2] = i__;
		goto L280;
	    }
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
/* L270: */
	}
L280:
	;
    }

/*     Compute quantities for convergence test */

    tmp = 1. / *ztz;
    *nrminv = template_blas_sqrt(tmp);
    *resid = absMACRO(*mingma) * *nrminv;
    *rqcorr = *mingma * tmp;


    return 0;

/*     End of DLAR1V */

} /* dlar1v_ */

#endif
