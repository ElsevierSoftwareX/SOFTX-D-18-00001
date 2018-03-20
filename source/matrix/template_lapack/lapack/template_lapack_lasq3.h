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
 

#ifndef TEMPLATE_LAPACK_LASQ3_HEADER
#define TEMPLATE_LAPACK_LASQ3_HEADER


#include "template_lapack_lasq4.h"
#include "template_lapack_lasq5.h"
#include "template_lapack_lasq6.h"


template<class Treal>
int template_lapack_lasq3(integer *i0, integer *n0, Treal *z__, 
	integer *pp, Treal *dmin__, Treal *sigma, Treal *desig, 
	 Treal *qmax, integer *nfail, integer *iter, integer *ndiv, 
	logical *ieee, integer *ttype, Treal *dmin1, Treal *dmin2, 
	Treal *dn, Treal *dn1, Treal *dn2, Treal *g, 
	Treal *tau)
{
    /* System generated locals */
    integer i__1;
    Treal d__1, d__2;


    /* Local variables */
    Treal s, t;
    integer j4, nn;
    Treal eps, tol;
    integer n0in, ipn4;
    Treal tol2, temp;


/*  -- LAPACK routine (version 3.2)                                    -- */

/*  -- Contributed by Osni Marques of the Lawrence Berkeley National   -- */
/*  -- Laboratory and Beresford Parlett of the Univ. of California at  -- */
/*  -- Berkeley                                                        -- */
/*  -- November 2008                                                   -- */

/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds. */
/*  In case of failure it changes shifts, and tries again until output */
/*  is positive. */

/*  Arguments */
/*  ========= */

/*  I0     (input) INTEGER */
/*         First index. */

/*  N0     (input) INTEGER */
/*         Last index. */

/*  Z      (input) DOUBLE PRECISION array, dimension ( 4*N ) */
/*         Z holds the qd array. */

/*  PP     (input/output) INTEGER */
/*         PP=0 for ping, PP=1 for pong. */
/*         PP=2 indicates that flipping was applied to the Z array */
/*         and that the initial tests for deflation should not be */
/*         performed. */

/*  DMIN   (output) DOUBLE PRECISION */
/*         Minimum value of d. */

/*  SIGMA  (output) DOUBLE PRECISION */
/*         Sum of shifts used in current segment. */

/*  DESIG  (input/output) DOUBLE PRECISION */
/*         Lower order part of SIGMA */

/*  QMAX   (input) DOUBLE PRECISION */
/*         Maximum value of q. */

/*  NFAIL  (output) INTEGER */
/*         Number of times shift was too big. */

/*  ITER   (output) INTEGER */
/*         Number of iterations. */

/*  NDIV   (output) INTEGER */
/*         Number of divisions. */

/*  IEEE   (input) LOGICAL */
/*         Flag for IEEE or non IEEE arithmetic (passed to DLASQ5). */

/*  TTYPE  (input/output) INTEGER */
/*         Shift type. */

/*  DMIN1, DMIN2, DN, DN1, DN2, G, TAU (input/output) DOUBLE PRECISION */
/*         These are passed as arguments in order to save their values */
/*         between calls to DLASQ3. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Function .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --z__;

    /* Function Body */
    n0in = *n0;
    eps = template_lapack_lamch("Precision", (Treal)0);
    tol = eps * 100.;
/* Computing 2nd power */
    d__1 = tol;
    tol2 = d__1 * d__1;

/*     Check for deflation. */

L10:

    if (*n0 < *i0) {
	return 0;
    }
    if (*n0 == *i0) {
	goto L20;
    }
    nn = (*n0 << 2) + *pp;
    if (*n0 == *i0 + 1) {
	goto L40;
    }

/*     Check whether E(N0-1) is negligible, 1 eigenvalue. */

    if (z__[nn - 5] > tol2 * (*sigma + z__[nn - 3]) && z__[nn - (*pp << 1) - 
	    4] > tol2 * z__[nn - 7]) {
	goto L30;
    }

L20:

    z__[(*n0 << 2) - 3] = z__[(*n0 << 2) + *pp - 3] + *sigma;
    --(*n0);
    goto L10;

/*     Check  whether E(N0-2) is negligible, 2 eigenvalues. */

L30:

    if (z__[nn - 9] > tol2 * *sigma && z__[nn - (*pp << 1) - 8] > tol2 * z__[
	    nn - 11]) {
	goto L50;
    }

L40:

    if (z__[nn - 3] > z__[nn - 7]) {
	s = z__[nn - 3];
	z__[nn - 3] = z__[nn - 7];
	z__[nn - 7] = s;
    }
    if (z__[nn - 5] > z__[nn - 3] * tol2) {
	t = (z__[nn - 7] - z__[nn - 3] + z__[nn - 5]) * .5;
	s = z__[nn - 3] * (z__[nn - 5] / t);
	if (s <= t) {
	    s = z__[nn - 3] * (z__[nn - 5] / (t * (template_blas_sqrt(s / t + 1.) + 1.)));
	} else {
	    s = z__[nn - 3] * (z__[nn - 5] / (t + template_blas_sqrt(t) * template_blas_sqrt(t + s)));
	}
	t = z__[nn - 7] + (s + z__[nn - 5]);
	z__[nn - 3] *= z__[nn - 7] / t;
	z__[nn - 7] = t;
    }
    z__[(*n0 << 2) - 7] = z__[nn - 7] + *sigma;
    z__[(*n0 << 2) - 3] = z__[nn - 3] + *sigma;
    *n0 += -2;
    goto L10;

L50:
    if (*pp == 2) {
	*pp = 0;
    }

/*     Reverse the qd-array, if warranted. */

    if (*dmin__ <= 0. || *n0 < n0in) {
	if (z__[(*i0 << 2) + *pp - 3] * 1.5 < z__[(*n0 << 2) + *pp - 3]) {
	  ipn4 = ( *i0 + *n0 ) << 2;
	  i__1 = ( *i0 + *n0 - 1 ) << 1;
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		temp = z__[j4 - 3];
		z__[j4 - 3] = z__[ipn4 - j4 - 3];
		z__[ipn4 - j4 - 3] = temp;
		temp = z__[j4 - 2];
		z__[j4 - 2] = z__[ipn4 - j4 - 2];
		z__[ipn4 - j4 - 2] = temp;
		temp = z__[j4 - 1];
		z__[j4 - 1] = z__[ipn4 - j4 - 5];
		z__[ipn4 - j4 - 5] = temp;
		temp = z__[j4];
		z__[j4] = z__[ipn4 - j4 - 4];
		z__[ipn4 - j4 - 4] = temp;
/* L60: */
	    }
	    if (*n0 - *i0 <= 4) {
		z__[(*n0 << 2) + *pp - 1] = z__[(*i0 << 2) + *pp - 1];
		z__[(*n0 << 2) - *pp] = z__[(*i0 << 2) - *pp];
	    }
/* Computing MIN */
	    d__1 = *dmin2, d__2 = z__[(*n0 << 2) + *pp - 1];
	    *dmin2 = minMACRO(d__1,d__2);
/* Computing MIN */
	    d__1 = z__[(*n0 << 2) + *pp - 1], d__2 = z__[(*i0 << 2) + *pp - 1]
		    , d__1 = minMACRO(d__1,d__2), d__2 = z__[(*i0 << 2) + *pp + 3];
	    z__[(*n0 << 2) + *pp - 1] = minMACRO(d__1,d__2);
/* Computing MIN */
	    d__1 = z__[(*n0 << 2) - *pp], d__2 = z__[(*i0 << 2) - *pp], d__1 =
		     minMACRO(d__1,d__2), d__2 = z__[(*i0 << 2) - *pp + 4];
	    z__[(*n0 << 2) - *pp] = minMACRO(d__1,d__2);
/* Computing MAX */
	    d__1 = *qmax, d__2 = z__[(*i0 << 2) + *pp - 3], d__1 = maxMACRO(d__1,
		    d__2), d__2 = z__[(*i0 << 2) + *pp + 1];
	    *qmax = maxMACRO(d__1,d__2);
	    *dmin__ = -0.;
	}
    }

/*     Choose a shift. */

    template_lapack_lasq4(i0, n0, &z__[1], pp, &n0in, dmin__, dmin1, dmin2, dn, dn1, dn2, 
	    tau, ttype, g);

/*     Call dqds until DMIN > 0. */

L70:

    template_lapack_lasq5(i0, n0, &z__[1], pp, tau, dmin__, dmin1, dmin2, dn, dn1, dn2, 
	    ieee);

    *ndiv += *n0 - *i0 + 2;
    ++(*iter);

/*     Check status. */

    if (*dmin__ >= 0. && *dmin1 > 0.) {

/*        Success. */

	goto L90;

    } else if (*dmin__ < 0. && *dmin1 > 0. && z__[( ( *n0 - 1 ) << 2) - *pp] < tol 
	    * (*sigma + *dn1) && absMACRO(*dn) < tol * *sigma) {

/*        Convergence hidden by negative DN. */

      z__[( ( *n0 - 1 ) << 2) - *pp + 2] = 0.;
	*dmin__ = 0.;
	goto L90;
    } else if (*dmin__ < 0.) {

/*        TAU too big. Select new TAU and try again. */

	++(*nfail);
	if (*ttype < -22) {

/*           Failed twice. Play it safe. */

	    *tau = 0.;
	} else if (*dmin1 > 0.) {

/*           Late failure. Gives excellent shift. */

	    *tau = (*tau + *dmin__) * (1. - eps * 2.);
	    *ttype += -11;
	} else {

/*           Early failure. Divide by 4. */

	    *tau *= .25;
	    *ttype += -12;
	}
	goto L70;
    } else if (template_lapack_isnan(dmin__)) {

/*        NaN. */

	if (*tau == 0.) {
	    goto L80;
	} else {
	    *tau = 0.;
	    goto L70;
	}
    } else {

/*        Possible underflow. Play it safe. */

	goto L80;
    }

/*     Risk of underflow. */

L80:
    template_lapack_lasq6(i0, n0, &z__[1], pp, dmin__, dmin1, dmin2, dn, dn1, dn2);
    *ndiv += *n0 - *i0 + 2;
    ++(*iter);
    *tau = 0.;

L90:
    if (*tau < *sigma) {
	*desig += *tau;
	t = *sigma + *desig;
	*desig -= t - *sigma;
    } else {
	t = *sigma + *tau;
	*desig = *sigma - (t - *tau) + *desig;
    }
    *sigma = t;

    return 0;

/*     End of DLASQ3 */

} /* dlasq3_ */

#endif
