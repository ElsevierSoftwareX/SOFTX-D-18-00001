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
 

#ifndef TEMPLATE_LAPACK_LASQ6_HEADER
#define TEMPLATE_LAPACK_LASQ6_HEADER

template<class Treal>
int template_lapack_lasq6(integer *i0, integer *n0, Treal *z__, 
	integer *pp, Treal *dmin__, Treal *dmin1, Treal *dmin2, 
	 Treal *dn, Treal *dnm1, Treal *dnm2)
{
    /* System generated locals */
    integer i__1;
    Treal d__1, d__2;

    /* Local variables */
    Treal d__;
    integer j4, j4p2;
    Treal emin, temp;
    Treal safmin;


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

/*  DLASQ6 computes one dqd (shift equal to zero) transform in */
/*  ping-pong form, with protection against underflow and overflow. */

/*  Arguments */
/*  ========= */

/*  I0    (input) INTEGER */
/*        First index. */

/*  N0    (input) INTEGER */
/*        Last index. */

/*  Z     (input) DOUBLE PRECISION array, dimension ( 4*N ) */
/*        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid */
/*        an extra argument. */

/*  PP    (input) INTEGER */
/*        PP=0 for ping, PP=1 for pong. */

/*  DMIN  (output) DOUBLE PRECISION */
/*        Minimum value of d. */

/*  DMIN1 (output) DOUBLE PRECISION */
/*        Minimum value of d, excluding D( N0 ). */

/*  DMIN2 (output) DOUBLE PRECISION */
/*        Minimum value of d, excluding D( N0 ) and D( N0-1 ). */

/*  DN    (output) DOUBLE PRECISION */
/*        d(N0), the last value of d. */

/*  DNM1  (output) DOUBLE PRECISION */
/*        d(N0-1). */

/*  DNM2  (output) DOUBLE PRECISION */
/*        d(N0-2). */

/*  ===================================================================== */

/*     .. Parameter .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Function .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --z__;

    /* Function Body */
    if (*n0 - *i0 - 1 <= 0) {
	return 0;
    }

    safmin = template_lapack_lamch("Safe minimum", (Treal)0);
    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4];
    *dmin__ = d__;

    if (*pp == 0) {
      i__1 = (*n0 - 3) << 2;
	for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
	    z__[j4 - 2] = d__ + z__[j4 - 1];
	    if (z__[j4 - 2] == 0.) {
		z__[j4] = 0.;
		d__ = z__[j4 + 1];
		*dmin__ = d__;
		emin = 0.;
	    } else if (safmin * z__[j4 + 1] < z__[j4 - 2] && safmin * z__[j4 
		    - 2] < z__[j4 + 1]) {
		temp = z__[j4 + 1] / z__[j4 - 2];
		z__[j4] = z__[j4 - 1] * temp;
		d__ *= temp;
	    } else {
		z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
		d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]);
	    }
	    *dmin__ = minMACRO(*dmin__,d__);
/* Computing MIN */
	    d__1 = emin, d__2 = z__[j4];
	    emin = minMACRO(d__1,d__2);
/* L10: */
	}
    } else {
      i__1 = ( *n0 - 3 ) << 2;
	for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
	    z__[j4 - 3] = d__ + z__[j4];
	    if (z__[j4 - 3] == 0.) {
		z__[j4 - 1] = 0.;
		d__ = z__[j4 + 2];
		*dmin__ = d__;
		emin = 0.;
	    } else if (safmin * z__[j4 + 2] < z__[j4 - 3] && safmin * z__[j4 
		    - 3] < z__[j4 + 2]) {
		temp = z__[j4 + 2] / z__[j4 - 3];
		z__[j4 - 1] = z__[j4] * temp;
		d__ *= temp;
	    } else {
		z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
		d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]);
	    }
	    *dmin__ = minMACRO(*dmin__,d__);
/* Computing MIN */
	    d__1 = emin, d__2 = z__[j4 - 1];
	    emin = minMACRO(d__1,d__2);
/* L20: */
	}
    }

/*     Unroll last two steps. */

    *dnm2 = d__;
    *dmin2 = *dmin__;
    j4 = ( ( *n0 - 2 ) << 2) - *pp;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm2 + z__[j4p2];
    if (z__[j4 - 2] == 0.) {
	z__[j4] = 0.;
	*dnm1 = z__[j4p2 + 2];
	*dmin__ = *dnm1;
	emin = 0.;
    } else if (safmin * z__[j4p2 + 2] < z__[j4 - 2] && safmin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
	temp = z__[j4p2 + 2] / z__[j4 - 2];
	z__[j4] = z__[j4p2] * temp;
	*dnm1 = *dnm2 * temp;
    } else {
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]);
    }
    *dmin__ = minMACRO(*dmin__,*dnm1);

    *dmin1 = *dmin__;
    j4 += 4;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm1 + z__[j4p2];
    if (z__[j4 - 2] == 0.) {
	z__[j4] = 0.;
	*dn = z__[j4p2 + 2];
	*dmin__ = *dn;
	emin = 0.;
    } else if (safmin * z__[j4p2 + 2] < z__[j4 - 2] && safmin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
	temp = z__[j4p2 + 2] / z__[j4 - 2];
	z__[j4] = z__[j4p2] * temp;
	*dn = *dnm1 * temp;
    } else {
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]);
    }
    *dmin__ = minMACRO(*dmin__,*dn);

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return 0;

/*     End of DLASQ6 */

} /* dlasq6_ */

#endif
