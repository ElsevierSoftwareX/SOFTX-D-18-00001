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
 

#ifndef TEMPLATE_LAPACK_LARRA_HEADER
#define TEMPLATE_LAPACK_LARRA_HEADER

template<class Treal>
int template_lapack_larra(const integer *n, Treal *d__, Treal *e, 
	Treal *e2, Treal *spltol, Treal *tnrm, integer *nsplit, 
	 integer *isplit, integer *info)
{
    /* System generated locals */
    integer i__1;
    Treal d__1, d__2;


    /* Local variables */
    integer i__;
    Treal tmp1, eabs;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  Compute the splitting points with threshold SPLTOL. */
/*  DLARRA sets any "small" off-diagonal elements to zero. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix. N > 0. */

/*  D       (input) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the N diagonal elements of the tridiagonal */
/*          matrix T. */

/*  E       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the first (N-1) entries contain the subdiagonal */
/*          elements of the tridiagonal matrix T; E(N) need not be set. */
/*          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT, */
/*          are set to zero, the other entries of E are untouched. */

/*  E2      (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the first (N-1) entries contain the SQUARES of the */
/*          subdiagonal elements of the tridiagonal matrix T; */
/*          E2(N) need not be set. */
/*          On exit, the entries E2( ISPLIT( I ) ), */
/*          1 <= I <= NSPLIT, have been set to zero */

/*  SPLTOL (input) DOUBLE PRECISION */
/*          The threshold for splitting. Two criteria can be used: */
/*          SPLTOL<0 : criterion based on absolute off-diagonal value */
/*          SPLTOL>0 : criterion that preserves relative accuracy */

/*  TNRM (input) DOUBLE PRECISION */
/*          The norm of the matrix. */

/*  NSPLIT  (output) INTEGER */
/*          The number of blocks T splits into. 1 <= NSPLIT <= N. */

/*  ISPLIT  (output) INTEGER array, dimension (N) */
/*          The splitting points, at which T breaks up into blocks. */
/*          The first block consists of rows/columns 1 to ISPLIT(1), */
/*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2), */
/*          etc., and the NSPLIT-th consists of rows/columns */
/*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. */


/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */

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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --isplit;
    --e2;
    --e;
    --d__;

    /* Function Body */
    *info = 0;
/*     Compute splitting points */
    *nsplit = 1;
    if (*spltol < 0.) {
/*        Criterion based on absolute off-diagonal value */
	tmp1 = absMACRO(*spltol) * *tnrm;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    eabs = (d__1 = e[i__], absMACRO(d__1));
	    if (eabs <= tmp1) {
		e[i__] = 0.;
		e2[i__] = 0.;
		isplit[*nsplit] = i__;
		++(*nsplit);
	    }
/* L9: */
	}
    } else {
/*        Criterion that guarantees relative accuracy */
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    eabs = (d__1 = e[i__], absMACRO(d__1));
	    if (eabs <= *spltol * template_blas_sqrt((d__1 = d__[i__], absMACRO(d__1))) * template_blas_sqrt((
		    d__2 = d__[i__ + 1], absMACRO(d__2)))) {
		e[i__] = 0.;
		e2[i__] = 0.;
		isplit[*nsplit] = i__;
		++(*nsplit);
	    }
/* L10: */
	}
    }
    isplit[*nsplit] = *n;
    return 0;

/*     End of DLARRA */

} /* dlarra_ */

#endif
