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
 

#ifndef TEMPLATE_LAPACK_LARRR_HEADER
#define TEMPLATE_LAPACK_LARRR_HEADER

template<class Treal>
int template_lapack_larrr(const integer *n, Treal *d__, Treal *e, 
	integer *info)
{
    /* System generated locals */
    integer i__1;
    Treal d__1;


    /* Local variables */
    integer i__;
    Treal eps, tmp, tmp2, rmin;
    Treal offdig, safmin;
    logical yesrel;
    Treal smlnum, offdig2;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */


/*  Purpose */
/*  ======= */

/*  Perform tests to decide whether the symmetric tridiagonal matrix T */
/*  warrants expensive computations which guarantee high relative accuracy */
/*  in the eigenvalues. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix. N > 0. */

/*  D       (input) DOUBLE PRECISION array, dimension (N) */
/*          The N diagonal elements of the tridiagonal matrix T. */

/*  E       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the first (N-1) entries contain the subdiagonal */
/*          elements of the tridiagonal matrix T; E(N) is set to ZERO. */

/*  INFO    (output) INTEGER */
/*          INFO = 0(default) : the matrix warrants computations preserving */
/*                              relative accuracy. */
/*          INFO = 1          : the matrix warrants computations guaranteeing */
/*                              only absolute accuracy. */

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

/*     As a default, do NOT go for relative-accuracy preserving computations. */
    /* Parameter adjustments */
    --e;
    --d__;

    /* Function Body */
    *info = 1;
    safmin = template_lapack_lamch("Safe minimum", (Treal)0);
    eps = template_lapack_lamch("Precision", (Treal)0);
    smlnum = safmin / eps;
    rmin = template_blas_sqrt(smlnum);
/*     Tests for relative accuracy */

/*     Test for scaled diagonal dominance */
/*     Scale the diagonal entries to one and check whether the sum of the */
/*     off-diagonals is less than one */

/*     The sdd relative error bounds have a 1/(1- 2*x) factor in them, */
/*     x = max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative */
/*     accuracy is promised.  In the notation of the code fragment below, */
/*     1/(1 - (OFFDIG + OFFDIG2)) is the condition number. */
/*     We don't think it is worth going into "sdd mode" unless the relative */
/*     condition number is reasonable, not 1/macheps. */
/*     The threshold should be compatible with other thresholds used in the */
/*     code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds */
/*     to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000 */
/*     instead of the current OFFDIG + OFFDIG2 < 1 */

    yesrel = TRUE_;
    offdig = 0.;
    tmp = template_blas_sqrt((absMACRO(d__[1])));
    if (tmp < rmin) {
	yesrel = FALSE_;
    }
    if (! yesrel) {
	goto L11;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	tmp2 = template_blas_sqrt((d__1 = d__[i__], absMACRO(d__1)));
	if (tmp2 < rmin) {
	    yesrel = FALSE_;
	}
	if (! yesrel) {
	    goto L11;
	}
	offdig2 = (d__1 = e[i__ - 1], absMACRO(d__1)) / (tmp * tmp2);
	if (offdig + offdig2 >= .999) {
	    yesrel = FALSE_;
	}
	if (! yesrel) {
	    goto L11;
	}
	tmp = tmp2;
	offdig = offdig2;
/* L10: */
    }
L11:
    if (yesrel) {
	*info = 0;
	return 0;
    } else {
    }


/*     *** MORE TO BE IMPLEMENTED *** */


/*     Test if the lower bidiagonal matrix L from T = L D L^T */
/*     (zero shift facto) is well conditioned */


/*     Test if the upper bidiagonal matrix U from T = U D U^T */
/*     (zero shift facto) is well conditioned. */
/*     In this case, the matrix needs to be flipped and, at the end */
/*     of the eigenvector computation, the flip needs to be applied */
/*     to the computed eigenvectors (and the support) */


    return 0;

/*     END OF DLARRR */

} /* dlarrr_ */

#endif
