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
 

#ifndef TEMPLATE_LAPACK_LALN2_HEADER
#define TEMPLATE_LAPACK_LALN2_HEADER


template<class Treal>
int template_lapack_laln2(const logical *ltrans, const integer *na, const integer *nw, 
	const Treal *smin, const Treal *ca, const Treal *a, const integer *lda, 
	const Treal *d1, const Treal *d2, const Treal *b, const integer *ldb, 
	const Treal *wr, const Treal *wi, Treal *x, const integer *ldx, 
	Treal *scale, Treal *xnorm, integer *info)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLALN2 solves a system of the form  (ca A - w D ) X = s B   
    or (ca A' - w D) X = s B   with possible scaling ("s") and   
    perturbation of A.  (A' means A-transpose.)   

    A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA   
    real diagonal matrix, w is a real or complex value, and X and B are   
    NA x 1 matrices -- real if w is real, complex if w is complex.  NA   
    may be 1 or 2.   

    If w is complex, X and B are represented as NA x 2 matrices,   
    the first column of each being the real part and the second   
    being the imaginary part.   

    "s" is a scaling factor (.LE. 1), computed by DLALN2, which is   
    so chosen that X can be computed without overflow.  X is further   
    scaled if necessary to assure that norm(ca A - w D)*norm(X) is less   
    than overflow.   

    If both singular values of (ca A - w D) are less than SMIN,   
    SMIN*identity will be used instead of (ca A - w D).  If only one   
    singular value is less than SMIN, one element of (ca A - w D) will be   
    perturbed enough to make the smallest singular value roughly SMIN.   
    If both singular values are at least SMIN, (ca A - w D) will not be   
    perturbed.  In any case, the perturbation will be at most some small   
    multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values   
    are computed by infinity-norm approximations, and thus will only be   
    correct to a factor of 2 or so.   

    Note: all input quantities are assumed to be smaller than overflow   
    by a reasonable factor.  (See BIGNUM.)   

    Arguments   
    ==========   

    LTRANS  (input) LOGICAL   
            =.TRUE.:  A-transpose will be used.   
            =.FALSE.: A will be used (not transposed.)   

    NA      (input) INTEGER   
            The size of the matrix A.  It may (only) be 1 or 2.   

    NW      (input) INTEGER   
            1 if "w" is real, 2 if "w" is complex.  It may only be 1   
            or 2.   

    SMIN    (input) DOUBLE PRECISION   
            The desired lower bound on the singular values of A.  This   
            should be a safe distance away from underflow or overflow,   
            say, between (underflow/machine precision) and  (machine   
            precision * overflow ).  (See BIGNUM and ULP.)   

    CA      (input) DOUBLE PRECISION   
            The coefficient c, which A is multiplied by.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,NA)   
            The NA x NA matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of A.  It must be at least NA.   

    D1      (input) DOUBLE PRECISION   
            The 1,1 element in the diagonal matrix D.   

    D2      (input) DOUBLE PRECISION   
            The 2,2 element in the diagonal matrix D.  Not used if NW=1.   

    B       (input) DOUBLE PRECISION array, dimension (LDB,NW)   
            The NA x NW matrix B (right-hand side).  If NW=2 ("w" is   
            complex), column 1 contains the real part of B and column 2   
            contains the imaginary part.   

    LDB     (input) INTEGER   
            The leading dimension of B.  It must be at least NA.   

    WR      (input) DOUBLE PRECISION   
            The real part of the scalar "w".   

    WI      (input) DOUBLE PRECISION   
            The imaginary part of the scalar "w".  Not used if NW=1.   

    X       (output) DOUBLE PRECISION array, dimension (LDX,NW)   
            The NA x NW matrix X (unknowns), as computed by DLALN2.   
            If NW=2 ("w" is complex), on exit, column 1 will contain   
            the real part of X and column 2 will contain the imaginary   
            part.   

    LDX     (input) INTEGER   
            The leading dimension of X.  It must be at least NA.   

    SCALE   (output) DOUBLE PRECISION   
            The scale factor that B must be multiplied by to insure   
            that overflow does not occur when computing X.  Thus,   
            (ca A - w D) X  will be SCALE*B, not B (ignoring   
            perturbations of A.)  It will be at most 1.   

    XNORM   (output) DOUBLE PRECISION   
            The infinity-norm of X, when X is regarded as an NA x NW   
            real matrix.   

    INFO    (output) INTEGER   
            An error flag.  It will be set to zero if no error occurs,   
            a negative number if an argument is in error, or a positive   
            number if  ca A - w D  had to be perturbed.   
            The possible values are:   
            = 0: No error occurred, and (ca A - w D) did not have to be   
                   perturbed.   
            = 1: (ca A - w D) had to be perturbed to make its smallest   
                 (or only) singular value greater than SMIN.   
            NOTE: In the interests of speed, this routine does not   
                  check the inputs for errors.   

   =====================================================================   

       Parameter adjustments */
    /* Initialized data */
     logical zswap[4] = { FALSE_,FALSE_,TRUE_,TRUE_ };
     logical rswap[4] = { FALSE_,TRUE_,FALSE_,TRUE_ };
     integer ipivot[16]	/* was [4][4] */ = { 1,2,3,4,2,1,4,3,3,4,1,2,
	    4,3,2,1 };
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset;
    Treal d__1, d__2, d__3, d__4, d__5, d__6;
     Treal equiv_0[4], equiv_1[4];
    /* Local variables */
     Treal bbnd, cmax, ui11r, ui12s, temp, ur11r, ur12s;
     integer j;
     Treal u22abs;
     integer icmax;
     Treal bnorm, cnorm, smini;
#define ci (equiv_0)
#define cr (equiv_1)
     Treal bignum, bi1, bi2, br1, br2, smlnum, xi1, xi2, xr1, xr2, 
	    ci21, ci22, cr21, cr22, li21, csi, ui11, lr21, ui12, ui22;
#define civ (equiv_0)
     Treal csr, ur11, ur12, ur22;
#define crv (equiv_1)
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define x_ref(a_1,a_2) x[(a_2)*x_dim1 + a_1]
#define ci_ref(a_1,a_2) ci[(a_2)*2 + a_1 - 3]
#define cr_ref(a_1,a_2) cr[(a_2)*2 + a_1 - 3]
#define ipivot_ref(a_1,a_2) ipivot[(a_2)*4 + a_1 - 5]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1 * 1;
    x -= x_offset;

    /* Function Body   

       Compute BIGNUM */

    smlnum = 2. * template_lapack_lamch("Safe minimum", (Treal)0);
    bignum = 1. / smlnum;
    smini = maxMACRO(*smin,smlnum);

/*     Don't check for input errors */

    *info = 0;

/*     Standard Initializations */

    *scale = 1.;

    if (*na == 1) {

/*        1 x 1  (i.e., scalar) system   C X = B */

	if (*nw == 1) {

/*           Real 1x1 system.   

             C = ca A - w D */

	    csr = *ca * a_ref(1, 1) - *wr * *d1;
	    cnorm = absMACRO(csr);

/*           If | C | < SMINI, use C = SMINI */

	    if (cnorm < smini) {
		csr = smini;
		cnorm = smini;
		*info = 1;
	    }

/*           Check scaling for  X = B / C */

	    bnorm = (d__1 = b_ref(1, 1), absMACRO(d__1));
	    if (cnorm < 1. && bnorm > 1.) {
		if (bnorm > bignum * cnorm) {
		    *scale = 1. / bnorm;
		}
	    }

/*           Compute X */

	    x_ref(1, 1) = b_ref(1, 1) * *scale / csr;
	    *xnorm = (d__1 = x_ref(1, 1), absMACRO(d__1));
	} else {

/*           Complex 1x1 system (w is complex)   

             C = ca A - w D */

	    csr = *ca * a_ref(1, 1) - *wr * *d1;
	    csi = -(*wi) * *d1;
	    cnorm = absMACRO(csr) + absMACRO(csi);

/*           If | C | < SMINI, use C = SMINI */

	    if (cnorm < smini) {
		csr = smini;
		csi = 0.;
		cnorm = smini;
		*info = 1;
	    }

/*           Check scaling for  X = B / C */

	    bnorm = (d__1 = b_ref(1, 1), absMACRO(d__1)) + (d__2 = b_ref(1, 2), 
		    absMACRO(d__2));
	    if (cnorm < 1. && bnorm > 1.) {
		if (bnorm > bignum * cnorm) {
		    *scale = 1. / bnorm;
		}
	    }

/*           Compute X */

	    d__1 = *scale * b_ref(1, 1);
	    d__2 = *scale * b_ref(1, 2);
	    template_lapack_ladiv(&d__1, &d__2, &csr, &csi, &x_ref(1, 1), &x_ref(1, 2));
	    *xnorm = (d__1 = x_ref(1, 1), absMACRO(d__1)) + (d__2 = x_ref(1, 2), 
		    absMACRO(d__2));
	}

    } else {

/*        2x2 System   

          Compute the real part of  C = ca A - w D  (or  ca A' - w D ) */

	cr_ref(1, 1) = *ca * a_ref(1, 1) - *wr * *d1;
	cr_ref(2, 2) = *ca * a_ref(2, 2) - *wr * *d2;
	if (*ltrans) {
	    cr_ref(1, 2) = *ca * a_ref(2, 1);
	    cr_ref(2, 1) = *ca * a_ref(1, 2);
	} else {
	    cr_ref(2, 1) = *ca * a_ref(2, 1);
	    cr_ref(1, 2) = *ca * a_ref(1, 2);
	}

	if (*nw == 1) {

/*           Real 2x2 system  (w is real)   

             Find the largest element in C */

	    cmax = 0.;
	    icmax = 0;

	    for (j = 1; j <= 4; ++j) {
		if ((d__1 = crv[j - 1], absMACRO(d__1)) > cmax) {
		    cmax = (d__1 = crv[j - 1], absMACRO(d__1));
		    icmax = j;
		}
/* L10: */
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

	    if (cmax < smini) {
/* Computing MAX */
		d__3 = (d__1 = b_ref(1, 1), absMACRO(d__1)), d__4 = (d__2 = b_ref(
			2, 1), absMACRO(d__2));
		bnorm = maxMACRO(d__3,d__4);
		if (smini < 1. && bnorm > 1.) {
		    if (bnorm > bignum * smini) {
			*scale = 1. / bnorm;
		    }
		}
		temp = *scale / smini;
		x_ref(1, 1) = temp * b_ref(1, 1);
		x_ref(2, 1) = temp * b_ref(2, 1);
		*xnorm = temp * bnorm;
		*info = 1;
		return 0;
	    }

/*           Gaussian elimination with complete pivoting. */

	    ur11 = crv[icmax - 1];
	    cr21 = crv[ipivot_ref(2, icmax) - 1];
	    ur12 = crv[ipivot_ref(3, icmax) - 1];
	    cr22 = crv[ipivot_ref(4, icmax) - 1];
	    ur11r = 1. / ur11;
	    lr21 = ur11r * cr21;
	    ur22 = cr22 - ur12 * lr21;

/*           If smaller pivot < SMINI, use SMINI */

	    if (absMACRO(ur22) < smini) {
		ur22 = smini;
		*info = 1;
	    }
	    if (rswap[icmax - 1]) {
		br1 = b_ref(2, 1);
		br2 = b_ref(1, 1);
	    } else {
		br1 = b_ref(1, 1);
		br2 = b_ref(2, 1);
	    }
	    br2 -= lr21 * br1;
/* Computing MAX */
	    d__2 = (d__1 = br1 * (ur22 * ur11r), absMACRO(d__1)), d__3 = absMACRO(br2);
	    bbnd = maxMACRO(d__2,d__3);
	    if (bbnd > 1. && absMACRO(ur22) < 1.) {
		if (bbnd >= bignum * absMACRO(ur22)) {
		    *scale = 1. / bbnd;
		}
	    }

	    xr2 = br2 * *scale / ur22;
	    xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
	    if (zswap[icmax - 1]) {
		x_ref(1, 1) = xr2;
		x_ref(2, 1) = xr1;
	    } else {
		x_ref(1, 1) = xr1;
		x_ref(2, 1) = xr2;
	    }
/* Computing MAX */
	    d__1 = absMACRO(xr1), d__2 = absMACRO(xr2);
	    *xnorm = maxMACRO(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

	    if (*xnorm > 1. && cmax > 1.) {
		if (*xnorm > bignum / cmax) {
		    temp = cmax / bignum;
		    x_ref(1, 1) = temp * x_ref(1, 1);
		    x_ref(2, 1) = temp * x_ref(2, 1);
		    *xnorm = temp * *xnorm;
		    *scale = temp * *scale;
		}
	    }
	} else {

/*           Complex 2x2 system  (w is complex)   

             Find the largest element in C */

	    ci_ref(1, 1) = -(*wi) * *d1;
	    ci_ref(2, 1) = 0.;
	    ci_ref(1, 2) = 0.;
	    ci_ref(2, 2) = -(*wi) * *d2;
	    cmax = 0.;
	    icmax = 0;

	    for (j = 1; j <= 4; ++j) {
		if ((d__1 = crv[j - 1], absMACRO(d__1)) + (d__2 = civ[j - 1], absMACRO(
			d__2)) > cmax) {
		    cmax = (d__1 = crv[j - 1], absMACRO(d__1)) + (d__2 = civ[j - 1]
			    , absMACRO(d__2));
		    icmax = j;
		}
/* L20: */
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

	    if (cmax < smini) {
/* Computing MAX */
		d__5 = (d__1 = b_ref(1, 1), absMACRO(d__1)) + (d__2 = b_ref(1, 2), 
			absMACRO(d__2)), d__6 = (d__3 = b_ref(2, 1), absMACRO(d__3)) + (
			d__4 = b_ref(2, 2), absMACRO(d__4));
		bnorm = maxMACRO(d__5,d__6);
		if (smini < 1. && bnorm > 1.) {
		    if (bnorm > bignum * smini) {
			*scale = 1. / bnorm;
		    }
		}
		temp = *scale / smini;
		x_ref(1, 1) = temp * b_ref(1, 1);
		x_ref(2, 1) = temp * b_ref(2, 1);
		x_ref(1, 2) = temp * b_ref(1, 2);
		x_ref(2, 2) = temp * b_ref(2, 2);
		*xnorm = temp * bnorm;
		*info = 1;
		return 0;
	    }

/*           Gaussian elimination with complete pivoting. */

	    ur11 = crv[icmax - 1];
	    ui11 = civ[icmax - 1];
	    cr21 = crv[ipivot_ref(2, icmax) - 1];
	    ci21 = civ[ipivot_ref(2, icmax) - 1];
	    ur12 = crv[ipivot_ref(3, icmax) - 1];
	    ui12 = civ[ipivot_ref(3, icmax) - 1];
	    cr22 = crv[ipivot_ref(4, icmax) - 1];
	    ci22 = civ[ipivot_ref(4, icmax) - 1];
	    if (icmax == 1 || icmax == 4) {

/*              Code when off-diagonals of pivoted C are real */

		if (absMACRO(ur11) > absMACRO(ui11)) {
		    temp = ui11 / ur11;
/* Computing 2nd power */
		    d__1 = temp;
		    ur11r = 1. / (ur11 * (d__1 * d__1 + 1.));
		    ui11r = -temp * ur11r;
		} else {
		    temp = ur11 / ui11;
/* Computing 2nd power */
		    d__1 = temp;
		    ui11r = -1. / (ui11 * (d__1 * d__1 + 1.));
		    ur11r = -temp * ui11r;
		}
		lr21 = cr21 * ur11r;
		li21 = cr21 * ui11r;
		ur12s = ur12 * ur11r;
		ui12s = ur12 * ui11r;
		ur22 = cr22 - ur12 * lr21;
		ui22 = ci22 - ur12 * li21;
	    } else {

/*              Code when diagonals of pivoted C are real */

		ur11r = 1. / ur11;
		ui11r = 0.;
		lr21 = cr21 * ur11r;
		li21 = ci21 * ur11r;
		ur12s = ur12 * ur11r;
		ui12s = ui12 * ur11r;
		ur22 = cr22 - ur12 * lr21 + ui12 * li21;
		ui22 = -ur12 * li21 - ui12 * lr21;
	    }
	    u22abs = absMACRO(ur22) + absMACRO(ui22);

/*           If smaller pivot < SMINI, use SMINI */

	    if (u22abs < smini) {
		ur22 = smini;
		ui22 = 0.;
		*info = 1;
	    }
	    if (rswap[icmax - 1]) {
		br2 = b_ref(1, 1);
		br1 = b_ref(2, 1);
		bi2 = b_ref(1, 2);
		bi1 = b_ref(2, 2);
	    } else {
		br1 = b_ref(1, 1);
		br2 = b_ref(2, 1);
		bi1 = b_ref(1, 2);
		bi2 = b_ref(2, 2);
	    }
	    br2 = br2 - lr21 * br1 + li21 * bi1;
	    bi2 = bi2 - li21 * br1 - lr21 * bi1;
/* Computing MAX */
	    d__1 = (absMACRO(br1) + absMACRO(bi1)) * (u22abs * (absMACRO(ur11r) + absMACRO(ui11r))
		    ), d__2 = absMACRO(br2) + absMACRO(bi2);
	    bbnd = maxMACRO(d__1,d__2);
	    if (bbnd > 1. && u22abs < 1.) {
		if (bbnd >= bignum * u22abs) {
		    *scale = 1. / bbnd;
		    br1 = *scale * br1;
		    bi1 = *scale * bi1;
		    br2 = *scale * br2;
		    bi2 = *scale * bi2;
		}
	    }

	    template_lapack_ladiv(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
	    xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
	    xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
	    if (zswap[icmax - 1]) {
		x_ref(1, 1) = xr2;
		x_ref(2, 1) = xr1;
		x_ref(1, 2) = xi2;
		x_ref(2, 2) = xi1;
	    } else {
		x_ref(1, 1) = xr1;
		x_ref(2, 1) = xr2;
		x_ref(1, 2) = xi1;
		x_ref(2, 2) = xi2;
	    }
/* Computing MAX */
	    d__1 = absMACRO(xr1) + absMACRO(xi1), d__2 = absMACRO(xr2) + absMACRO(xi2);
	    *xnorm = maxMACRO(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

	    if (*xnorm > 1. && cmax > 1.) {
		if (*xnorm > bignum / cmax) {
		    temp = cmax / bignum;
		    x_ref(1, 1) = temp * x_ref(1, 1);
		    x_ref(2, 1) = temp * x_ref(2, 1);
		    x_ref(1, 2) = temp * x_ref(1, 2);
		    x_ref(2, 2) = temp * x_ref(2, 2);
		    *xnorm = temp * *xnorm;
		    *scale = temp * *scale;
		}
	    }
	}
    }

    return 0;

/*     End of DLALN2 */

} /* dlaln2_ */

#undef ipivot_ref
#undef cr_ref
#undef ci_ref
#undef x_ref
#undef b_ref
#undef a_ref
#undef crv
#undef civ
#undef cr
#undef ci


#endif
