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
 

#ifndef TEMPLATE_LAPACK_LASV2_HEADER
#define TEMPLATE_LAPACK_LASV2_HEADER


template<class Treal>
int template_lapack_lasv2(const Treal *f, const Treal *g, const Treal *h__, 
	Treal *ssmin, Treal *ssmax, Treal *snr, Treal *
	csr, Treal *snl, Treal *csl)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASV2 computes the singular value decomposition of a 2-by-2   
    triangular matrix   
       [  F   G  ]   
       [  0   H  ].   
    On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the   
    smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and   
    right singular vectors for abs(SSMAX), giving the decomposition   

       [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]   
       [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].   

    Arguments   
    =========   

    F       (input) DOUBLE PRECISION   
            The (1,1) element of the 2-by-2 matrix.   

    G       (input) DOUBLE PRECISION   
            The (1,2) element of the 2-by-2 matrix.   

    H       (input) DOUBLE PRECISION   
            The (2,2) element of the 2-by-2 matrix.   

    SSMIN   (output) DOUBLE PRECISION   
            abs(SSMIN) is the smaller singular value.   

    SSMAX   (output) DOUBLE PRECISION   
            abs(SSMAX) is the larger singular value.   

    SNL     (output) DOUBLE PRECISION   
    CSL     (output) DOUBLE PRECISION   
            The vector (CSL, SNL) is a unit left singular vector for the   
            singular value abs(SSMAX).   

    SNR     (output) DOUBLE PRECISION   
    CSR     (output) DOUBLE PRECISION   
            The vector (CSR, SNR) is a unit right singular vector for the   
            singular value abs(SSMAX).   

    Further Details   
    ===============   

    Any input parameter may be aliased with any output parameter.   

    Barring over/underflow and assuming a guard digit in subtraction, all   
    output quantities are correct to within a few units in the last   
    place (ulps).   

    In IEEE arithmetic, the code works correctly if one matrix element is   
    infinite.   

    Overflow will not occur unless the largest singular value itself   
    overflows or is within a few ulps of overflow. (On machines with   
    partial overflow, like the Cray, overflow may occur if the largest   
    singular value is within a factor of 2 of overflow.)   

    Underflow is harmless if underflow is gradual. Otherwise, results   
    may correspond to a matrix modified by perturbations of size near   
    the underflow threshold.   

   ===================================================================== */
    /* Table of constant values */
     Treal c_b3 = 2.;
     Treal c_b4 = 1.;
    
    /* System generated locals */
    Treal d__1;
    /* Local variables */
     integer pmax;
     Treal temp;
     logical swap;
     Treal a, d__, l, m, r__, s, t, tsign, fa, ga, ha;
     Treal ft, gt, ht, mm;
     logical gasmal;
     Treal tt, clt, crt, slt, srt;

     /* Initialization added by Elias to get rid of compiler warnings. */
     tsign = 0;


    ft = *f;
    fa = absMACRO(ft);
    ht = *h__;
    ha = absMACRO(*h__);

/*     PMAX points to the maximum absolute element of matrix   
         PMAX = 1 if F largest in absolute values   
         PMAX = 2 if G largest in absolute values   
         PMAX = 3 if H largest in absolute values */

    pmax = 1;
    swap = ha > fa;
    if (swap) {
	pmax = 3;
	temp = ft;
	ft = ht;
	ht = temp;
	temp = fa;
	fa = ha;
	ha = temp;

/*        Now FA .ge. HA */

    }
    gt = *g;
    ga = absMACRO(gt);
    if (ga == 0.) {

/*        Diagonal matrix */

	*ssmin = ha;
	*ssmax = fa;
	clt = 1.;
	crt = 1.;
	slt = 0.;
	srt = 0.;
    } else {
	gasmal = TRUE_;
	if (ga > fa) {
	    pmax = 2;
	    if (fa / ga < template_lapack_lamch("EPS", (Treal)0)) {

/*              Case of very large GA */

		gasmal = FALSE_;
		*ssmax = ga;
		if (ha > 1.) {
		    *ssmin = fa / (ga / ha);
		} else {
		    *ssmin = fa / ga * ha;
		}
		clt = 1.;
		slt = ht / gt;
		srt = 1.;
		crt = ft / gt;
	    }
	}
	if (gasmal) {

/*           Normal case */

	    d__ = fa - ha;
	    if (d__ == fa) {

/*              Copes with infinite F or H */

		l = 1.;
	    } else {
		l = d__ / fa;
	    }

/*           Note that 0 .le. L .le. 1 */

	    m = gt / ft;

/*           Note that abs(M) .le. 1/macheps */

	    t = 2. - l;

/*           Note that T .ge. 1 */

	    mm = m * m;
	    tt = t * t;
	    s = template_blas_sqrt(tt + mm);

/*           Note that 1 .le. S .le. 1 + 1/macheps */

	    if (l == 0.) {
		r__ = absMACRO(m);
	    } else {
		r__ = template_blas_sqrt(l * l + mm);
	    }

/*           Note that 0 .le. R .le. 1 + 1/macheps */

	    a = (s + r__) * .5;

/*           Note that 1 .le. A .le. 1 + abs(M) */

	    *ssmin = ha / a;
	    *ssmax = fa * a;
	    if (mm == 0.) {

/*              Note that M is very tiny */

		if (l == 0.) {
		    t = template_lapack_d_sign(&c_b3, &ft) * template_lapack_d_sign(&c_b4, &gt);
		} else {
		    t = gt / template_lapack_d_sign(&d__, &ft) + m / t;
		}
	    } else {
		t = (m / (s + t) + m / (r__ + l)) * (a + 1.);
	    }
	    l = template_blas_sqrt(t * t + 4.);
	    crt = 2. / l;
	    srt = t / l;
	    clt = (crt + srt * m) / a;
	    slt = ht / ft * srt / a;
	}
    }
    if (swap) {
	*csl = srt;
	*snl = crt;
	*csr = slt;
	*snr = clt;
    } else {
	*csl = clt;
	*snl = slt;
	*csr = crt;
	*snr = srt;
    }

/*     Correct signs of SSMAX and SSMIN */

    if (pmax == 1) {
	tsign = template_lapack_d_sign(&c_b4, csr) * template_lapack_d_sign(&c_b4, csl) * template_lapack_d_sign(&c_b4, f);
    }
    if (pmax == 2) {
	tsign = template_lapack_d_sign(&c_b4, snr) * template_lapack_d_sign(&c_b4, csl) * template_lapack_d_sign(&c_b4, g);
    }
    if (pmax == 3) {
	tsign = template_lapack_d_sign(&c_b4, snr) * template_lapack_d_sign(&c_b4, snl) * template_lapack_d_sign(&c_b4, h__);
    }
    *ssmax = template_lapack_d_sign(ssmax, &tsign);
    d__1 = tsign * template_lapack_d_sign(&c_b4, f) * template_lapack_d_sign(&c_b4, h__);
    *ssmin = template_lapack_d_sign(ssmin, &d__1);
    return 0;

/*     End of DLASV2 */

} /* dlasv2_ */

#endif
