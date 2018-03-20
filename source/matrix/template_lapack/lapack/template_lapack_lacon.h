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
 

#ifndef TEMPLATE_LAPACK_LACON_HEADER
#define TEMPLATE_LAPACK_LACON_HEADER


template<class Treal>
int template_lapack_lacon(const integer *n, Treal *v, Treal *x, 
	integer *isgn, Treal *est, integer *kase)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLACON estimates the 1-norm of a square, real matrix A.   
    Reverse communication is used for evaluating matrix-vector products.   

    Arguments   
    =========   

    N      (input) INTEGER   
           The order of the matrix.  N >= 1.   

    V      (workspace) DOUBLE PRECISION array, dimension (N)   
           On the final return, V = A*W,  where  EST = norm(V)/norm(W)   
           (W is not returned).   

    X      (input/output) DOUBLE PRECISION array, dimension (N)   
           On an intermediate return, X should be overwritten by   
                 A * X,   if KASE=1,   
                 A' * X,  if KASE=2,   
           and DLACON must be re-called with all the other parameters   
           unchanged.   

    ISGN   (workspace) INTEGER array, dimension (N)   

    EST    (output) DOUBLE PRECISION   
           An estimate (a lower bound) for norm(A).   

    KASE   (input/output) INTEGER   
           On the initial call to DLACON, KASE should be 0.   
           On an intermediate return, KASE will be 1 or 2, indicating   
           whether X should be overwritten by A * X  or A' * X.   
           On the final return from DLACON, KASE will again be 0.   

    Further Details   
    ======= =======   

    Contributed by Nick Higham, University of Manchester.   
    Originally named SONEST, dated March 16, 1988.   

    Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of   
    a real or complex matrix, with applications to condition estimation",   
    ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.   

    =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     Treal c_b11 = 1.;
    
    /* System generated locals */
    integer i__1;
    Treal d__1;
    /* Builtin functions */
    double d_sign(Treal *, Treal *);
    integer i_dnnt(Treal *);
    /* Local variables */
     integer iter;
     Treal temp;
     integer jump, i__, j;
     integer jlast;
     Treal altsgn, estold;


    --isgn;
    --x;
    --v;

    /* Function Body */
    if (*kase == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = 1. / (Treal) (*n);
/* L10: */
	}
	*kase = 1;
	jump = 1;
	return 0;
    }

    switch (jump) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L110;
	case 5:  goto L140;
    }

/*     ................ ENTRY   (JUMP = 1)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

L20:
    if (*n == 1) {
	v[1] = x[1];
	*est = absMACRO(v[1]);
/*        ... QUIT */
	goto L150;
    }
    *est = template_blas_asum(n, &x[1], &c__1);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = d_sign(&c_b11, &x[i__]);
	isgn[i__] = i_dnnt(&x[i__]);
/* L30: */
    }
    *kase = 2;
    jump = 2;
    return 0;

/*     ................ ENTRY   (JUMP = 2)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

L40:
    j = template_blas_idamax(n, &x[1], &c__1);
    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

L50:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
/* L60: */
    }
    x[j] = 1.;
    *kase = 1;
    jump = 3;
    return 0;

/*     ................ ENTRY   (JUMP = 3)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L70:
    template_blas_copy(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = template_blas_asum(n, &v[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = d_sign(&c_b11, &x[i__]);
	if (i_dnnt(&d__1) != isgn[i__]) {
	    goto L90;
	}
/* L80: */
    }
/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    goto L120;

L90:
/*     TEST FOR CYCLING. */
    if (*est <= estold) {
	goto L120;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = d_sign(&c_b11, &x[i__]);
	isgn[i__] = i_dnnt(&x[i__]);
/* L100: */
    }
    *kase = 2;
    jump = 4;
    return 0;

/*     ................ ENTRY   (JUMP = 4)   
       X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

L110:
    jlast = j;
    j = template_blas_idamax(n, &x[1], &c__1);
    if (x[jlast] != (d__1 = x[j], absMACRO(d__1)) && iter < 5) {
	++iter;
	goto L50;
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

L120:
    altsgn = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = altsgn * ((Treal) (i__ - 1) / (Treal) (*n - 1) + 
		1.);
	altsgn = -altsgn;
/* L130: */
    }
    *kase = 1;
    jump = 5;
    return 0;

/*     ................ ENTRY   (JUMP = 5)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L140:
    temp = template_blas_asum(n, &x[1], &c__1) / (Treal) (*n * 3) * 2.;
    if (temp > *est) {
	template_blas_copy(n, &x[1], &c__1, &v[1], &c__1);
	*est = temp;
    }

L150:
    *kase = 0;
    return 0;

/*     End of DLACON */

} /* dlacon_ */

#endif
