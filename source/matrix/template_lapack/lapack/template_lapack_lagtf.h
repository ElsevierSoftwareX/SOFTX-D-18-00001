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
 

#ifndef TEMPLATE_LAPACK_LAGTF_HEADER
#define TEMPLATE_LAPACK_LAGTF_HEADER


template<class Treal>
int template_lapack_lagtf(const integer *n, Treal *a, const Treal *lambda, 
	Treal *b, Treal *c__, const Treal *tol, Treal *d__, 
	integer *in, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n   
    tridiagonal matrix and lambda is a scalar, as   

       T - lambda*I = PLU,   

    where P is a permutation matrix, L is a unit lower tridiagonal matrix   
    with at most one non-zero sub-diagonal elements per column and U is   
    an upper triangular matrix with at most two non-zero super-diagonal   
    elements per column.   

    The factorization is obtained by Gaussian elimination with partial   
    pivoting and implicit row scaling.   

    The parameter LAMBDA is included in the routine so that DLAGTF may   
    be used, in conjunction with DLAGTS, to obtain eigenvectors of T by   
    inverse iteration.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix T.   

    A       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, A must contain the diagonal elements of T.   

            On exit, A is overwritten by the n diagonal elements of the   
            upper triangular matrix U of the factorization of T.   

    LAMBDA  (input) DOUBLE PRECISION   
            On entry, the scalar lambda.   

    B       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, B must contain the (n-1) super-diagonal elements of   
            T.   

            On exit, B is overwritten by the (n-1) super-diagonal   
            elements of the matrix U of the factorization of T.   

    C       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, C must contain the (n-1) sub-diagonal elements of   
            T.   

            On exit, C is overwritten by the (n-1) sub-diagonal elements   
            of the matrix L of the factorization of T.   

    TOL     (input) DOUBLE PRECISION   
            On entry, a relative tolerance used to indicate whether or   
            not the matrix (T - lambda*I) is nearly singular. TOL should   
            normally be chose as approximately the largest relative error   
            in the elements of T. For example, if the elements of T are   
            correct to about 4 significant figures, then TOL should be   
            set to about 5*10**(-4). If TOL is supplied as less than eps,   
            where eps is the relative machine precision, then the value   
            eps is used in place of TOL.   

    D       (output) DOUBLE PRECISION array, dimension (N-2)   
            On exit, D is overwritten by the (n-2) second super-diagonal   
            elements of the matrix U of the factorization of T.   

    IN      (output) INTEGER array, dimension (N)   
            On exit, IN contains details of the permutation matrix P. If   
            an interchange occurred at the kth step of the elimination,   
            then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)   
            returns the smallest positive integer j such that   

               abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,   

            where norm( A(j) ) denotes the sum of the absolute values of   
            the jth row of the matrix A. If no such j exists then IN(n)   
            is returned as zero. If IN(n) is returned as positive, then a   
            diagonal element of U is small, indicating that   
            (T - lambda*I) is singular or nearly singular,   

    INFO    (output) INTEGER   
            = 0   : successful exit   
            .lt. 0: if INFO = -k, the kth argument had an illegal value   

   =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer i__1;
    Treal d__1, d__2;
    /* Local variables */
     Treal temp, mult;
     integer k;
     Treal scale1, scale2;
     Treal tl;
     Treal eps, piv1, piv2;

    --in;
    --d__;
    --c__;
    --b;
    --a;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	template_blas_erbla("LAGTF ", &i__1);
	return 0;
    }

    if (*n == 0) {
	return 0;
    }

    a[1] -= *lambda;
    in[*n] = 0;
    if (*n == 1) {
	if (a[1] == 0.) {
	    in[1] = 1;
	}
	return 0;
    }

    eps = template_lapack_lamch("Epsilon", (Treal)0);

    tl = maxMACRO(*tol,eps);
    scale1 = absMACRO(a[1]) + absMACRO(b[1]);
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	a[k + 1] -= *lambda;
	scale2 = (d__1 = c__[k], absMACRO(d__1)) + (d__2 = a[k + 1], absMACRO(d__2));
	if (k < *n - 1) {
	    scale2 += (d__1 = b[k + 1], absMACRO(d__1));
	}
	if (a[k] == 0.) {
	    piv1 = 0.;
	} else {
	    piv1 = (d__1 = a[k], absMACRO(d__1)) / scale1;
	}
	if (c__[k] == 0.) {
	    in[k] = 0;
	    piv2 = 0.;
	    scale1 = scale2;
	    if (k < *n - 1) {
		d__[k] = 0.;
	    }
	} else {
	    piv2 = (d__1 = c__[k], absMACRO(d__1)) / scale2;
	    if (piv2 <= piv1) {
		in[k] = 0;
		scale1 = scale2;
		c__[k] /= a[k];
		a[k + 1] -= c__[k] * b[k];
		if (k < *n - 1) {
		    d__[k] = 0.;
		}
	    } else {
		in[k] = 1;
		mult = a[k] / c__[k];
		a[k] = c__[k];
		temp = a[k + 1];
		a[k + 1] = b[k] - mult * temp;
		if (k < *n - 1) {
		    d__[k] = b[k + 1];
		    b[k + 1] = -mult * d__[k];
		}
		b[k] = temp;
		c__[k] = mult;
	    }
	}
	if (maxMACRO(piv1,piv2) <= tl && in[*n] == 0) {
	    in[*n] = k;
	}
/* L10: */
    }
    if ((d__1 = a[*n], absMACRO(d__1)) <= scale1 * tl && in[*n] == 0) {
	in[*n] = *n;
    }

    return 0;

/*     End of DLAGTF */

} /* dlagtf_ */

#endif
