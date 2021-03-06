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
 

#ifndef TEMPLATE_LAPACK_LATRS_HEADER
#define TEMPLATE_LAPACK_LATRS_HEADER


template<class Treal>
int template_lapack_latrs(const char *uplo, const char *trans, const char *diag, const char *
	normin, const integer *n, const Treal *a, const integer *lda, Treal *x, 
	Treal *scale, Treal *cnorm, integer *info)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1992   


    Purpose   
    =======   

    DLATRS solves one of the triangular systems   

       A *x = s*b  or  A'*x = s*b   

    with scaling to prevent overflow.  Here A is an upper or lower   
    triangular matrix, A' denotes the transpose of A, x and b are   
    n-element vectors, and s is a scaling factor, usually less than   
    or equal to 1, chosen so that the components of x will be less than   
    the overflow threshold.  If the unscaled problem will not cause   
    overflow, the Level 2 BLAS routine DTRSV is called.  If the matrix A   
    is singular (A(j,j) = 0 for some j), then s is set to 0 and a   
    non-trivial solution to A*x = 0 is returned.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular.   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    TRANS   (input) CHARACTER*1   
            Specifies the operation applied to A.   
            = 'N':  Solve A * x = s*b  (No transpose)   
            = 'T':  Solve A'* x = s*b  (Transpose)   
            = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    NORMIN  (input) CHARACTER*1   
            Specifies whether CNORM has been set or not.   
            = 'Y':  CNORM contains the column norms on entry   
            = 'N':  CNORM is not set on entry.  On exit, the norms will   
                    be computed and stored in CNORM.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The triangular matrix A.  If UPLO = 'U', the leading n by n   
            upper triangular part of the array A contains the upper   
            triangular matrix, and the strictly lower triangular part of   
            A is not referenced.  If UPLO = 'L', the leading n by n lower   
            triangular part of the array A contains the lower triangular   
            matrix, and the strictly upper triangular part of A is not   
            referenced.  If DIAG = 'U', the diagonal elements of A are   
            also not referenced and are assumed to be 1.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max (1,N).   

    X       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the right hand side b of the triangular system.   
            On exit, X is overwritten by the solution vector x.   

    SCALE   (output) DOUBLE PRECISION   
            The scaling factor s for the triangular system   
               A * x = s*b  or  A'* x = s*b.   
            If SCALE = 0, the matrix A is singular or badly scaled, and   
            the vector x is an exact or approximate solution to A*x = 0.   

    CNORM   (input or output) DOUBLE PRECISION array, dimension (N)   

            If NORMIN = 'Y', CNORM is an input argument and CNORM(j)   
            contains the norm of the off-diagonal part of the j-th column   
            of A.  If TRANS = 'N', CNORM(j) must be greater than or equal   
            to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)   
            must be greater than or equal to the 1-norm.   

            If NORMIN = 'N', CNORM is an output argument and CNORM(j)   
            returns the 1-norm of the offdiagonal part of the j-th column   
            of A.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -k, the k-th argument had an illegal value   

    Further Details   
    ======= =======   

    A rough bound on x is computed; if that is less than overflow, DTRSV   
    is called, otherwise, specific code is used which checks for possible   
    overflow or divide-by-zero at every operation.   

    A columnwise scheme is used for solving A*x = b.  The basic algorithm   
    if A is lower triangular is   

         x[1:n] := b[1:n]   
         for j = 1, ..., n   
              x(j) := x(j) / A(j,j)   
              x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]   
         end   

    Define bounds on the components of x after j iterations of the loop:   
       M(j) = bound on x[1:j]   
       G(j) = bound on x[j+1:n]   
    Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.   

    Then for iteration j+1 we have   
       M(j+1) <= G(j) / | A(j+1,j+1) |   
       G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |   
              <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )   

    where CNORM(j+1) is greater than or equal to the infinity-norm of   
    column j+1 of A, not counting the diagonal.  Hence   

       G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )   
                    1<=i<=j   
    and   

       |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )   
                                     1<=i< j   

    Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTRSV if the   
    reciprocal of the largest M(j), j=1,..,n, is larger than   
    max(underflow, 1/overflow).   

    The bound on x(j) is also used to determine when a step in the   
    columnwise method can be performed without fear of overflow.  If   
    the computed bound is greater than a large constant, x is scaled to   
    prevent overflow, but if the bound overflows, x is set to 0, x(j) to   
    1, and scale to 0, and a non-trivial solution to A*x = 0 is found.   

    Similarly, a row-wise scheme is used to solve A'*x = b.  The basic   
    algorithm for A upper triangular is   

         for j = 1, ..., n   
              x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)   
         end   

    We simultaneously compute two bounds   
         G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j   
         M(j) = bound on x(i), 1<=i<=j   

    The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we   
    add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.   
    Then the bound on x(j) is   

         M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |   

              <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )   
                        1<=i<=j   

    and we can safely call DTRSV if 1/M(n) and 1/G(n) are both greater   
    than max(underflow, 1/overflow).   

    =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     Treal c_b36 = .5;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    Treal d__1, d__2, d__3;
    /* Local variables */
     integer jinc;
     Treal xbnd;
     integer imax;
     Treal tmax, tjjs, xmax, grow, sumj;
     integer i__, j;
     Treal tscal, uscal;
     integer jlast;
     logical upper;
     Treal xj;
     Treal bignum;
     logical notran;
     integer jfirst;
     Treal smlnum;
     logical nounit;
     Treal rec, tjj;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
    --cnorm;

    /* Function Body */
    *info = 0;
    upper = template_blas_lsame(uplo, "U");
    notran = template_blas_lsame(trans, "N");
    nounit = template_blas_lsame(diag, "N");

/*     Test the input parameters. */

    if (! upper && ! template_blas_lsame(uplo, "L")) {
	*info = -1;
    } else if (! notran && ! template_blas_lsame(trans, "T") && ! 
	    template_blas_lsame(trans, "C")) {
	*info = -2;
    } else if (! nounit && ! template_blas_lsame(diag, "U")) {
	*info = -3;
    } else if (! template_blas_lsame(normin, "Y") && ! template_blas_lsame(normin,
	     "N")) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("LATRS ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine machine dependent parameters to control overflow. */

    smlnum = template_lapack_lamch("Safe minimum", (Treal)0) / template_lapack_lamch("Precision", (Treal)0);
    bignum = 1. / smlnum;
    *scale = 1.;

    if (template_blas_lsame(normin, "N")) {

/*        Compute the 1-norm of each column, not including the diagonal. */

	if (upper) {

/*           A is upper triangular. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 1;
		cnorm[j] = template_blas_asum(&i__2, &a_ref(1, j), &c__1);
/* L10: */
	    }
	} else {

/*           A is lower triangular. */

	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		cnorm[j] = template_blas_asum(&i__2, &a_ref(j + 1, j), &c__1);
/* L20: */
	    }
	    cnorm[*n] = 0.;
	}
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is   
       greater than BIGNUM. */

    imax = template_blas_idamax(n, &cnorm[1], &c__1);
    tmax = cnorm[imax];
    if (tmax <= bignum) {
	tscal = 1.;
    } else {
	tscal = 1. / (smlnum * tmax);
	dscal_(n, &tscal, &cnorm[1], &c__1);
    }

/*     Compute a bound on the computed solution vector to see if the   
       Level 2 BLAS routine DTRSV can be used. */

    j = template_blas_idamax(n, &x[1], &c__1);
    xmax = (d__1 = x[j], absMACRO(d__1));
    xbnd = xmax;
    if (notran) {

/*        Compute the growth in A * x = b. */

	if (upper) {
	    jfirst = *n;
	    jlast = 1;
	    jinc = -1;
	} else {
	    jfirst = 1;
	    jlast = *n;
	    jinc = 1;
	}

	if (tscal != 1.) {
	    grow = 0.;
	    goto L50;
	}

	if (nounit) {

/*           A is non-unit triangular.   

             Compute GROW = 1/G(j) and XBND = 1/M(j).   
             Initially, G(0) = max{x(i), i=1,...,n}. */

	    grow = 1. / maxMACRO(xbnd,smlnum);
	    xbnd = grow;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

		if (grow <= smlnum) {
		    goto L50;
		}

/*              M(j) = G(j-1) / absMACRO(A(j,j)) */

		tjj = (d__1 = a_ref(j, j), absMACRO(d__1));
/* Computing MIN */
		d__1 = xbnd, d__2 = minMACRO(1.,tjj) * grow;
		xbnd = minMACRO(d__1,d__2);
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / absMACRO(A(j,j)) ) */

		    grow *= tjj / (tjj + cnorm[j]);
		} else {

/*                 G(j) could overflow, set GROW to 0. */

		    grow = 0.;
		}
/* L30: */
	    }
	    grow = xbnd;
	} else {

/*           A is unit triangular.   

             Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.   

   Computing MIN */
	    d__1 = 1., d__2 = 1. / maxMACRO(xbnd,smlnum);
	    grow = minMACRO(d__1,d__2);
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

		if (grow <= smlnum) {
		    goto L50;
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

		grow *= 1. / (cnorm[j] + 1.);
/* L40: */
	    }
	}
L50:

	;
    } else {

/*        Compute the growth in A' * x = b. */

	if (upper) {
	    jfirst = 1;
	    jlast = *n;
	    jinc = 1;
	} else {
	    jfirst = *n;
	    jlast = 1;
	    jinc = -1;
	}

	if (tscal != 1.) {
	    grow = 0.;
	    goto L80;
	}

	if (nounit) {

/*           A is non-unit triangular.   

             Compute GROW = 1/G(j) and XBND = 1/M(j).   
             Initially, M(0) = max{x(i), i=1,...,n}. */

	    grow = 1. / maxMACRO(xbnd,smlnum);
	    xbnd = grow;
	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

		if (grow <= smlnum) {
		    goto L80;
		}

/*              G(j) = maxMACRO( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

		xj = cnorm[j] + 1.;
/* Computing MIN */
		d__1 = grow, d__2 = xbnd / xj;
		grow = minMACRO(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / absMACRO(A(j,j)) */

		tjj = (d__1 = a_ref(j, j), absMACRO(d__1));
		if (xj > tjj) {
		    xbnd *= tjj / xj;
		}
/* L60: */
	    }
	    grow = minMACRO(grow,xbnd);
	} else {

/*           A is unit triangular.   

             Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.   

   Computing MIN */
	    d__1 = 1., d__2 = 1. / maxMACRO(xbnd,smlnum);
	    grow = minMACRO(d__1,d__2);
	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

		if (grow <= smlnum) {
		    goto L80;
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

		xj = cnorm[j] + 1.;
		grow /= xj;
/* L70: */
	    }
	}
L80:
	;
    }

    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on   
          elements of X is not too small. */

	template_blas_trsv(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1);
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal to   
             BIGNUM in absolute value. */

	    *scale = bignum / xmax;
	    dscal_(n, scale, &x[1], &c__1);
	    xmax = bignum;
	}

	if (notran) {

/*           Solve A * x = b */

	    i__1 = jlast;
	    i__2 = jinc;
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

		xj = (d__1 = x[j], absMACRO(d__1));
		if (nounit) {
		    tjjs = a_ref(j, j) * tscal;
		} else {
		    tjjs = tscal;
		    if (tscal == 1.) {
			goto L100;
		    }
		}
		tjj = absMACRO(tjjs);
		if (tjj > smlnum) {

/*                    absMACRO(A(j,j)) > SMLNUM: */

		    if (tjj < 1.) {
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

			    rec = 1. / xj;
			    dscal_(n, &rec, &x[1], &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    x[j] /= tjjs;
		    xj = (d__1 = x[j], absMACRO(d__1));
		} else if (tjj > 0.) {

/*                    0 < absMACRO(A(j,j)) <= SMLNUM: */

		    if (xj > tjj * bignum) {

/*                       Scale x by (1/absMACRO(x(j)))*absMACRO(A(j,j))*BIGNUM   
                         to avoid overflow when dividing by A(j,j). */

			rec = tjj * bignum / xj;
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when   
                            multiplying x(j) times column j. */

			    rec /= cnorm[j];
			}
			dscal_(n, &rec, &x[1], &c__1);
			*scale *= rec;
			xmax *= rec;
		    }
		    x[j] /= tjjs;
		    xj = (d__1 = x[j], absMACRO(d__1));
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and   
                      scale = 0, and compute a solution to A*x = 0. */

		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			x[i__] = 0.;
/* L90: */
		    }
		    x[j] = 1.;
		    xj = 1.;
		    *scale = 0.;
		    xmax = 0.;
		}
L100:

/*              Scale x if necessary to avoid overflow when adding a   
                multiple of column j of A. */

		if (xj > 1.) {
		    rec = 1. / xj;
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*absMACRO(x(j))). */

			rec *= .5;
			dscal_(n, &rec, &x[1], &c__1);
			*scale *= rec;
		    }
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

		    dscal_(n, &c_b36, &x[1], &c__1);
		    *scale *= .5;
		}

		if (upper) {
		    if (j > 1) {

/*                    Compute the update   
                         x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

			i__3 = j - 1;
			d__1 = -x[j] * tscal;
			daxpy_(&i__3, &d__1, &a_ref(1, j), &c__1, &x[1], &
				c__1);
			i__3 = j - 1;
			i__ = template_blas_idamax(&i__3, &x[1], &c__1);
			xmax = (d__1 = x[i__], absMACRO(d__1));
		    }
		} else {
		    if (j < *n) {

/*                    Compute the update   
                         x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

			i__3 = *n - j;
			d__1 = -x[j] * tscal;
			daxpy_(&i__3, &d__1, &a_ref(j + 1, j), &c__1, &x[j + 
				1], &c__1);
			i__3 = *n - j;
			i__ = j + template_blas_idamax(&i__3, &x[j + 1], &c__1);
			xmax = (d__1 = x[i__], absMACRO(d__1));
		    }
		}
/* L110: */
	    }

	} else {

/*           Solve A' * x = b */

	    i__2 = jlast;
	    i__1 = jinc;
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k).   
                                      k<>j */

		xj = (d__1 = x[j], absMACRO(d__1));
		uscal = tscal;
		rec = 1. / maxMACRO(xmax,1.);
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

		    rec *= .5;
		    if (nounit) {
			tjjs = a_ref(j, j) * tscal;
		    } else {
			tjjs = tscal;
		    }
		    tjj = absMACRO(tjjs);
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1.   

   Computing MIN */
			d__1 = 1., d__2 = rec * tjj;
			rec = minMACRO(d__1,d__2);
			uscal /= tjjs;
		    }
		    if (rec < 1.) {
			dscal_(n, &rec, &x[1], &c__1);
			*scale *= rec;
			xmax *= rec;
		    }
		}

		sumj = 0.;
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot product is 1,   
                   call DDOT to perform the dot product. */

		    if (upper) {
			i__3 = j - 1;
			sumj = ddot_(&i__3, &a_ref(1, j), &c__1, &x[1], &c__1)
				;
		    } else if (j < *n) {
			i__3 = *n - j;
			sumj = ddot_(&i__3, &a_ref(j + 1, j), &c__1, &x[j + 1]
				, &c__1);
		    }
		} else {

/*                 Otherwise, use in-line code for the dot product. */

		    if (upper) {
			i__3 = j - 1;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    sumj += a_ref(i__, j) * uscal * x[i__];
/* L120: */
			}
		    } else if (j < *n) {
			i__3 = *n;
			for (i__ = j + 1; i__ <= i__3; ++i__) {
			    sumj += a_ref(i__, j) * uscal * x[i__];
/* L130: */
			}
		    }
		}

		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)   
                   was not used to scale the dotproduct. */

		    x[j] -= sumj;
		    xj = (d__1 = x[j], absMACRO(d__1));
		    if (nounit) {
			tjjs = a_ref(j, j) * tscal;
		    } else {
			tjjs = tscal;
			if (tscal == 1.) {
			    goto L150;
			}
		    }

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

		    tjj = absMACRO(tjjs);
		    if (tjj > smlnum) {

/*                       absMACRO(A(j,j)) > SMLNUM: */

			if (tjj < 1.) {
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/absMACRO(x(j)). */

				rec = 1. / xj;
				dscal_(n, &rec, &x[1], &c__1);
				*scale *= rec;
				xmax *= rec;
			    }
			}
			x[j] /= tjjs;
		    } else if (tjj > 0.) {

/*                       0 < absMACRO(A(j,j)) <= SMLNUM: */

			if (xj > tjj * bignum) {

/*                          Scale x by (1/absMACRO(x(j)))*absMACRO(A(j,j))*BIGNUM. */

			    rec = tjj * bignum / xj;
			    dscal_(n, &rec, &x[1], &c__1);
			    *scale *= rec;
			    xmax *= rec;
			}
			x[j] /= tjjs;
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and   
                         scale = 0, and compute a solution to A'*x = 0. */

			i__3 = *n;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    x[i__] = 0.;
/* L140: */
			}
			x[j] = 1.;
			*scale = 0.;
			xmax = 0.;
		    }
L150:
		    ;
		} else {

/*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot   
                   product has already been divided by 1/A(j,j). */

		    x[j] = x[j] / tjjs - sumj;
		}
/* Computing MAX */
		d__2 = xmax, d__3 = (d__1 = x[j], absMACRO(d__1));
		xmax = maxMACRO(d__2,d__3);
/* L160: */
	    }
	}
	*scale /= tscal;
    }

/*     Scale the column norms by 1/TSCAL for return. */

    if (tscal != 1.) {
	d__1 = 1. / tscal;
	dscal_(n, &d__1, &cnorm[1], &c__1);
    }

    return 0;

/*     End of DLATRS */

} /* dlatrs_ */

#undef a_ref


#endif
