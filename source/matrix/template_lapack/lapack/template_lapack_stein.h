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
 

#ifndef TEMPLATE_LAPACK_STEIN_HEADER
#define TEMPLATE_LAPACK_STEIN_HEADER


template<class Treal>
int template_lapack_stein(const integer *n, const Treal *d__, const Treal *e, 
	const integer *m, const Treal *w, const integer *iblock, const integer *isplit, 
	Treal *z__, const integer *ldz, Treal *work, integer *iwork, 
	integer *ifail, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEIN computes the eigenvectors of a real symmetric tridiagonal   
    matrix T corresponding to specified eigenvalues, using inverse   
    iteration.   

    The maximum number of iterations allowed for each eigenvector is   
    specified by an internal parameter MAXITS (currently set to 5).   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the tridiagonal matrix T.   

    E       (input) DOUBLE PRECISION array, dimension (N)   
            The (n-1) subdiagonal elements of the tridiagonal matrix   
            T, in elements 1 to N-1.  E(N) need not be set.   

    M       (input) INTEGER   
            The number of eigenvectors to be found.  0 <= M <= N.   

    W       (input) DOUBLE PRECISION array, dimension (N)   
            The first M elements of W contain the eigenvalues for   
            which eigenvectors are to be computed.  The eigenvalues   
            should be grouped by split-off block and ordered from   
            smallest to largest within the block.  ( The output array   
            W from DSTEBZ with ORDER = 'B' is expected here. )   

    IBLOCK  (input) INTEGER array, dimension (N)   
            The submatrix indices associated with the corresponding   
            eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to   
            the first submatrix from the top, =2 if W(i) belongs to   
            the second submatrix, etc.  ( The output array IBLOCK   
            from DSTEBZ is expected here. )   

    ISPLIT  (input) INTEGER array, dimension (N)   
            The splitting points, at which T breaks up into submatrices.   
            The first submatrix consists of rows/columns 1 to   
            ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1   
            through ISPLIT( 2 ), etc.   
            ( The output array ISPLIT from DSTEBZ is expected here. )   

    Z       (output) DOUBLE PRECISION array, dimension (LDZ, M)   
            The computed eigenvectors.  The eigenvector associated   
            with the eigenvalue W(i) is stored in the i-th column of   
            Z.  Any vector which fails to converge is set to its current   
            iterate after MAXITS iterations.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= max(1,N).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    IFAIL   (output) INTEGER array, dimension (M)   
            On normal exit, all elements of IFAIL are zero.   
            If one or more eigenvectors fail to converge after   
            MAXITS iterations, then their indices are stored in   
            array IFAIL.   

    INFO    (output) INTEGER   
            = 0: successful exit.   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, then i eigenvectors failed to converge   
                 in MAXITS iterations.  Their indices are stored in   
                 array IFAIL.   

    Internal Parameters   
    ===================   

    MAXITS  INTEGER, default = 5   
            The maximum number of iterations performed.   

    EXTRA   INTEGER, default = 2   
            The number of iterations performed after norm growth   
            criterion is satisfied, should be at least 1.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__2 = 2;
     integer c__1 = 1;
     integer c_n1 = -1;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    Treal d__1, d__2, d__3, d__4, d__5;
    /* Local variables */
     integer jblk, nblk;
     integer jmax;
     integer i__, j;
     integer iseed[4], gpind, iinfo;
     integer b1;
     integer j1;
     Treal ortol;
     integer indrv1, indrv2, indrv3, indrv4, indrv5, bn;
     Treal xj;
     integer nrmchk;
     integer blksiz;
     Treal onenrm, dtpcrt, pertol, scl, eps, sep, nrm, tol;
     integer its;
     Treal xjm, ztr, eps1;
#define z___ref(a_1,a_2) z__[(a_2)*z_dim1 + a_1]


    --d__;
    --e;
    --w;
    --iblock;
    --isplit;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1 * 1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;

    /* Initialization added by Elias to get rid of compiler warnings. */
    ortol = dtpcrt = xjm = onenrm = gpind = 0;
    /* Function Body */
    *info = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ifail[i__] = 0;
/* L10: */
    }

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0 || *m > *n) {
	*info = -4;
    } else if (*ldz < maxMACRO(1,*n)) {
	*info = -9;
    } else {
	i__1 = *m;
	for (j = 2; j <= i__1; ++j) {
	    if (iblock[j] < iblock[j - 1]) {
		*info = -6;
		goto L30;
	    }
	    if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
		*info = -5;
		goto L30;
	    }
/* L20: */
	}
L30:
	;
    }

    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("STEIN ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    } else if (*n == 1) {
	z___ref(1, 1) = 1.;
	return 0;
    }

/*     Get machine constants. */

    eps = template_lapack_lamch("Precision", (Treal)0);

/*     Initialize seed for random number generator DLARNV. */

    for (i__ = 1; i__ <= 4; ++i__) {
	iseed[i__ - 1] = 1;
/* L40: */
    }

/*     Initialize pointers. */

    indrv1 = 0;
    indrv2 = indrv1 + *n;
    indrv3 = indrv2 + *n;
    indrv4 = indrv3 + *n;
    indrv5 = indrv4 + *n;

/*     Compute eigenvectors of matrix blocks. */

    j1 = 1;
    i__1 = iblock[*m];
    for (nblk = 1; nblk <= i__1; ++nblk) {

/*        Find starting and ending indices of block nblk. */

	if (nblk == 1) {
	    b1 = 1;
	} else {
	    b1 = isplit[nblk - 1] + 1;
	}
	bn = isplit[nblk];
	blksiz = bn - b1 + 1;
	if (blksiz == 1) {
	    goto L60;
	}
	gpind = b1;

/*        Compute reorthogonalization criterion and stopping criterion. */

	onenrm = (d__1 = d__[b1], absMACRO(d__1)) + (d__2 = e[b1], absMACRO(d__2));
/* Computing MAX */
	d__3 = onenrm, d__4 = (d__1 = d__[bn], absMACRO(d__1)) + (d__2 = e[bn - 1],
		 absMACRO(d__2));
	onenrm = maxMACRO(d__3,d__4);
	i__2 = bn - 1;
	for (i__ = b1 + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__4 = onenrm, d__5 = (d__1 = d__[i__], absMACRO(d__1)) + (d__2 = e[
		    i__ - 1], absMACRO(d__2)) + (d__3 = e[i__], absMACRO(d__3));
	    onenrm = maxMACRO(d__4,d__5);
/* L50: */
	}
	ortol = onenrm * .001;

	dtpcrt = template_blas_sqrt(.1 / blksiz);

/*        Loop through eigenvalues of block nblk. */

L60:
	jblk = 0;
	i__2 = *m;
	for (j = j1; j <= i__2; ++j) {
	    if (iblock[j] != nblk) {
		j1 = j;
		goto L160;
	    }
	    ++jblk;
	    xj = w[j];

/*           Skip all the work if the block size is one. */

	    if (blksiz == 1) {
		work[indrv1 + 1] = 1.;
		goto L120;
	    }

/*           If eigenvalues j and j-1 are too close, add a relatively   
             small perturbation. */

	    if (jblk > 1) {
		eps1 = (d__1 = eps * xj, absMACRO(d__1));
		pertol = eps1 * 10.;
		sep = xj - xjm;
		if (sep < pertol) {
		    xj = xjm + pertol;
		}
	    }

	    its = 0;
	    nrmchk = 0;

/*           Get random starting vector. */

	    template_lapack_larnv(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

/*           Copy the matrix T so it won't be destroyed in factorization. */

	    template_blas_copy(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
	    i__3 = blksiz - 1;
	    template_blas_copy(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
	    i__3 = blksiz - 1;
	    template_blas_copy(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

/*           Compute LU factors with partial pivoting  ( PT = LU ) */

	    tol = 0.;
	    template_lapack_lagtf(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
		    indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

/*           Update iteration count. */

L70:
	    ++its;
	    if (its > 5) {
		goto L100;
	    }

/*           Normalize and scale the righthand side vector Pb.   

   Computing MAX */
	    d__2 = eps, d__3 = (d__1 = work[indrv4 + blksiz], absMACRO(d__1));
	    scl = blksiz * onenrm * maxMACRO(d__2,d__3) / template_blas_asum(&blksiz, &work[
		    indrv1 + 1], &c__1);
	    template_blas_scal(&blksiz, &scl, &work[indrv1 + 1], &c__1);

/*           Solve the system LU = Pb. */

	    template_lapack_lagts(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
		    work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
		    indrv1 + 1], &tol, &iinfo);

/*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are   
             close enough. */

	    if (jblk == 1) {
		goto L90;
	    }
	    if ((d__1 = xj - xjm, absMACRO(d__1)) > ortol) {
		gpind = j;
	    }
	    if (gpind != j) {
		i__3 = j - 1;
		for (i__ = gpind; i__ <= i__3; ++i__) {
		    ztr = -template_blas_dot(&blksiz, &work[indrv1 + 1], &c__1, &z___ref(
			    b1, i__), &c__1);
		    template_blas_axpy(&blksiz, &ztr, &z___ref(b1, i__), &c__1, &work[
			    indrv1 + 1], &c__1);
/* L80: */
		}
	    }

/*           Check the infinity norm of the iterate. */

L90:
	    jmax = template_blas_idamax(&blksiz, &work[indrv1 + 1], &c__1);
	    nrm = (d__1 = work[indrv1 + jmax], absMACRO(d__1));

/*           Continue for additional iterations after norm reaches   
             stopping criterion. */

	    if (nrm < dtpcrt) {
		goto L70;
	    }
	    ++nrmchk;
	    if (nrmchk < 3) {
		goto L70;
	    }

	    goto L110;

/*           If stopping criterion was not satisfied, update info and   
             store eigenvector number in array ifail. */

L100:
	    ++(*info);
	    ifail[*info] = j;

/*           Accept iterate as jth eigenvector. */

L110:
	    scl = 1. / template_blas_nrm2(&blksiz, &work[indrv1 + 1], &c__1);
	    jmax = template_blas_idamax(&blksiz, &work[indrv1 + 1], &c__1);
	    if (work[indrv1 + jmax] < 0.) {
		scl = -scl;
	    }
	    template_blas_scal(&blksiz, &scl, &work[indrv1 + 1], &c__1);
L120:
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		z___ref(i__, j) = 0.;
/* L130: */
	    }
	    i__3 = blksiz;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		z___ref(b1 + i__ - 1, j) = work[indrv1 + i__];
/* L140: */
	    }

/*           Save the shift to check eigenvalue spacing at next   
             iteration. */

	    xjm = xj;

/* L150: */
	}
L160:
	;
    }

    return 0;

/*     End of DSTEIN */

} /* dstein_ */

#undef z___ref


#endif
