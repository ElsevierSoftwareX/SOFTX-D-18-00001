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
 

#ifndef TEMPLATE_LAPACK_STEVX_HEADER
#define TEMPLATE_LAPACK_STEVX_HEADER


template<class Treal>
int template_lapack_stevx(const char *jobz, const char *range, const integer *n, Treal *
	d__, Treal *e, const Treal *vl, const Treal *vu, const integer *il, 
	const integer *iu, const Treal *abstol, integer *m, Treal *w, 
	Treal *z__, const integer *ldz, Treal *work, integer *iwork, 
	integer *ifail, integer *info)
{
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DSTEVX computes selected eigenvalues and, optionally, eigenvectors   
    of a real symmetric tridiagonal matrix A.  Eigenvalues and   
    eigenvectors can be selected by specifying either a range of values   
    or a range of indices for the desired eigenvalues.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    RANGE   (input) CHARACTER*1   
            = 'A': all eigenvalues will be found.   
            = 'V': all eigenvalues in the half-open interval (VL,VU]   
                   will be found.   
            = 'I': the IL-th through IU-th eigenvalues will be found.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.   
            On exit, D may be multiplied by a constant factor chosen   
            to avoid over/underflow in computing the eigenvalues.   

    E       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix A in elements 1 to N-1 of E; E(N) need not be set.   
            On exit, E may be multiplied by a constant factor chosen   
            to avoid over/underflow in computing the eigenvalues.   

    VL      (input) DOUBLE PRECISION   
    VU      (input) DOUBLE PRECISION   
            If RANGE='V', the lower and upper bounds of the interval to   
            be searched for eigenvalues. VL < VU.   
            Not referenced if RANGE = 'A' or 'I'.   

    IL      (input) INTEGER   
    IU      (input) INTEGER   
            If RANGE='I', the indices (in ascending order) of the   
            smallest and largest eigenvalues to be returned.   
            1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.   
            Not referenced if RANGE = 'A' or 'V'.   

    ABSTOL  (input) DOUBLE PRECISION   
            The absolute error tolerance for the eigenvalues.   
            An approximate eigenvalue is accepted as converged   
            when it is determined to lie in an interval [a,b]   
            of width less than or equal to   

                    ABSTOL + EPS *   max( |a|,|b| ) ,   

            where EPS is the machine precision.  If ABSTOL is less   
            than or equal to zero, then  EPS*|T|  will be used in   
            its place, where |T| is the 1-norm of the tridiagonal   
            matrix.   

            Eigenvalues will be computed most accurately when ABSTOL is   
            set to twice the underflow threshold 2*DLAMCH('S'), not zero.   
            If this routine returns with INFO>0, indicating that some   
            eigenvectors did not converge, try setting ABSTOL to   
            2*DLAMCH('S').   

            See "Computing Small Singular Values of Bidiagonal Matrices   
            with Guaranteed High Relative Accuracy," by Demmel and   
            Kahan, LAPACK Working Note #3.   

    M       (output) INTEGER   
            The total number of eigenvalues found.  0 <= M <= N.   
            If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            The first M elements contain the selected eigenvalues in   
            ascending order.   

    Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )   
            If JOBZ = 'V', then if INFO = 0, the first M columns of Z   
            contain the orthonormal eigenvectors of the matrix A   
            corresponding to the selected eigenvalues, with the i-th   
            column of Z holding the eigenvector associated with W(i).   
            If an eigenvector fails to converge (INFO > 0), then that   
            column of Z contains the latest approximation to the   
            eigenvector, and the index of the eigenvector is returned   
            in IFAIL.  If JOBZ = 'N', then Z is not referenced.   
            Note: the user must ensure that at least max(1,M) columns are   
            supplied in the array Z; if RANGE = 'V', the exact value of M   
            is not known in advance and an upper bound must be used.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= max(1,N).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)   

    IWORK   (workspace) INTEGER array, dimension (5*N)   

    IFAIL   (output) INTEGER array, dimension (N)   
            If JOBZ = 'V', then if INFO = 0, the first M elements of   
            IFAIL are zero.  If INFO > 0, then IFAIL contains the   
            indices of the eigenvectors that failed to converge.   
            If JOBZ = 'N', then IFAIL is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, then i eigenvectors failed to converge.   
                  Their indices are stored in array IFAIL.   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    Treal d__1, d__2;
    /* Local variables */
     integer imax;
     Treal rmin, rmax, tnrm;
     integer itmp1, i__, j;
     Treal sigma;
     char order[1];
     logical wantz;
     integer jj;
     logical alleig, indeig;
     integer iscale, indibl;
     logical valeig;
     Treal safmin;
     Treal bignum;
     integer indisp;
     integer indiwo;
     integer indwrk;
     integer nsplit;
     Treal smlnum, eps, vll, vuu, tmp1;
#define z___ref(a_1,a_2) z__[(a_2)*z_dim1 + a_1]


    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1 * 1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;

    /* Function Body */
    wantz = template_blas_lsame(jobz, "V");
    alleig = template_blas_lsame(range, "A");
    valeig = template_blas_lsame(range, "V");
    indeig = template_blas_lsame(range, "I");

    *info = 0;
    if (! (wantz || template_blas_lsame(jobz, "N"))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else {
	if (valeig) {
	    if (*n > 0 && *vu <= *vl) {
		*info = -7;
	    }
	} else if (indeig) {
	    if (*il < 1 || *il > maxMACRO(1,*n)) {
		*info = -8;
	    } else if (*iu < minMACRO(*n,*il) || *iu > *n) {
		*info = -9;
	    }
	}
    }
    if (*info == 0) {
      if (*ldz < 1 || (wantz && *ldz < *n) ) {
	    *info = -14;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("STEVX ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *m = 0;
    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = d__[1];
	} else {
	    if (*vl < d__[1] && *vu >= d__[1]) {
		*m = 1;
		w[1] = d__[1];
	    }
	}
	if (wantz) {
	    z___ref(1, 1) = 1.;
	}
	return 0;
    }

/*     Get machine constants. */

    safmin = template_lapack_lamch("Safe minimum", (Treal)0);
    eps = template_lapack_lamch("Precision", (Treal)0);
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = template_blas_sqrt(smlnum);
/* Computing MIN */
    d__1 = template_blas_sqrt(bignum), d__2 = 1. / template_blas_sqrt(template_blas_sqrt(safmin));
    rmax = minMACRO(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

    iscale = 0;
    if (valeig) {
	vll = *vl;
	vuu = *vu;
    } else {
	vll = 0.;
	vuu = 0.;
    }
    tnrm = template_lapack_lanst("M", n, &d__[1], &e[1]);
    if (tnrm > 0. && tnrm < rmin) {
	iscale = 1;
	sigma = rmin / tnrm;
    } else if (tnrm > rmax) {
	iscale = 1;
	sigma = rmax / tnrm;
    }
    if (iscale == 1) {
	template_blas_scal(n, &sigma, &d__[1], &c__1);
	i__1 = *n - 1;
	template_blas_scal(&i__1, &sigma, &e[1], &c__1);
	if (valeig) {
	    vll = *vl * sigma;
	    vuu = *vu * sigma;
	}
    }

/*     If all eigenvalues are desired and ABSTOL is less than zero, then   
       call DSTERF or SSTEQR.  If this fails for some eigenvalue, then   
       try DSTEBZ. */

    if ((alleig || (indeig && *il == 1 && *iu == *n) ) && *abstol <= 0.) {
	template_blas_copy(n, &d__[1], &c__1, &w[1], &c__1);
	i__1 = *n - 1;
	template_blas_copy(&i__1, &e[1], &c__1, &work[1], &c__1);
	indwrk = *n + 1;
	if (! wantz) {
	    template_lapack_sterf(n, &w[1], &work[1], info);
	} else {
	    template_lapack_steqr("I", n, &w[1], &work[1], &z__[z_offset], ldz, &work[
		    indwrk], info);
	    if (*info == 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ifail[i__] = 0;
/* L10: */
		}
	    }
	}
	if (*info == 0) {
	    *m = *n;
	    goto L20;
	}
	*info = 0;
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */

    if (wantz) {
	*(unsigned char *)order = 'B';
    } else {
	*(unsigned char *)order = 'E';
    }
    indwrk = 1;
    indibl = 1;
    indisp = indibl + *n;
    indiwo = indisp + *n;
    template_lapack_stebz(range, order, n, &vll, &vuu, il, iu, abstol, &d__[1], &e[1], m, &
	    nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[indwrk], &
	    iwork[indiwo], info);

    if (wantz) {
	template_lapack_stein(n, &d__[1], &e[1], m, &w[1], &iwork[indibl], &iwork[indisp], &
		z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &ifail[1], 
		info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

L20:
    if (iscale == 1) {
	if (*info == 0) {
	    imax = *m;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	template_blas_scal(&imax, &d__1, &w[1], &c__1);
    }

/*     If eigenvalues are not in order, then sort them, along with   
       eigenvectors. */

    if (wantz) {
	i__1 = *m - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp1 = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp1) {
		    i__ = jj;
		    tmp1 = w[jj];
		}
/* L30: */
	    }

	    if (i__ != 0) {
		itmp1 = iwork[indibl + i__ - 1];
		w[i__] = w[j];
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
		w[j] = tmp1;
		iwork[indibl + j - 1] = itmp1;
		template_blas_swap(n, &z___ref(1, i__), &c__1, &z___ref(1, j), &c__1);
		if (*info != 0) {
		    itmp1 = ifail[i__];
		    ifail[i__] = ifail[j];
		    ifail[j] = itmp1;
		}
	    }
/* L40: */
	}
    }

    return 0;

/*     End of DSTEVX */

} /* dstevx_ */

#undef z___ref


#endif
