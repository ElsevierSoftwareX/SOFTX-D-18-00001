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
 

#ifndef TEMPLATE_LAPACK_GGBAK_HEADER
#define TEMPLATE_LAPACK_GGBAK_HEADER


template<class Treal>
int template_lapack_ggbak(const char *job, const char *side, const integer *n, const integer *ilo, 
	const integer *ihi, const Treal *lscale, const Treal *rscale, const integer *m, 
	Treal *v, const integer *ldv, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGBAK forms the right or left eigenvectors of a real generalized   
    eigenvalue problem A*x = lambda*B*x, by backward transformation on   
    the computed eigenvectors of the balanced pair of matrices output by   
    DGGBAL.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the type of backward transformation required:   
            = 'N':  do nothing, return immediately;   
            = 'P':  do backward transformation for permutation only;   
            = 'S':  do backward transformation for scaling only;   
            = 'B':  do backward transformations for both permutation and   
                    scaling.   
            JOB must be the same as the argument JOB supplied to DGGBAL.   

    SIDE    (input) CHARACTER*1   
            = 'R':  V contains right eigenvectors;   
            = 'L':  V contains left eigenvectors.   

    N       (input) INTEGER   
            The number of rows of the matrix V.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            The integers ILO and IHI determined by DGGBAL.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    LSCALE  (input) DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and/or scaling factors applied   
            to the left side of A and B, as returned by DGGBAL.   

    RSCALE  (input) DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and/or scaling factors applied   
            to the right side of A and B, as returned by DGGBAL.   

    M       (input) INTEGER   
            The number of columns of the matrix V.  M >= 0.   

    V       (input/output) DOUBLE PRECISION array, dimension (LDV,M)   
            On entry, the matrix of right or left eigenvectors to be   
            transformed, as returned by DTGEVC.   
            On exit, V is overwritten by the transformed eigenvectors.   

    LDV     (input) INTEGER   
            The leading dimension of the matrix V. LDV >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    See R.C. Ward, Balancing the generalized eigenvalue problem,   
                   SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.   

    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    /* System generated locals */
    integer v_dim1, v_offset, i__1;
    /* Local variables */
     integer i__, k;
     logical leftv;
     logical rightv;
#define v_ref(a_1,a_2) v[(a_2)*v_dim1 + a_1]

    --lscale;
    --rscale;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1 * 1;
    v -= v_offset;

    /* Function Body */
    rightv = template_blas_lsame(side, "R");
    leftv = template_blas_lsame(side, "L");

    *info = 0;
    if (! template_blas_lsame(job, "N") && ! template_blas_lsame(job, "P") && ! template_blas_lsame(job, "S") 
	    && ! template_blas_lsame(job, "B")) {
	*info = -1;
    } else if (! rightv && ! leftv) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1) {
	*info = -4;
    } else if (*ihi < *ilo || *ihi > maxMACRO(1,*n)) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*ldv < maxMACRO(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("GGBAK ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }
    if (*m == 0) {
	return 0;
    }
    if (template_blas_lsame(job, "N")) {
	return 0;
    }

    if (*ilo == *ihi) {
	goto L30;
    }

/*     Backward balance */

    if (template_blas_lsame(job, "S") || template_blas_lsame(job, "B")) {

/*        Backward transformation on right eigenvectors */

	if (rightv) {
	    i__1 = *ihi;
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
		template_blas_scal(m, &rscale[i__], &v_ref(i__, 1), ldv);
/* L10: */
	    }
	}

/*        Backward transformation on left eigenvectors */

	if (leftv) {
	    i__1 = *ihi;
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
		template_blas_scal(m, &lscale[i__], &v_ref(i__, 1), ldv);
/* L20: */
	    }
	}
    }

/*     Backward permutation */

L30:
    if (template_blas_lsame(job, "P") || template_blas_lsame(job, "B")) {

/*        Backward permutation on right eigenvectors */

	if (rightv) {
	    if (*ilo == 1) {
		goto L50;
	    }

	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
		k = (integer) rscale[i__];
		if (k == i__) {
		    goto L40;
		}
		template_blas_swap(m, &v_ref(i__, 1), ldv, &v_ref(k, 1), ldv);
L40:
		;
	    }

L50:
	    if (*ihi == *n) {
		goto L70;
	    }
	    i__1 = *n;
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
		k = (integer) rscale[i__];
		if (k == i__) {
		    goto L60;
		}
		template_blas_swap(m, &v_ref(i__, 1), ldv, &v_ref(k, 1), ldv);
L60:
		;
	    }
	}

/*        Backward permutation on left eigenvectors */

L70:
	if (leftv) {
	    if (*ilo == 1) {
		goto L90;
	    }
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
		k = (integer) lscale[i__];
		if (k == i__) {
		    goto L80;
		}
		template_blas_swap(m, &v_ref(i__, 1), ldv, &v_ref(k, 1), ldv);
L80:
		;
	    }

L90:
	    if (*ihi == *n) {
		goto L110;
	    }
	    i__1 = *n;
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
		k = (integer) lscale[i__];
		if (k == i__) {
		    goto L100;
		}
		template_blas_swap(m, &v_ref(i__, 1), ldv, &v_ref(k, 1), ldv);
L100:
		;
	    }
	}
    }

L110:

    return 0;

/*     End of DGGBAK */

} /* dggbak_ */

#undef v_ref


#endif
