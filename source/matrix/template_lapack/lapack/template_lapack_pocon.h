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
 

#ifndef TEMPLATE_LAPACK_POCON_HEADER
#define TEMPLATE_LAPACK_POCON_HEADER


template<class Treal>
int template_lapack_pocon(const char *uplo, const integer *n, const Treal *a, const integer *
	lda, const Treal *anorm, Treal *rcond, Treal *work, integer *
	iwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPOCON estimates the reciprocal of the condition number (in the   
    1-norm) of a real symmetric positive definite matrix using the   
    Cholesky factorization A = U**T*U or A = L*L**T computed by DPOTRF.   

    An estimate is obtained for norm(inv(A)), and the reciprocal of the   
    condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The triangular factor U or L from the Cholesky factorization   
            A = U**T*U or A = L*L**T, as computed by DPOTRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    ANORM   (input) DOUBLE PRECISION   
            The 1-norm (or infinity-norm) of the symmetric matrix A.   

    RCOND   (output) DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an   
            estimate of the 1-norm of inv(A) computed in this routine.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    Treal d__1;
    /* Local variables */
     integer kase;
     Treal scale;
     logical upper;
     integer ix;
     Treal scalel;
     Treal scaleu;
     Treal ainvnm;
     char normin[1];
     Treal smlnum;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --work;
    --iwork;

    /* Function Body */
    *info = 0;
    upper = template_blas_lsame(uplo, "U");
    if (! upper && ! template_blas_lsame(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -4;
    } else if (*anorm < 0.) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("POCON ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return 0;
    } else if (*anorm == 0.) {
	return 0;
    }

    smlnum = template_lapack_lamch("Safe minimum", (Treal)0);

/*     Estimate the 1-norm of inv(A). */

    kase = 0;
    *(unsigned char *)normin = 'N';
L10:
    template_lapack_lacon(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase);
    if (kase != 0) {
	if (upper) {

/*           Multiply by inv(U'). */

	    template_lapack_latrs("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset],
		     lda, &work[1], &scalel, &work[(*n << 1) + 1], info);
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */

	    template_lapack_latrs("Upper", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &scaleu, &work[(*n << 1) + 1], 
		    info);
	} else {

/*           Multiply by inv(L). */

	    template_lapack_latrs("Lower", "No transpose", "Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &scalel, &work[(*n << 1) + 1], 
		    info);
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L'). */

	    template_lapack_latrs("Lower", "Transpose", "Non-unit", normin, n, &a[a_offset],
		     lda, &work[1], &scaleu, &work[(*n << 1) + 1], info);
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

	scale = scalel * scaleu;
	if (scale != 1.) {
	    ix = template_blas_idamax(n, &work[1], &c__1);
	    if (scale < (d__1 = work[ix], absMACRO(d__1)) * smlnum || scale == 0.) 
		    {
		goto L20;
	    }
	    template_lapack_rscl(n, &scale, &work[1], &c__1);
	}
	goto L10;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

L20:
    return 0;

/*     End of DPOCON */

} /* dpocon_ */

#endif
