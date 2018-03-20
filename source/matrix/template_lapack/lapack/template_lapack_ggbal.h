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
 

#ifndef TEMPLATE_LAPACK_GGBAL_HEADER
#define TEMPLATE_LAPACK_GGBAL_HEADER


template<class Treal>
int template_lapack_ggbal(const char *job, const integer *n, Treal *a, const integer *
	lda, Treal *b, const integer *ldb, integer *ilo, integer *ihi, 
	Treal *lscale, Treal *rscale, Treal *work, integer *
	info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGBAL balances a pair of general real matrices (A,B).  This   
    involves, first, permuting A and B by similarity transformations to   
    isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N   
    elements on the diagonal; and second, applying a diagonal similarity   
    transformation to rows and columns ILO to IHI to make the rows   
    and columns as close in norm as possible. Both steps are optional.   

    Balancing may reduce the 1-norm of the matrices, and improve the   
    accuracy of the computed eigenvalues and/or eigenvectors in the   
    generalized eigenvalue problem A*x = lambda*B*x.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the operations to be performed on A and B:   
            = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0   
                    and RSCALE(I) = 1.0 for i = 1,...,N.   
            = 'P':  permute only;   
            = 'S':  scale only;   
            = 'B':  both permute and scale.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the input matrix A.   
            On exit,  A is overwritten by the balanced matrix.   
            If JOB = 'N', A is not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,N).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the input matrix B.   
            On exit,  B is overwritten by the balanced matrix.   
            If JOB = 'N', B is not referenced.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,N).   

    ILO     (output) INTEGER   
    IHI     (output) INTEGER   
            ILO and IHI are set to integers such that on exit   
            A(i,j) = 0 and B(i,j) = 0 if i > j and   
            j = 1,...,ILO-1 or i = IHI+1,...,N.   
            If JOB = 'N' or 'S', ILO = 1 and IHI = N.   

    LSCALE  (output) DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and scaling factors applied   
            to the left side of A and B.  If P(j) is the index of the   
            row interchanged with row j, and D(j)   
            is the scaling factor applied to row j, then   
              LSCALE(j) = P(j)    for J = 1,...,ILO-1   
                        = D(j)    for J = ILO,...,IHI   
                        = P(j)    for J = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    RSCALE  (output) DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and scaling factors applied   
            to the right side of A and B.  If P(j) is the index of the   
            column interchanged with column j, and D(j)   
            is the scaling factor applied to column j, then   
              LSCALE(j) = P(j)    for J = 1,...,ILO-1   
                        = D(j)    for J = ILO,...,IHI   
                        = P(j)    for J = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    See R.C. WARD, Balancing the generalized eigenvalue problem,   
                   SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.   

    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
     Treal c_b34 = 10.;
     Treal c_b70 = .5;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    Treal d__1, d__2, d__3;
    /* Local variables */
     integer lcab;
     Treal beta, coef;
     integer irab, lrab;
     Treal basl, cmax;
     Treal coef2, coef5;
     integer i__, j, k, l, m;
     Treal gamma, t, alpha;
     Treal sfmin, sfmax;
     integer iflow;
     integer kount, jc;
     Treal ta, tb, tc;
     integer ir, it;
     Treal ew;
     integer nr;
     Treal pgamma;
     integer lsfmin, lsfmax, ip1, jp1, lm1;
     Treal cab, rab, ewc, cor, sum;
     integer nrp2, icab;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    --lscale;
    --rscale;
    --work;

    /* Initialization added by Elias to get rid of compiler warnings. */
    pgamma = 0;
    /* Function Body */
    *info = 0;
    if (! template_blas_lsame(job, "N") && ! template_blas_lsame(job, "P") && ! template_blas_lsame(job, "S") 
	    && ! template_blas_lsame(job, "B")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < maxMACRO(1,*n)) {
	*info = -4;
    } else if (*ldb < maxMACRO(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	template_blas_erbla("GGBAL ", &i__1);
	return 0;
    }

    k = 1;
    l = *n;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (template_blas_lsame(job, "N")) {
	*ilo = 1;
	*ihi = *n;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lscale[i__] = 1.;
	    rscale[i__] = 1.;
/* L10: */
	}
	return 0;
    }

    if (k == l) {
	*ilo = 1;
	*ihi = 1;
	lscale[1] = 1.;
	rscale[1] = 1.;
	return 0;
    }

    if (template_blas_lsame(job, "S")) {
	goto L190;
    }

    goto L30;

/*     Permute the matrices A and B to isolate the eigenvalues.   

       Find row with one nonzero in columns 1 through L */

L20:
    l = lm1;
    if (l != 1) {
	goto L30;
    }

    rscale[1] = 1.;
    lscale[1] = 1.;
    goto L190;

L30:
    lm1 = l - 1;
    for (i__ = l; i__ >= 1; --i__) {
	i__1 = lm1;
	for (j = 1; j <= i__1; ++j) {
	    jp1 = j + 1;
	    if (a_ref(i__, j) != 0. || b_ref(i__, j) != 0.) {
		goto L50;
	    }
/* L40: */
	}
	j = l;
	goto L70;

L50:
	i__1 = l;
	for (j = jp1; j <= i__1; ++j) {
	    if (a_ref(i__, j) != 0. || b_ref(i__, j) != 0.) {
		goto L80;
	    }
/* L60: */
	}
	j = jp1 - 1;

L70:
	m = l;
	iflow = 1;
	goto L160;
L80:
	;
    }
    goto L100;

/*     Find column with one nonzero in rows K through N */

L90:
    ++k;

L100:
    i__1 = l;
    for (j = k; j <= i__1; ++j) {
	i__2 = lm1;
	for (i__ = k; i__ <= i__2; ++i__) {
	    ip1 = i__ + 1;
	    if (a_ref(i__, j) != 0. || b_ref(i__, j) != 0.) {
		goto L120;
	    }
/* L110: */
	}
	i__ = l;
	goto L140;
L120:
	i__2 = l;
	for (i__ = ip1; i__ <= i__2; ++i__) {
	    if (a_ref(i__, j) != 0. || b_ref(i__, j) != 0.) {
		goto L150;
	    }
/* L130: */
	}
	i__ = ip1 - 1;
L140:
	m = k;
	iflow = 2;
	goto L160;
L150:
	;
    }
    goto L190;

/*     Permute rows M and I */

L160:
    lscale[m] = (Treal) i__;
    if (i__ == m) {
	goto L170;
    }
    i__1 = *n - k + 1;
    template_blas_swap(&i__1, &a_ref(i__, k), lda, &a_ref(m, k), lda);
    i__1 = *n - k + 1;
    template_blas_swap(&i__1, &b_ref(i__, k), ldb, &b_ref(m, k), ldb);

/*     Permute columns M and J */

L170:
    rscale[m] = (Treal) j;
    if (j == m) {
	goto L180;
    }
    template_blas_swap(&l, &a_ref(1, j), &c__1, &a_ref(1, m), &c__1);
    template_blas_swap(&l, &b_ref(1, j), &c__1, &b_ref(1, m), &c__1);

L180:
    switch (iflow) {
	case 1:  goto L20;
	case 2:  goto L90;
    }

L190:
    *ilo = k;
    *ihi = l;

    if (*ilo == *ihi) {
	return 0;
    }

    if (template_blas_lsame(job, "P")) {
	return 0;
    }

/*     Balance the submatrix in rows ILO to IHI. */

    nr = *ihi - *ilo + 1;
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	rscale[i__] = 0.;
	lscale[i__] = 0.;

	work[i__] = 0.;
	work[i__ + *n] = 0.;
	work[i__ + (*n << 1)] = 0.;
	work[i__ + *n * 3] = 0.;
	work[i__ + (*n << 2)] = 0.;
	work[i__ + *n * 5] = 0.;
/* L200: */
    }

/*     Compute right side vector in resulting linear equations */

    basl = template_blas_lg10(&c_b34);
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	i__2 = *ihi;
	for (j = *ilo; j <= i__2; ++j) {
	    tb = b_ref(i__, j);
	    ta = a_ref(i__, j);
	    if (ta == 0.) {
		goto L210;
	    }
	    d__1 = absMACRO(ta);
	    ta = template_blas_lg10(&d__1) / basl;
L210:
	    if (tb == 0.) {
		goto L220;
	    }
	    d__1 = absMACRO(tb);
	    tb = template_blas_lg10(&d__1) / basl;
L220:
	    work[i__ + (*n << 2)] = work[i__ + (*n << 2)] - ta - tb;
	    work[j + *n * 5] = work[j + *n * 5] - ta - tb;
/* L230: */
	}
/* L240: */
    }

    coef = 1. / (Treal) (nr << 1);
    coef2 = coef * coef;
    coef5 = coef2 * .5;
    nrp2 = nr + 2;
    beta = 0.;
    it = 1;

/*     Start generalized conjugate gradient iteration */

L250:

    gamma = template_blas_dot(&nr, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1) + template_blas_dot(&nr, &work[*ilo + *n * 5], &c__1, &work[*ilo + *
	    n * 5], &c__1);

    ew = 0.;
    ewc = 0.;
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	ew += work[i__ + (*n << 2)];
	ewc += work[i__ + *n * 5];
/* L260: */
    }

/* Computing 2nd power */
    d__1 = ew;
/* Computing 2nd power */
    d__2 = ewc;
/* Computing 2nd power */
    d__3 = ew - ewc;
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
    if (gamma == 0.) {
	goto L350;
    }
    if (it != 1) {
	beta = gamma / pgamma;
    }
    t = coef5 * (ewc - ew * 3.);
    tc = coef5 * (ew - ewc * 3.);

    template_blas_scal(&nr, &beta, &work[*ilo], &c__1);
    template_blas_scal(&nr, &beta, &work[*ilo + *n], &c__1);

    template_blas_axpy(&nr, &coef, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + *n], &
	    c__1);
    template_blas_axpy(&nr, &coef, &work[*ilo + *n * 5], &c__1, &work[*ilo], &c__1);

    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	work[i__] += tc;
	work[i__ + *n] += t;
/* L270: */
    }

/*     Apply matrix to vector */

    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	kount = 0;
	sum = 0.;
	i__2 = *ihi;
	for (j = *ilo; j <= i__2; ++j) {
	    if (a_ref(i__, j) == 0.) {
		goto L280;
	    }
	    ++kount;
	    sum += work[j];
L280:
	    if (b_ref(i__, j) == 0.) {
		goto L290;
	    }
	    ++kount;
	    sum += work[j];
L290:
	    ;
	}
	work[i__ + (*n << 1)] = (Treal) kount * work[i__ + *n] + sum;
/* L300: */
    }

    i__1 = *ihi;
    for (j = *ilo; j <= i__1; ++j) {
	kount = 0;
	sum = 0.;
	i__2 = *ihi;
	for (i__ = *ilo; i__ <= i__2; ++i__) {
	    if (a_ref(i__, j) == 0.) {
		goto L310;
	    }
	    ++kount;
	    sum += work[i__ + *n];
L310:
	    if (b_ref(i__, j) == 0.) {
		goto L320;
	    }
	    ++kount;
	    sum += work[i__ + *n];
L320:
	    ;
	}
	work[j + *n * 3] = (Treal) kount * work[j] + sum;
/* L330: */
    }

    sum = template_blas_dot(&nr, &work[*ilo + *n], &c__1, &work[*ilo + (*n << 1)], &c__1) 
	    + template_blas_dot(&nr, &work[*ilo], &c__1, &work[*ilo + *n * 3], &c__1);
    alpha = gamma / sum;

/*     Determine correction to current iteration */

    cmax = 0.;
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	cor = alpha * work[i__ + *n];
	if (absMACRO(cor) > cmax) {
	    cmax = absMACRO(cor);
	}
	lscale[i__] += cor;
	cor = alpha * work[i__];
	if (absMACRO(cor) > cmax) {
	    cmax = absMACRO(cor);
	}
	rscale[i__] += cor;
/* L340: */
    }
    if (cmax < .5) {
	goto L350;
    }

    d__1 = -alpha;
    template_blas_axpy(&nr, &d__1, &work[*ilo + (*n << 1)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1);
    d__1 = -alpha;
    template_blas_axpy(&nr, &d__1, &work[*ilo + *n * 3], &c__1, &work[*ilo + *n * 5], &
	    c__1);

    pgamma = gamma;
    ++it;
    if (it <= nrp2) {
	goto L250;
    }

/*     End generalized conjugate gradient iteration */

L350:
    sfmin = template_lapack_lamch("S", (Treal)0);
    sfmax = 1. / sfmin;
    lsfmin = (integer) (template_blas_lg10(&sfmin) / basl + 1.);
    lsfmax = (integer) (template_blas_lg10(&sfmax) / basl);
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	i__2 = *n - *ilo + 1;
	irab = template_blas_idamax(&i__2, &a_ref(i__, *ilo), lda);
	rab = (d__1 = a_ref(i__, irab + *ilo - 1), absMACRO(d__1));
	i__2 = *n - *ilo + 1;
	irab = template_blas_idamax(&i__2, &b_ref(i__, *ilo), lda);
/* Computing MAX */
	d__2 = rab, d__3 = (d__1 = b_ref(i__, irab + *ilo - 1), absMACRO(d__1));
	rab = maxMACRO(d__2,d__3);
	d__1 = rab + sfmin;
	lrab = (integer) (template_blas_lg10(&d__1) / basl + 1.);
	ir = (integer) (lscale[i__] + template_lapack_d_sign(&c_b70, &lscale[i__]));
/* Computing MIN */
	i__2 = maxMACRO(ir,lsfmin), i__2 = minMACRO(i__2,lsfmax), i__3 = lsfmax - lrab;
	ir = minMACRO(i__2,i__3);
	lscale[i__] = template_lapack_pow_di(&c_b34, &ir);
	icab = template_blas_idamax(ihi, &a_ref(1, i__), &c__1);
	cab = (d__1 = a_ref(icab, i__), absMACRO(d__1));
	icab = template_blas_idamax(ihi, &b_ref(1, i__), &c__1);
/* Computing MAX */
	d__2 = cab, d__3 = (d__1 = b_ref(icab, i__), absMACRO(d__1));
	cab = maxMACRO(d__2,d__3);
	d__1 = cab + sfmin;
	lcab = (integer) (template_blas_lg10(&d__1) / basl + 1.);
	jc = (integer) (rscale[i__] + template_lapack_d_sign(&c_b70, &rscale[i__]));
/* Computing MIN */
	i__2 = maxMACRO(jc,lsfmin), i__2 = minMACRO(i__2,lsfmax), i__3 = lsfmax - lcab;
	jc = minMACRO(i__2,i__3);
	rscale[i__] = template_lapack_pow_di(&c_b34, &jc);
/* L360: */
    }

/*     Row scaling of matrices A and B */

    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	i__2 = *n - *ilo + 1;
	template_blas_scal(&i__2, &lscale[i__], &a_ref(i__, *ilo), lda);
	i__2 = *n - *ilo + 1;
	template_blas_scal(&i__2, &lscale[i__], &b_ref(i__, *ilo), ldb);
/* L370: */
    }

/*     Column scaling of matrices A and B */

    i__1 = *ihi;
    for (j = *ilo; j <= i__1; ++j) {
	template_blas_scal(ihi, &rscale[j], &a_ref(1, j), &c__1);
	template_blas_scal(ihi, &rscale[j], &b_ref(1, j), &c__1);
/* L380: */
    }

    return 0;

/*     End of DGGBAL */

} /* dggbal_ */

#undef b_ref
#undef a_ref


#endif
