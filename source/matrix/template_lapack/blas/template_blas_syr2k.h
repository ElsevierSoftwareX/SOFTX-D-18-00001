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
 

#ifndef TEMPLATE_BLAS_SYR2K_HEADER
#define TEMPLATE_BLAS_SYR2K_HEADER


template<class Treal>
int template_blas_syr2k(const char *uplo, const char *trans, const integer *n, 
			const integer *k, const Treal *alpha, const Treal *a, 
			const integer *lda, const Treal *b, const integer *ldb, 
			const Treal *beta, Treal *c__, const integer *ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    /* Local variables */
     integer info;
     Treal temp1, temp2;
     integer i__, j, l;
     integer nrowa;
     logical upper;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
/*  Purpose   
    =======   
    DSYR2K  performs one of the symmetric rank 2k operations   
       C := alpha*A*B' + alpha*B*A' + beta*C,   
    or   
       C := alpha*A'*B + alpha*B'*A + beta*C,   
    where  alpha and beta  are scalars, C is an  n by n  symmetric matrix   
    and  A and B  are  n by k  matrices  in the  first  case  and  k by n   
    matrices in the second case.   
    Parameters   
    ==========   
    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower   
             triangular  part  of the  array  C  is to be  referenced  as   
             follows:   
                UPLO = 'U' or 'u'   Only the  upper triangular part of  C   
                                    is to be referenced.   
                UPLO = 'L' or 'l'   Only the  lower triangular part of  C   
                                    is to be referenced.   
             Unchanged on exit.   
    TRANS  - CHARACTER*1.   
             On entry,  TRANS  specifies the operation to be performed as   
             follows:   
                TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +   
                                          beta*C.   
                TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +   
                                          beta*C.   
                TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +   
                                          beta*C.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry,  N specifies the order of the matrix C.  N must be   
             at least zero.   
             Unchanged on exit.   
    K      - INTEGER.   
             On entry with  TRANS = 'N' or 'n',  K  specifies  the number   
             of  columns  of the  matrices  A and B,  and on  entry  with   
             TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number   
             of rows of the matrices  A and B.  K must be at least  zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is   
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k   
             part of the array  A  must contain the matrix  A,  otherwise   
             the leading  k by n  part of the array  A  must contain  the   
             matrix A.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'   
             then  LDA must be at least  max( 1, n ), otherwise  LDA must   
             be at least  max( 1, k ).   
             Unchanged on exit.   
    B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is   
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k   
             part of the array  B  must contain the matrix  B,  otherwise   
             the leading  k by n  part of the array  B  must contain  the   
             matrix B.   
             Unchanged on exit.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'   
             then  LDB must be at least  max( 1, n ), otherwise  LDB must   
             be at least  max( 1, k ).   
             Unchanged on exit.   
    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta.   
             Unchanged on exit.   
    C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry  with  UPLO = 'U' or 'u',  the leading  n by n   
             upper triangular part of the array C must contain the upper   
             triangular part  of the  symmetric matrix  and the strictly   
             lower triangular part of C is not referenced.  On exit, the   
             upper triangular part of the array  C is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry  with  UPLO = 'L' or 'l',  the leading  n by n   
             lower triangular part of the array C must contain the lower   
             triangular part  of the  symmetric matrix  and the strictly   
             upper triangular part of C is not referenced.  On exit, the   
             lower triangular part of the array  C is overwritten by the   
             lower triangular part of the updated matrix.   
    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared   
             in  the  calling  (sub)  program.   LDC  must  be  at  least   
             max( 1, n ).   
             Unchanged on exit.   
    Level 3 Blas routine.   
    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    /* Function Body */
    if (template_blas_lsame(trans, "N")) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = template_blas_lsame(uplo, "U");
    info = 0;
    if (! upper && ! template_blas_lsame(uplo, "L")) {
	info = 1;
    } else if (! template_blas_lsame(trans, "N") && ! template_blas_lsame(trans, 
	    "T") && ! template_blas_lsame(trans, "C")) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < maxMACRO(1,nrowa)) {
	info = 7;
    } else if (*ldb < maxMACRO(1,nrowa)) {
	info = 9;
    } else if (*ldc < maxMACRO(1,*n)) {
	info = 12;
    }
    if (info != 0) {
	template_blas_erbla("SYR2K ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0 || ( (*alpha == 0. || *k == 0) && *beta == 1. ) ) {
	return 0;
    }
/*     And when  alpha.eq.zero. */
    if (*alpha == 0.) {
	if (upper) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }
/*     Start the operations. */
    if (template_blas_lsame(trans, "N")) {
/*        Form  C := alpha*A*B' + alpha*B*A' + C. */
	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.;
/* L90: */
		    }
		} else if (*beta != 1.) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L100: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a_ref(j, l) != 0. || b_ref(j, l) != 0.) {
			temp1 = *alpha * b_ref(j, l);
			temp2 = *alpha * a_ref(j, l);
			i__3 = j;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + a_ref(i__, l) 
				    * temp1 + b_ref(i__, l) * temp2;
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.;
/* L140: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L150: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a_ref(j, l) != 0. || b_ref(j, l) != 0.) {
			temp1 = *alpha * b_ref(j, l);
			temp2 = *alpha * a_ref(j, l);
			i__3 = *n;
			for (i__ = j; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + a_ref(i__, l) 
				    * temp1 + b_ref(i__, l) * temp2;
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {
/*        Form  C := alpha*A'*B + alpha*B'*A + C. */
	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp1 = 0.;
		    temp2 = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp1 += a_ref(l, i__) * b_ref(l, j);
			temp2 += b_ref(l, i__) * a_ref(l, j);
/* L190: */
		    }
		    if (*beta == 0.) {
			c___ref(i__, j) = *alpha * temp1 + *alpha * temp2;
		    } else {
			c___ref(i__, j) = *beta * c___ref(i__, j) + *alpha * 
				temp1 + *alpha * temp2;
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    temp1 = 0.;
		    temp2 = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp1 += a_ref(l, i__) * b_ref(l, j);
			temp2 += b_ref(l, i__) * a_ref(l, j);
/* L220: */
		    }
		    if (*beta == 0.) {
			c___ref(i__, j) = *alpha * temp1 + *alpha * temp2;
		    } else {
			c___ref(i__, j) = *beta * c___ref(i__, j) + *alpha * 
				temp1 + *alpha * temp2;
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }
    return 0;
/*     End of DSYR2K. */
} /* dsyr2k_ */
#undef c___ref
#undef b_ref
#undef a_ref

#endif
