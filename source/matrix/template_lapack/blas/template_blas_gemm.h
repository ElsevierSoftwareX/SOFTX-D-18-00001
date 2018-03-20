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
 

#ifndef TEMPLATE_BLAS_GEMM_HEADER
#define TEMPLATE_BLAS_GEMM_HEADER

#include "template_blas_common.h"

template<class Treal>
int template_blas_gemm(const char *transa, const char *transb, const integer *m, const integer *
	n, const integer *k, const Treal *alpha, const Treal *a, const integer *lda, 
	const Treal *b, const integer *ldb, const Treal *beta, Treal *c__, 
	const integer *ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    /* Local variables */
     integer info;
     logical nota, notb;
     Treal temp;
     integer i__, j, l;
     integer nrowa, nrowb;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
/*  Purpose   
    =======   
    DGEMM  performs one of the matrix-matrix operations   
       C := alpha*op( A )*op( B ) + beta*C,   
    where  op( X ) is one of   
       op( X ) = X   or   op( X ) = X',   
    alpha and beta are scalars, and A, B and C are matrices, with op( A )   
    an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.   
    Parameters   
    ==========   
    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in   
             the matrix multiplication as follows:   
                TRANSA = 'N' or 'n',  op( A ) = A.   
                TRANSA = 'T' or 't',  op( A ) = A'.   
                TRANSA = 'C' or 'c',  op( A ) = A'.   
             Unchanged on exit.   
    TRANSB - CHARACTER*1.   
             On entry, TRANSB specifies the form of op( B ) to be used in   
             the matrix multiplication as follows:   
                TRANSB = 'N' or 'n',  op( B ) = B.   
                TRANSB = 'T' or 't',  op( B ) = B'.   
                TRANSB = 'C' or 'c',  op( B ) = B'.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry,  M  specifies  the number  of rows  of the  matrix   
             op( A )  and of the  matrix  C.  M  must  be at least  zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry,  N  specifies the number  of columns of the matrix   
             op( B ) and the number of columns of the matrix C. N must be   
             at least zero.   
             Unchanged on exit.   
    K      - INTEGER.   
             On entry,  K  specifies  the number of columns of the matrix   
             op( A ) and the number of rows of the matrix op( B ). K must   
             be at least  zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is   
             k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.   
             Before entry with  TRANSA = 'N' or 'n',  the leading  m by k   
             part of the array  A  must contain the matrix  A,  otherwise   
             the leading  k by m  part of the array  A  must contain  the   
             matrix A.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. When  TRANSA = 'N' or 'n' then   
             LDA must be at least  max( 1, m ), otherwise  LDA must be at   
             least  max( 1, k ).   
             Unchanged on exit.   
    B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is   
             n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.   
             Before entry with  TRANSB = 'N' or 'n',  the leading  k by n   
             part of the array  B  must contain the matrix  B,  otherwise   
             the leading  n by k  part of the array  B  must contain  the   
             matrix B.   
             Unchanged on exit.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in the calling (sub) program. When  TRANSB = 'N' or 'n' then   
             LDB must be at least  max( 1, k ), otherwise  LDB must be at   
             least  max( 1, n ).   
             Unchanged on exit.   
    BETA   - DOUBLE PRECISION.   
             On entry,  BETA  specifies the scalar  beta.  When  BETA  is   
             supplied as zero then C need not be set on input.   
             Unchanged on exit.   
    C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry, the leading  m by n  part of the array  C must   
             contain the matrix  C,  except when  beta  is zero, in which   
             case C need not be set on entry.   
             On exit, the array  C  is overwritten by the  m by n  matrix   
             ( alpha*op( A )*op( B ) + beta*C ).   
    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared   
             in  the  calling  (sub)  program.   LDC  must  be  at  least   
             max( 1, m ).   
             Unchanged on exit.   
    Level 3 Blas routine.   
    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   
       Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not   
       transposed and set  NROWA, NCOLA and  NROWB  as the number of rows   
       and  columns of  A  and the  number of  rows  of  B  respectively.   
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
    nota = template_blas_lsame(transa, "N");
    notb = template_blas_lsame(transb, "N");
    if (nota) {
	nrowa = *m;
    } else {
	nrowa = *k;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }
/*     Test the input parameters. */
    info = 0;
    if (! nota && ! template_blas_lsame(transa, "C") && ! template_blas_lsame(
	    transa, "T")) {
	info = 1;
    } else if (! notb && ! template_blas_lsame(transb, "C") && ! 
	    template_blas_lsame(transb, "T")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < maxMACRO(1,nrowa)) {
	info = 8;
    } else if (*ldb < maxMACRO(1,nrowb)) {
	info = 10;
    } else if (*ldc < maxMACRO(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	template_blas_erbla("DGEMM ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*m == 0 || *n == 0 || ( (*alpha == 0. || *k == 0) && *beta == 1.) ) {
	return 0;
    }
/*     And if  alpha.eq.zero. */
    if (*alpha == 0.) {
	if (*beta == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = *beta * c___ref(i__, j);
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }
/*     Start the operations. */
    if (notb) {
	if (nota) {
/*           Form  C := alpha*A*B + beta*C. */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.;
/* L50: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b_ref(l, j) != 0.) {
			temp = *alpha * b_ref(l, j);
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + temp * a_ref(
				    i__, l);
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else {
/*           Form  C := alpha*A'*B + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a_ref(l, i__) * b_ref(l, j);
/* L100: */
		    }
		    if (*beta == 0.) {
			c___ref(i__, j) = *alpha * temp;
		    } else {
			c___ref(i__, j) = *alpha * temp + *beta * c___ref(i__,
				 j);
		    }
/* L110: */
		}
/* L120: */
	    }
	}
    } else {
	if (nota) {
/*           Form  C := alpha*A*B' + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.;
/* L130: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L140: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b_ref(j, l) != 0.) {
			temp = *alpha * b_ref(j, l);
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + temp * a_ref(
				    i__, l);
/* L150: */
			}
		    }
/* L160: */
		}
/* L170: */
	    }
	} else {
/*           Form  C := alpha*A'*B' + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a_ref(l, i__) * b_ref(j, l);
/* L180: */
		    }
		    if (*beta == 0.) {
			c___ref(i__, j) = *alpha * temp;
		    } else {
			c___ref(i__, j) = *alpha * temp + *beta * c___ref(i__,
				 j);
		    }
/* L190: */
		}
/* L200: */
	    }
	}
    }
    return 0;
/*     End of DGEMM . */
} /* dgemm_ */
#undef c___ref
#undef b_ref
#undef a_ref

#endif
