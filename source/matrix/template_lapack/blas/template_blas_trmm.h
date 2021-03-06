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
 

#ifndef TEMPLATE_BLAS_TRMM_HEADER
#define TEMPLATE_BLAS_TRMM_HEADER

#include "template_blas_common.h"

template<class Treal>
int template_blas_trmm(const char *side, const char *uplo, const char *transa, const char *diag, 
	const integer *m, const integer *n, const Treal *alpha, const Treal *a, const integer *
	lda, Treal *b, const integer *ldb)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
     integer info;
     Treal temp;
     integer i__, j, k;
     logical lside;
     integer nrowa;
     logical upper;
     logical nounit;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
/*  Purpose   
    =======   
    DTRMM  performs one of the matrix-matrix operations   
       B := alpha*op( A )*B,   or   B := alpha*B*op( A ),   
    where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or   
    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of   
       op( A ) = A   or   op( A ) = A'.   
    Parameters   
    ==========   
    SIDE   - CHARACTER*1.   
             On entry,  SIDE specifies whether  op( A ) multiplies B from   
             the left or right as follows:   
                SIDE = 'L' or 'l'   B := alpha*op( A )*B.   
                SIDE = 'R' or 'r'   B := alpha*B*op( A ).   
             Unchanged on exit.   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix A is an upper or   
             lower triangular matrix as follows:   
                UPLO = 'U' or 'u'   A is an upper triangular matrix.   
                UPLO = 'L' or 'l'   A is a lower triangular matrix.   
             Unchanged on exit.   
    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in   
             the matrix multiplication as follows:   
                TRANSA = 'N' or 'n'   op( A ) = A.   
                TRANSA = 'T' or 't'   op( A ) = A'.   
                TRANSA = 'C' or 'c'   op( A ) = A'.   
             Unchanged on exit.   
    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit triangular   
             as follows:   
                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry, M specifies the number of rows of B. M must be at   
             least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of B.  N must be   
             at least zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry,  ALPHA specifies the scalar  alpha. When  alpha is   
             zero then  A is not referenced and  B need not be set before   
             entry.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m   
             when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.   
             Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k   
             upper triangular part of the array  A must contain the upper   
             triangular matrix  and the strictly lower triangular part of   
             A is not referenced.   
             Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k   
             lower triangular part of the array  A must contain the lower   
             triangular matrix  and the strictly upper triangular part of   
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u',  the diagonal elements of   
             A  are not referenced either,  but are assumed to be  unity.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then   
             LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'   
             then LDA must be at least max( 1, n ).   
             Unchanged on exit.   
    B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
             Before entry,  the leading  m by n part of the array  B must   
             contain the matrix  B,  and  on exit  is overwritten  by the   
             transformed matrix.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in  the  calling  (sub)  program.   LDB  must  be  at  least   
             max( 1, m ).   
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
    /* Function Body */
    lside = template_blas_lsame(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = template_blas_lsame(diag, "N");
    upper = template_blas_lsame(uplo, "U");
    info = 0;
    if (! lside && ! template_blas_lsame(side, "R")) {
	info = 1;
    } else if (! upper && ! template_blas_lsame(uplo, "L")) {
	info = 2;
    } else if (! template_blas_lsame(transa, "N") && ! template_blas_lsame(transa,
	     "T") && ! template_blas_lsame(transa, "C")) {
	info = 3;
    } else if (! template_blas_lsame(diag, "U") && ! template_blas_lsame(diag, 
	    "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < maxMACRO(1,nrowa)) {
	info = 9;
    } else if (*ldb < maxMACRO(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	template_blas_erbla("TRMM  ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0) {
	return 0;
    }
/*     And when  alpha.eq.zero. */
    if (*alpha == 0.) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b_ref(i__, j) = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }
/*     Start the operations. */
    if (lside) {
	if (template_blas_lsame(transa, "N")) {
/*           Form  B := alpha*A*B. */
	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b_ref(k, j) != 0.) {
			    temp = *alpha * b_ref(k, j);
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * a_ref(
					i__, k);
/* L30: */
			    }
			    if (nounit) {
				temp *= a_ref(k, k);
			    }
			    b_ref(k, j) = temp;
			}
/* L40: */
		    }
/* L50: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = *m; k >= 1; --k) {
			if (b_ref(k, j) != 0.) {
			    temp = *alpha * b_ref(k, j);
			    b_ref(k, j) = temp;
			    if (nounit) {
				b_ref(k, j) = b_ref(k, j) * a_ref(k, k);
			    }
			    i__2 = *m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * a_ref(
					i__, k);
/* L60: */
			    }
			}
/* L70: */
		    }
/* L80: */
		}
	    }
	} else {
/*           Form  B := alpha*A'*B. */
	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = b_ref(i__, j);
			if (nounit) {
			    temp *= a_ref(i__, i__);
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a_ref(k, i__) * b_ref(k, j);
/* L90: */
			}
			b_ref(i__, j) = *alpha * temp;
/* L100: */
		    }
/* L110: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b_ref(i__, j);
			if (nounit) {
			    temp *= a_ref(i__, i__);
			}
			i__3 = *m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a_ref(k, i__) * b_ref(k, j);
/* L120: */
			}
			b_ref(i__, j) = *alpha * temp;
/* L130: */
		    }
/* L140: */
		}
	    }
	}
    } else {
	if (template_blas_lsame(transa, "N")) {
/*           Form  B := alpha*B*A. */
	    if (upper) {
		for (j = *n; j >= 1; --j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b_ref(i__, j) = temp * b_ref(i__, j);
/* L150: */
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if (a_ref(k, j) != 0.) {
			    temp = *alpha * a_ref(k, j);
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L160: */
			    }
			}
/* L170: */
		    }
/* L180: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b_ref(i__, j) = temp * b_ref(i__, j);
/* L190: */
		    }
		    i__2 = *n;
		    for (k = j + 1; k <= i__2; ++k) {
			if (a_ref(k, j) != 0.) {
			    temp = *alpha * a_ref(k, j);
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L200: */
			    }
			}
/* L210: */
		    }
/* L220: */
		}
	    }
	} else {
/*           Form  B := alpha*B*A'. */
	    if (upper) {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (a_ref(j, k) != 0.) {
			    temp = *alpha * a_ref(j, k);
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L230: */
			    }
			}
/* L240: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(k, k);
		    }
		    if (temp != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b_ref(i__, k) = temp * b_ref(i__, k);
/* L250: */
			}
		    }
/* L260: */
		}
	    } else {
		for (k = *n; k >= 1; --k) {
		    i__1 = *n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (a_ref(j, k) != 0.) {
			    temp = *alpha * a_ref(j, k);
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L270: */
			    }
			}
/* L280: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(k, k);
		    }
		    if (temp != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b_ref(i__, k) = temp * b_ref(i__, k);
/* L290: */
			}
		    }
/* L300: */
		}
	    }
	}
    }
    return 0;
/*     End of DTRMM . */
} /* dtrmm_ */
#undef b_ref
#undef a_ref

#endif
