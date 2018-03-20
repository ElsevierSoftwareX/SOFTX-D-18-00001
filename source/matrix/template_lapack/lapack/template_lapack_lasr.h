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
 

#ifndef TEMPLATE_LAPACK_LASR_HEADER
#define TEMPLATE_LAPACK_LASR_HEADER


template<class Treal>
int template_lapack_lasr(const char *side, const char *pivot, const char *direct, const integer *m,
	 const integer *n, const Treal *c__, const Treal *s, Treal *a, const integer *
	lda)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASR   performs the transformation   

       A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )   

       A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )   

    where A is an m by n real matrix and P is an orthogonal matrix,   
    consisting of a sequence of plane rotations determined by the   
    parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'   
    and z = n when SIDE = 'R' or 'r' ):   

    When  DIRECT = 'F' or 'f'  ( Forward sequence ) then   

       P = P( z - 1 )*...*P( 2 )*P( 1 ),   

    and when DIRECT = 'B' or 'b'  ( Backward sequence ) then   

       P = P( 1 )*P( 2 )*...*P( z - 1 ),   

    where  P( k ) is a plane rotation matrix for the following planes:   

       when  PIVOT = 'V' or 'v'  ( Variable pivot ),   
          the plane ( k, k + 1 )   

       when  PIVOT = 'T' or 't'  ( Top pivot ),   
          the plane ( 1, k + 1 )   

       when  PIVOT = 'B' or 'b'  ( Bottom pivot ),   
          the plane ( k, z )   

    c( k ) and s( k )  must contain the  cosine and sine that define the   
    matrix  P( k ).  The two by two plane rotation part of the matrix   
    P( k ), R( k ), is assumed to be of the form   

       R( k ) = (  c( k )  s( k ) ).   
                ( -s( k )  c( k ) )   

    This version vectorises across rows of the array A when SIDE = 'L'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            Specifies whether the plane rotation matrix P is applied to   
            A on the left or the right.   
            = 'L':  Left, compute A := P*A   
            = 'R':  Right, compute A:= A*P'   

    DIRECT  (input) CHARACTER*1   
            Specifies whether P is a forward or backward sequence of   
            plane rotations.   
            = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )   
            = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )   

    PIVOT   (input) CHARACTER*1   
            Specifies the plane for which P(k) is a plane rotation   
            matrix.   
            = 'V':  Variable pivot, the plane (k,k+1)   
            = 'T':  Top pivot, the plane (1,k+1)   
            = 'B':  Bottom pivot, the plane (k,z)   

    M       (input) INTEGER   
            The number of rows of the matrix A.  If m <= 1, an immediate   
            return is effected.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  If n <= 1, an   
            immediate return is effected.   

    C, S    (input) DOUBLE PRECISION arrays, dimension   
                    (M-1) if SIDE = 'L'   
                    (N-1) if SIDE = 'R'   
            c(k) and s(k) contain the cosine and sine that define the   
            matrix P(k).  The two by two plane rotation part of the   
            matrix P(k), R(k), is assumed to be of the form   
            R( k ) = (  c( k )  s( k ) ).   
                     ( -s( k )  c( k ) )   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            The m by n matrix A.  On exit, A is overwritten by P*A if   
            SIDE = 'R' or by A*P' if SIDE = 'L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
     integer info;
     Treal temp;
     integer i__, j;
     Treal ctemp, stemp;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! (template_blas_lsame(side, "L") || template_blas_lsame(side, "R"))) {
	info = 1;
    } else if (! (template_blas_lsame(pivot, "V") || template_blas_lsame(pivot, 
	    "T") || template_blas_lsame(pivot, "B"))) {
	info = 2;
    } else if (! (template_blas_lsame(direct, "F") || template_blas_lsame(direct, 
	    "B"))) {
	info = 3;
    } else if (*m < 0) {
	info = 4;
    } else if (*n < 0) {
	info = 5;
    } else if (*lda < maxMACRO(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	template_blas_erbla("LASR  ", &info);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }
    if (template_blas_lsame(side, "L")) {

/*        Form  P * A */

	if (template_blas_lsame(pivot, "V")) {
	    if (template_blas_lsame(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(j + 1, i__);
			    a_ref(j + 1, i__) = ctemp * temp - stemp * a_ref(
				    j, i__);
			    a_ref(j, i__) = stemp * temp + ctemp * a_ref(j, 
				    i__);
/* L10: */
			}
		    }
/* L20: */
		}
	    } else if (template_blas_lsame(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(j + 1, i__);
			    a_ref(j + 1, i__) = ctemp * temp - stemp * a_ref(
				    j, i__);
			    a_ref(j, i__) = stemp * temp + ctemp * a_ref(j, 
				    i__);
/* L30: */
			}
		    }
/* L40: */
		}
	    }
	} else if (template_blas_lsame(pivot, "T")) {
	    if (template_blas_lsame(direct, "F")) {
		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(j, i__);
			    a_ref(j, i__) = ctemp * temp - stemp * a_ref(1, 
				    i__);
			    a_ref(1, i__) = stemp * temp + ctemp * a_ref(1, 
				    i__);
/* L50: */
			}
		    }
/* L60: */
		}
	    } else if (template_blas_lsame(direct, "B")) {
		for (j = *m; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(j, i__);
			    a_ref(j, i__) = ctemp * temp - stemp * a_ref(1, 
				    i__);
			    a_ref(1, i__) = stemp * temp + ctemp * a_ref(1, 
				    i__);
/* L70: */
			}
		    }
/* L80: */
		}
	    }
	} else if (template_blas_lsame(pivot, "B")) {
	    if (template_blas_lsame(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(j, i__);
			    a_ref(j, i__) = stemp * a_ref(*m, i__) + ctemp * 
				    temp;
			    a_ref(*m, i__) = ctemp * a_ref(*m, i__) - stemp * 
				    temp;
/* L90: */
			}
		    }
/* L100: */
		}
	    } else if (template_blas_lsame(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(j, i__);
			    a_ref(j, i__) = stemp * a_ref(*m, i__) + ctemp * 
				    temp;
			    a_ref(*m, i__) = ctemp * a_ref(*m, i__) - stemp * 
				    temp;
/* L110: */
			}
		    }
/* L120: */
		}
	    }
	}
    } else if (template_blas_lsame(side, "R")) {

/*        Form A * P' */

	if (template_blas_lsame(pivot, "V")) {
	    if (template_blas_lsame(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(i__, j + 1);
			    a_ref(i__, j + 1) = ctemp * temp - stemp * a_ref(
				    i__, j);
			    a_ref(i__, j) = stemp * temp + ctemp * a_ref(i__, 
				    j);
/* L130: */
			}
		    }
/* L140: */
		}
	    } else if (template_blas_lsame(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(i__, j + 1);
			    a_ref(i__, j + 1) = ctemp * temp - stemp * a_ref(
				    i__, j);
			    a_ref(i__, j) = stemp * temp + ctemp * a_ref(i__, 
				    j);
/* L150: */
			}
		    }
/* L160: */
		}
	    }
	} else if (template_blas_lsame(pivot, "T")) {
	    if (template_blas_lsame(direct, "F")) {
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(i__, j);
			    a_ref(i__, j) = ctemp * temp - stemp * a_ref(i__, 
				    1);
			    a_ref(i__, 1) = stemp * temp + ctemp * a_ref(i__, 
				    1);
/* L170: */
			}
		    }
/* L180: */
		}
	    } else if (template_blas_lsame(direct, "B")) {
		for (j = *n; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(i__, j);
			    a_ref(i__, j) = ctemp * temp - stemp * a_ref(i__, 
				    1);
			    a_ref(i__, 1) = stemp * temp + ctemp * a_ref(i__, 
				    1);
/* L190: */
			}
		    }
/* L200: */
		}
	    }
	} else if (template_blas_lsame(pivot, "B")) {
	    if (template_blas_lsame(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a_ref(i__, j);
			    a_ref(i__, j) = stemp * a_ref(i__, *n) + ctemp * 
				    temp;
			    a_ref(i__, *n) = ctemp * a_ref(i__, *n) - stemp * 
				    temp;
/* L210: */
			}
		    }
/* L220: */
		}
	    } else if (template_blas_lsame(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a_ref(i__, j);
			    a_ref(i__, j) = stemp * a_ref(i__, *n) + ctemp * 
				    temp;
			    a_ref(i__, *n) = ctemp * a_ref(i__, *n) - stemp * 
				    temp;
/* L230: */
			}
		    }
/* L240: */
		}
	    }
	}
    }

    return 0;

/*     End of DLASR */

} /* dlasr_ */

#undef a_ref


#endif
