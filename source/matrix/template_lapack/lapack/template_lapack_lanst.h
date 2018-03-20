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
 

#ifndef TEMPLATE_LAPACK_LANST_HEADER
#define TEMPLATE_LAPACK_LANST_HEADER


template<class Treal>
Treal template_lapack_lanst(const char *norm, const integer *n, const Treal *d__, const Treal *e)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLANST  returns the value of the one norm,  or the Frobenius norm, or   
    the  infinity norm,  or the  element of  largest absolute value  of a   
    real symmetric tridiagonal matrix A.   

    Description   
    ===========   

    DLANST returns the value   

       DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum),   
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
    normF  denotes the  Frobenius norm of a matrix (square root of sum of   
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANST as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANST is   
            set to zero.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of A.   

    E       (input) DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) sub-diagonal or super-diagonal elements of A.   

    =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
     integer c__1 = 1;
    
    /* System generated locals */
    integer i__1;
    Treal ret_val, d__1, d__2, d__3, d__4, d__5;
    /* Local variables */
     integer i__;
     Treal scale;
     Treal anorm;
     Treal sum;


    --e;
    --d__;

    /* Initialization added by Elias to get rid of compiler warnings. */
    anorm = 0;
    /* Function Body */
    if (*n <= 0) {
	anorm = 0.;
    } else if (template_blas_lsame(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	anorm = (d__1 = d__[*n], absMACRO(d__1));
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = d__[i__], absMACRO(d__1));
	    anorm = maxMACRO(d__2,d__3);
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = e[i__], absMACRO(d__1));
	    anorm = maxMACRO(d__2,d__3);
/* L10: */
	}
    } else if (template_blas_lsame(norm, "O") || *(unsigned char *)
	    norm == '1' || template_blas_lsame(norm, "I")) {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = absMACRO(d__[1]);
	} else {
/* Computing MAX */
	    d__3 = absMACRO(d__[1]) + absMACRO(e[1]), d__4 = (d__1 = e[*n - 1], absMACRO(
		    d__1)) + (d__2 = d__[*n], absMACRO(d__2));
	    anorm = maxMACRO(d__3,d__4);
	    i__1 = *n - 1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__4 = anorm, d__5 = (d__1 = d__[i__], absMACRO(d__1)) + (d__2 = e[
			i__], absMACRO(d__2)) + (d__3 = e[i__ - 1], absMACRO(d__3));
		anorm = maxMACRO(d__4,d__5);
/* L20: */
	    }
	}
    } else if (template_blas_lsame(norm, "F") || template_blas_lsame(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	if (*n > 1) {
	    i__1 = *n - 1;
	    template_lapack_lassq(&i__1, &e[1], &c__1, &scale, &sum);
	    sum *= 2;
	}
	template_lapack_lassq(n, &d__[1], &c__1, &scale, &sum);
	anorm = scale * template_blas_sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of DLANST */

} /* dlanst_ */

#endif
