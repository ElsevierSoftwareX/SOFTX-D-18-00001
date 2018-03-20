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
 

#ifndef TEMPLATE_LAPACK_LASWP_HEADER
#define TEMPLATE_LAPACK_LASWP_HEADER


template<class Treal>
int template_lapack_laswp(const integer *n, Treal *a, const integer *lda, const integer 
	*k1, const integer *k2, const integer *ipiv, const integer *incx)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DLASWP performs a series of row interchanges on the matrix A.   
    One row interchange is initiated for each of rows K1 through K2 of A.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of columns of the matrix A.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the matrix of column dimension N to which the row   
            interchanges will be applied.   
            On exit, the permuted matrix.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   

    K1      (input) INTEGER   
            The first element of IPIV for which a row interchange will   
            be done.   

    K2      (input) INTEGER   
            The last element of IPIV for which a row interchange will   
            be done.   

    IPIV    (input) INTEGER array, dimension (M*abs(INCX))   
            The vector of pivot indices.  Only the elements in positions   
            K1 through K2 of IPIV are accessed.   
            IPIV(K) = L implies rows K and L are to be interchanged.   

    INCX    (input) INTEGER   
            The increment between successive values of IPIV.  If IPIV   
            is negative, the pivots are applied in reverse order.   

    Further Details   
    ===============   

    Modified by   
     R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA   

   =====================================================================   


       Interchange row I with row IPIV(I) for each of rows K1 through K2.   

       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
     Treal temp;
     integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    if (*incx > 0) {
	ix0 = *k1;
	i1 = *k1;
	i2 = *k2;
	inc = 1;
    } else if (*incx < 0) {
	ix0 = (1 - *k2) * *incx + 1;
	i1 = *k2;
	i2 = *k1;
	inc = -1;
    } else {
	return 0;
    }

    n32 = *n / 32 << 5;
    if (n32 != 0) {
	i__1 = n32;
	for (j = 1; j <= i__1; j += 32) {
	    ix = ix0;
	    i__2 = i2;
	    i__3 = inc;
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
		    {
		ip = ipiv[ix];
		if (ip != i__) {
		    i__4 = j + 31;
		    for (k = j; k <= i__4; ++k) {
			temp = a_ref(i__, k);
			a_ref(i__, k) = a_ref(ip, k);
			a_ref(ip, k) = temp;
/* L10: */
		    }
		}
		ix += *incx;
/* L20: */
	    }
/* L30: */
	}
    }
    if (n32 != *n) {
	++n32;
	ix = ix0;
	i__1 = i2;
	i__3 = inc;
	for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
	    ip = ipiv[ix];
	    if (ip != i__) {
		i__2 = *n;
		for (k = n32; k <= i__2; ++k) {
		    temp = a_ref(i__, k);
		    a_ref(i__, k) = a_ref(ip, k);
		    a_ref(ip, k) = temp;
/* L40: */
		}
	    }
	    ix += *incx;
/* L50: */
	}
    }

    return 0;

/*     End of DLASWP */

} /* dlaswp_ */

#undef a_ref


#endif
