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
 

#ifndef TEMPLATE_LAPACK_LAMCH_HEADER
#define TEMPLATE_LAPACK_LAMCH_HEADER


#include <stdio.h>
#include <iostream>
#include <limits>



template<class Treal>
Treal template_lapack_d_sign(const Treal *a, const Treal *b)
{
  Treal x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}



// FIXME: ON THE NEXT LINE IS A HARD-CODED CONSTANT THAT NEEDS MORE DIGITS FOR QUAD PRECISION
#define log10e 0.43429448190325182765
template<class Treal>
Treal template_blas_lg10(Treal *x)
{
  return( log10e * template_blas_log(*x) );
}








template<class Treal>
int template_lapack_lassq(const integer *n, const Treal *x, const integer *incx, 
	Treal *scale, Treal *sumsq)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DLASSQ  returns the values  scl  and  smsq  such that   

       ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,   

    where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is   
    assumed to be non-negative and  scl  returns the value   

       scl = max( scale, abs( x( i ) ) ).   

    scale and sumsq must be supplied in SCALE and SUMSQ and   
    scl and smsq are overwritten on SCALE and SUMSQ respectively.   

    The routine makes only one pass through the vector x.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of elements to be used from the vector X.   

    X       (input) DOUBLE PRECISION array, dimension (N)   
            The vector for which a scaled sum of squares is computed.   
               x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.   

    INCX    (input) INTEGER   
            The increment between successive values of the vector X.   
            INCX > 0.   

    SCALE   (input/output) DOUBLE PRECISION   
            On entry, the value  scale  in the equation above.   
            On exit, SCALE is overwritten with  scl , the scaling factor   
            for the sum of squares.   

    SUMSQ   (input/output) DOUBLE PRECISION   
            On entry, the value  sumsq  in the equation above.   
            On exit, SUMSQ is overwritten with  smsq , the basic sum of   
            squares from which  scl  has been factored out.   

   =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer i__1, i__2;
    Treal d__1;
    /* Local variables */
     Treal absxi;
     integer ix;

    --x;

    /* Function Body */
    if (*n > 0) {
	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != 0.) {
   	        absxi = (d__1 = x[ix], absMACRO(d__1));
		if (*scale < absxi) {
/* Computing 2nd power */
		    d__1 = *scale / absxi;
		    *sumsq = *sumsq * (d__1 * d__1) + 1;
		    *scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / *scale;
		    *sumsq += d__1 * d__1;
		}
	    }
/* L10: */
	}
    }
    return 0;

/*     End of DLASSQ */

} /* dlassq_ */





template<class Treal>
double template_lapack_pow_di(Treal *ap, integer *bp)
{
  Treal pow, x;
  integer n;
  unsigned long u;

  pow = 1;
  x = *ap;
  n = *bp;

  if(n != 0)
    {
      if(n < 0)
	{
	  n = -n;
	  x = 1/x;
	}
      for(u = n; ; )
	{
	  if(u & 01)
	    pow *= x;
	  if(u >>= 1)
	    x *= x;
	  else
	    break;
	}
    }
  return(pow);
}




template<class Treal>
Treal template_lapack_lamch(const char *cmach, Treal dummyReal)
{

/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMCH determines double precision machine parameters.   

    Arguments   
    =========   

    CMACH   (input) CHARACTER*1   
            Specifies the value to be returned by DLAMCH:   
            = 'E' or 'e',   DLAMCH := eps   
            = 'S' or 's ,   DLAMCH := sfmin   
            = 'B' or 'b',   DLAMCH := base   
            = 'P' or 'p',   DLAMCH := eps*base   
            = 'N' or 'n',   DLAMCH := t   
            = 'R' or 'r',   DLAMCH := rnd   
            = 'M' or 'm',   DLAMCH := emin   
            = 'U' or 'u',   DLAMCH := rmin   
            = 'L' or 'l',   DLAMCH := emax   
            = 'O' or 'o',   DLAMCH := rmax   

            where   

            eps   = relative machine precision   
            sfmin = safe minimum, such that 1/sfmin does not overflow   
            base  = base of the machine   
            prec  = eps*base   
            t     = number of (base) digits in the mantissa   
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
            emin  = minimum exponent before (gradual) underflow   
            rmin  = underflow threshold - base**(emin-1)   
            emax  = largest exponent before overflow   
            rmax  = overflow threshold  - (base**emax)*(1-eps)   

   ===================================================================== 
*/

   Treal rmach, ret_val;

   /* Initialization added by Elias to get rid of compiler warnings. */
   rmach = 0;
  if (template_blas_lsame(cmach, "E")) { /* Epsilon */
    rmach = template_blas_get_machine_epsilon<Treal>();
  } else if (template_blas_lsame(cmach, "S")) { /* Safe minimum */
    rmach = template_blas_get_num_limit_min<Treal>();
  } else if (template_blas_lsame(cmach, "B")) { /* Base */
    /* Assume "base" is 2 */
    rmach = 2.0;
  } else if (template_blas_lsame(cmach, "P")) { /* Precision */
    /* Assume "base" is 2 */
    rmach = 2.0 * template_blas_get_machine_epsilon<Treal>();
  } else if (template_blas_lsame(cmach, "N")) {
    std::cout << "ERROR in template_lapack_lamch: case N not implemented." << std::endl;
    throw "ERROR in template_lapack_lamch: case N not implemented.";
  } else if (template_blas_lsame(cmach, "R")) {
    std::cout << "ERROR in template_lapack_lamch: case R not implemented." << std::endl;
    throw "ERROR in template_lapack_lamch: case R not implemented.";
  } else if (template_blas_lsame(cmach, "M")) {
    std::cout << "ERROR in template_lapack_lamch: case M not implemented." << std::endl;
    throw "ERROR in template_lapack_lamch: case M not implemented.";
  } else if (template_blas_lsame(cmach, "U")) {
    std::cout << "ERROR in template_lapack_lamch: case U not implemented." << std::endl;
    throw "ERROR in template_lapack_lamch: case U not implemented.";
  } else if (template_blas_lsame(cmach, "L")) {
    std::cout << "ERROR in template_lapack_lamch: case L not implemented." << std::endl;
    throw "ERROR in template_lapack_lamch: case L not implemented.";
  } else if (template_blas_lsame(cmach, "O")) {
    std::cout << "ERROR in template_lapack_lamch: case O not implemented." << std::endl;
    throw "ERROR in template_lapack_lamch: case O not implemented.";
  }

  ret_val = rmach;
  return ret_val;

  /*     End of DLAMCH */

} /* dlamch_ */



#endif
