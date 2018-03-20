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
 

#ifndef TEMPLATE_LAPACK_LARNV_HEADER
#define TEMPLATE_LAPACK_LARNV_HEADER


template<class Treal>
int template_lapack_larnv(const integer *idist, integer *iseed, const integer *n, 
	Treal *x)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARNV returns a vector of n random real numbers from a uniform or   
    normal distribution.   

    Arguments   
    =========   

    IDIST   (input) INTEGER   
            Specifies the distribution of the random numbers:   
            = 1:  uniform (0,1)   
            = 2:  uniform (-1,1)   
            = 3:  normal (0,1)   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array   
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    N       (input) INTEGER   
            The number of random numbers to be generated.   

    X       (output) DOUBLE PRECISION array, dimension (N)   
            The generated random numbers.   

    Further Details   
    ===============   

    This routine calls the auxiliary routine DLARUV to generate random   
    real numbers from a uniform (0,1) distribution, in batches of up to   
    128 using vectorisable code. The Box-Muller method is used to   
    transform numbers from a uniform to a normal distribution.   

    =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer i__1, i__2, i__3;
    /* Local variables */
     integer i__;
     Treal u[128];
     integer il, iv;
     integer il2;

    --x;
    --iseed;

    /* Function Body */
    i__1 = *n;
    for (iv = 1; iv <= i__1; iv += 64) {
/* Computing MIN */
	i__2 = 64, i__3 = *n - iv + 1;
	il = minMACRO(i__2,i__3);
	if (*idist == 3) {
	    il2 = il << 1;
	} else {
	    il2 = il;
	}

/*        Call DLARUV to generate IL2 numbers from a uniform (0,1)   
          distribution (IL2 <= LV) */

	dlaruv_(&iseed[1], &il2, u);

	if (*idist == 1) {

/*           Copy generated numbers */

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1];
/* L10: */
	    }
	} else if (*idist == 2) {

/*           Convert generated numbers to uniform (-1,1) distribution */

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
/* L20: */
	    }
	} else if (*idist == 3) {

/*           Convert generated numbers to normal (0,1) distribution */

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = template_blas_sqrt(template_blas_log(u[(i__ << 1) - 2]) * -2.) * template_blas_cos(u[(
			i__ << 1) - 1] * 6.2831853071795864769252867663);
/* L30: */
	    }
	}
/* L40: */
    }
    return 0;

/*     End of DLARNV */

} /* dlarnv_ */

#endif
