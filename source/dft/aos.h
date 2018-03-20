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

#if !defined(_AOS_H_)
#define _AOS_H_

/** @file aos.h Blocked version of orbtial evaluation routines. */
#include "realtype.h"
#include "basisinfo.h"

/** Limit for the number of grid point batch length. Should not be too
 * short because the loop overhead will grow too large, nor too long
 * because we run out of cache then. */

#define DFT_MAX_BLLEN 192

/** Computes values of basis functions at specified points in
    space. Only b.fs specified by iblcks[nblcnt] are computed.

    @param nvclen number of points to consider. IT must be smaller
    than DFT_MAX_BLLEN.

    @param gao matrix of computed b.f, values. Dimension is in C
    convention: [nderivatives][nvclen], where nderivatives is 1 for
    nder==0, 4 for nder==1.

    @param coor point coordinates.

    @param nblcnt so many continous blocks of basis functions will be computed.

    @param iblcks start and end indices of the b.fs shells. Computed shells
    are [iblcks[0], iblcks[1]).

    @param nder whether orbital derivatives are to be computed as
    well. Allowed values are 0 (no derivatives) and 1 (values and
    first order derivatives).

    @param bis structure describing the basis functions to be evaluated.

*/
void dft_get_orbs(int nvclen, ergo_real *gao,
		  const ergo_real (*coor)[3],
		  int nblcnt, int (*iblcks)[2],
		  int nder,
		  const BasisInfoStruct& bis);

#endif /* _AOS_H_ */
