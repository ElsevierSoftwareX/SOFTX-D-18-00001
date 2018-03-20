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

#ifndef RHOMAT_HEADER
#define RHOMAT_HEADER 1

#include "grid_matrix.h"

/** @file rho-mat.h Density and gradient evaluation interface. */

void
getrho_blocked_lda(int nbast, const real * dmat, const real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen, real *rho);

inline void
getrho_blocked_lda(int nbast, const Dft::FullMatrix& m, const real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen, real *rho)
{
  getrho_blocked_lda(nbast, m.mat, gao, nblocks, iblocks,
                     ldaib, tmp, nvclen, rho);
}


void
getrho_blocked_gga(int nbast, const real * dmat, const real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen,
                   real *rho, real (*grad)[3]);

inline void
getrho_blocked_gga(int nbast, const Dft::FullMatrix& dmat, const real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen,
                   real *rho, real (*grad)[3])
{
  getrho_blocked_gga(nbast, dmat.mat, gao, nblocks, iblocks,
                     ldaib, tmp, nvclen, rho, grad);
}

void
getexp_blocked_lda(int nbast, const real * dmat, const real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen, real *rho);
void
getexp_blocked_gga(int nbast, const real * dmat, const real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen,
                   real (*rgrad)[4]);

#endif
