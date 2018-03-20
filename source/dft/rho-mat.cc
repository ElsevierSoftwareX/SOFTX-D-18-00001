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

/** @file rho-mat.cc Functions for density and gradient evaluation.
    The density can be evaluated at entire batches of grid points.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#include <stdio.h>

#include "realtype.h"
typedef ergo_real real;

#include "mat_gblas.h"
#include "rho-mat.h"

#if !defined(restrict)
#define restrict
#endif

/** helper function for zeroing only used blocks of orbitals. Selected
    values for the first index are looped over, and all allowed values
    of the second index.

    @param tmp the matrix[nbast][nvclen]

    @param nblocks pointer to an integer containing number of nonzero
    orbital blocks.
    @param iblocks a set of nblocks integer pairs [a,b) defining the
    range of first index to be zeroed.
    @param ldaib not used
    @param nvclen batch length - and the second dimension of tmp.
 */
static void
zeroorbs(real *tmp, const int *nblocks, const int (*iblocks)[2],
         int ldaib, int nvclen)
{
    /* DIMENSION TMP(NVCLEN,NBAST),NBLOCKS(NSYM),IBLOCKS(2,LDAIB,NSYM) */
    int ibl, idx, k; 
    for(ibl=0; ibl<nblocks[0]; ibl++)
        for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
            real * tmpi = tmp + idx*nvclen;
            for(k=0; k<nvclen; k++) tmpi[k] = 0.0;
        }
}

/** Computes the expectation value <o|dmat|o'> for a symmetric matrix
  and given set of precomputed orbital values gao. Sparsity of gao as
  determined with help of nblocks and iblocks is used to reduce the
  computational effort.
  @param dmat full square symmetric matrix
  @param nbast size of dmat
  @param gao orbital matrix[nbast][nvclen]. Set values are determied
  by nblocks and iblocks, other values shall not be accessed.
  @param nblocks number of nonzero row blocks in gao
  @param iblocks ranges [a,b) of nonzero blocks in gao.
  @param ldaib not used
  @param tmp temporary matrix [nbast][nvclen]
  @param nvclen batch length - number of columns in gao.
  @param rho the vector[nvclen] where the computed expectation values
   will be stored.
*/
void
getrho_blocked_lda(int nbast, const real * dmat, const real * restrict gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen, real *rho)
{
/*
      DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
      DIMENSION NBLOCKS(NSYM), IBLOCKS(2,LDAIB,NSYM), RHO(NVCLEN)
      DIMENSION TMP(NVCLEN,NBAST)
*/
    int ibl, idx, jbl, k;

    /* dzero(TMP,NVCLEN*NBAST) */
    zeroorbs(tmp, nblocks, iblocks, ldaib, nvclen);

    /* only first symmetry */
#if USE_BLAS_IN_XC
    for(ibl=0; ibl<nblocks[0]; ibl++) {
      static const ergo_real ONER = 1.0;
      static const ergo_real HALF = 0.5;
      int cColumns, sumCols; 
      for(jbl=0; jbl<ibl; jbl++) {
	cColumns = iblocks[jbl][1]-iblocks[jbl][0];
	sumCols = iblocks[ibl][1]-iblocks[ibl][0];
	mat::gemm("N", "N", &nvclen, &cColumns, &sumCols, &ONER,
		  gao+iblocks[ibl][0]*nvclen, &nvclen,
		  dmat+iblocks[jbl][0]*nbast + iblocks[ibl][0], &nbast, &ONER,
		  tmp + iblocks[jbl][0]*nvclen, &nvclen);
      }	
      cColumns = iblocks[ibl][1]-iblocks[ibl][0];
      sumCols = iblocks[ibl][1]-iblocks[ibl][0];
      mat::symm("R", "U", &nvclen, &cColumns, &HALF,
		dmat+iblocks[ibl][0]*nbast + iblocks[ibl][0], &nbast,
		gao+iblocks[ibl][0]*nvclen, &nvclen, &ONER,
		tmp + iblocks[ibl][0]*nvclen, &nvclen);
    }
#else /* USE_BLAS_IN_XC */
    int jdx;
    real d;
    for(ibl=0; ibl<nblocks[0]; ibl++)
        /* print *,"block ", IBL,IBLOCKS(1,IBL,ISYM),IBLOCKS(2,IBL,ISYM) */
        for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
            real * restrict tmpj;
            const real * restrict gaoi = gao + idx*nvclen;
            for(jbl=0; jbl<nblocks[0]; jbl++) {
                int jtop = iblocks[jbl][1] > idx ? idx : iblocks[jbl][1];
                for(jdx=iblocks[jbl][0]; jdx<jtop; jdx++) {
                  d = dmat[idx +jdx*nbast];
                    tmpj = tmp + jdx*nvclen;
                    for(k=0; k<nvclen; k++)
                        tmpj[k] += gaoi[k]*d;
                }
            }
            tmpj = tmp + idx*nvclen;
            d = dmat[idx +idx*nbast]*0.5;
            for(k=0; k<nvclen; k++)
                tmpj[k] += gaoi[k]*d;
        }
#endif /* USE_BLAS_IN_XC */
    memset(rho, 0, nvclen*sizeof(real));
    for(ibl=0; ibl<nblocks[0]; ibl++)
        for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
            const real * restrict gaoi = gao + idx*nvclen;
            const real * restrict tmpi = tmp + idx*nvclen;
            for(k=0; k<nvclen; k++)
                rho[k] += gaoi[k]*tmpi[k]*2.0;
        }
}

/** Computes the expectation value <o|dmat|o'> and its derivatives for
  a symmetric matrix and given set of precomputed orbital values and
  their cartesian derivatives gao. Sparsity of gao as determined with
  help of nblocks and iblocks is used to reduce the computational
  effort.

  @param dmat full square symmetric matrix
  @param nbast size of dmat

  @param gao orbital matrix[4][nbast][nvclen]. First block [0][][]
  contains orbital values. Subsequent blocks - orbital derivatives wrt
  x,y, and z coordinates. Set values are determied by nblocks and
  iblocks, other values shall not be accessed.

  @param nblocks number of nonzero row blocks in gao
  @param iblocks ranges [a,b) of nonzero blocks in gao.
  @param ldaib not used
  @param tmp temporary matrix [nbast][nvclen]
  @param nvclen batch length - number of columns in gao.
  @param rho the vector[nvclen] where the computed expectation values
  will be stored.
  @param grad a vector of triples where the computed gradient values
   will be stored.  
*/
void
getrho_blocked_gga(int nbast, const real * dmat, const real * restrict gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen,
                   real *rho, real (*grad)[3])
{
/*
      DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
      DIMENSION NBLOCKS(NSYM), IBLOCKS(2,LDAIB,NSYM), RHO(NVCLEN)
      DIMENSION TMP(NVCLEN,NBAST)
*/
    int ibl, idx, jbl, k;

    /* dzero(TMP,NVCLEN*NBAST) */
    zeroorbs(tmp, nblocks, iblocks, ldaib, nvclen);
    /* only first symmetry */
#if USE_BLAS_IN_XC == 1
    for(ibl=0; ibl<nblocks[0]; ibl++)
      for(jbl=0; jbl<nblocks[0]; jbl++) {
	static const ergo_real ONER = 1.0;
	int cColumns = iblocks[jbl][1]-iblocks[jbl][0];
	int sumCols = iblocks[ibl][1]-iblocks[ibl][0];
	mat::gemm("N", "N", &nvclen, &cColumns, &sumCols, &ONER,
		  gao+iblocks[ibl][0]*nvclen, &nvclen,
		  dmat+iblocks[jbl][0]*nbast + iblocks[ibl][0],
		  &nbast, &ONER,
		  tmp + iblocks[jbl][0]*nvclen, &nvclen);
      }
#else /* USE_BLAS_IN_XC */
    int jdx;
    real d;
    for(ibl=0; ibl<nblocks[0]; ibl++)
        /* print *,"block ", IBL,IBLOCKS(1,IBL,ISYM),IBLOCKS(2,IBL,ISYM) */
        for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
            const real * restrict gaoi = gao + idx*nvclen;
            for(jbl=0; jbl<nblocks[0]; jbl++) {
                for(jdx=iblocks[jbl][0]; jdx<iblocks[jbl][1]; jdx++) {
                    real * restrict tmpj = tmp + jdx*nvclen;
                    d = dmat[idx +jdx*nbast];
                    for(k=0; k<nvclen; k++)
                        tmpj[k] += gaoi[k]*d;
                }
            }
        }
#endif /* USE_BLAS_IN_XC */

    memset(rho,  0,   nvclen*sizeof(real));
    memset(grad, 0, 3*nvclen*sizeof(real));
    for(ibl=0; ibl<nblocks[0]; ibl++)
        for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
            const real * restrict gaoi = gao + idx*nvclen;
            const real * restrict tmpi = tmp + idx*nvclen;
            for(k=0; k<nvclen; k++) {
                rho[k]     += gaoi[k]*tmpi[k];
                grad[k][0] += gaoi[k + nvclen*nbast]*tmpi[k]*2;
                grad[k][1] += gaoi[k + nvclen*nbast*2]*tmpi[k]*2;
                grad[k][2] += gaoi[k + nvclen*nbast*3]*tmpi[k]*2;
            }
        }
}

/** Computes the expectation value <o|dmat|o'> for a nonsymmetric matrix
  and given set of precomputed orbital values gao. Sparsity of gao as
  determined with help of nblocks and iblocks is used to reduce the
  computational effort.
  @param dmat full square symmetric matrix
  @param nbast size of dmat
  @param gao orbital matrix[nbast][nvclen]. Set values are determied
  by nblocks and iblocks, other values shall not be accessed.
  @param nblocks number of nonzero row blocks in gao
  @param iblocks ranges [a,b) of nonzero blocks in gao.
  @param ldaib not used
  @param tmp temporary matrix [nbast][nvclen]
  @param nvclen batch length - number of columns in gao.
  @param rho the vector[nvclen] where the computed expectation values
   will be stored.
*/
void
getexp_blocked_lda(int nbast, const real * dmat, const real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen, real *rho)
{
  /*
    DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
    DIMENSION NBLOCKS(NSYM), IBLOCKS(2,LDAIB,NSYM), RHO(NVCLEN)
    DIMENSION TMP(NVCLEN,NBAST)
  */
  int ibl, idx, jbl, jdx, k;
  real d;

  /* dzero(TMP,NVCLEN*NBAST) */
  zeroorbs(tmp, nblocks, iblocks, ldaib, nvclen);

  /* only first symmetry */
  for(ibl=0; ibl<nblocks[0]; ibl++)
    /* print *,"block ", IBL,IBLOCKS(1,IBL,ISYM),IBLOCKS(2,IBL,ISYM) */
    for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
      const real * gaoi = gao + idx*nvclen;
      for(jbl=0; jbl<nblocks[0]; jbl++) {
        for(jdx=iblocks[jbl][0]; jdx<iblocks[jbl][1]; jdx++) {
          real *tmpj = tmp + jdx*nvclen;
          d = dmat[idx +jdx*nbast];
          for(k=0; k<nvclen; k++)
            tmpj[k] += gaoi[k]*d;
        }
      }
    }

  memset(rho, 0, nvclen*sizeof(real));
  for(ibl=0; ibl<nblocks[0]; ibl++)
    for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
      const real * gaoi = gao + idx*nvclen;
      const real * tmpi = tmp + idx*nvclen;
      for(k=0; k<nvclen; k++)
        rho[k] += gaoi[k]*tmpi[k];
    }
}

/** Computes the expectation value <o|dmat|o'> and its derivatives for
  a nonsymmetric matrix and given set of precomputed orbital values
  and their cartesian derivatives gao. Sparsity of gao as determined
  with help of nblocks and iblocks is used to reduce the computational
  effort.

  @param dmat full square symmetric matrix
  @param nbast size of dmat

  @param gao orbital matrix[4][nbast][nvclen]. First block [0][][]
  contains orbital values. Subsequent blocks - orbital derivatives wrt
  x,y, and z coordinates. Set values are determied by nblocks and
  iblocks, other values shall not be accessed.

  @param nblocks number of nonzero row blocks in gao
  @param iblocks ranges [a,b) of nonzero blocks in gao.
  @param ldaib not used
  @param tmp temporary matrix [nbast][nvclen]
  @param nvclen batch length - number of columns in gao.
  @param rgrad a vector of quartets where the computed expectation
  values and gradient values will be stored.  */
void
getexp_blocked_gga(int nbast, const real * dmat, const real * gao,
                   const int* nblocks, const int (*iblocks)[2],
                   int ldaib, real *tmp, int nvclen,
                   real (*rgrad)[4])
{
/*
      DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
      DIMENSION NBLOCKS(NSYM), IBLOCKS(2,LDAIB,NSYM), RHO(NVCLEN)
      DIMENSION TMP(NVCLEN,NBAST)
*/
    int ibl, idx, jbl, jdx, k;
    real d;

    /* dzero(TMP,NVCLEN*NBAST) */
    zeroorbs(tmp, nblocks, iblocks, ldaib, nvclen);

    /* only first symmetry */
    for(ibl=0; ibl<nblocks[0]; ibl++)
        /* print *,"block ", IBL,IBLOCKS(1,IBL,ISYM),IBLOCKS(2,IBL,ISYM) */
        for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
            const real * gaoi = gao + idx*nvclen;
            for(jbl=0; jbl<nblocks[0]; jbl++) {
                for(jdx=iblocks[jbl][0]; jdx<iblocks[jbl][1]; jdx++) {
                    real *tmpj = tmp + jdx*nvclen;
                    d = dmat[idx +jdx*nbast] + dmat[jdx +idx*nbast];
                    for(k=0; k<nvclen; k++)
                        tmpj[k] += gaoi[k]*d;
                }
            }
        }

    memset(rgrad, 0, 4*nvclen*sizeof(real));
    for(ibl=0; ibl<nblocks[0]; ibl++)
      for(idx=iblocks[ibl][0]; idx<iblocks[ibl][1]; idx++) {
        const real * gaoi = gao + idx*nvclen;
        const real * tmpi = tmp + idx*nvclen;
        for(k=0; k<nvclen; k++) {
          rgrad[k][0] += gaoi[k                 ]*tmpi[k]*0.5;
          rgrad[k][1] += gaoi[k + nvclen*nbast  ]*tmpi[k];
          rgrad[k][2] += gaoi[k + nvclen*nbast*2]*tmpi[k];
          rgrad[k][3] += gaoi[k + nvclen*nbast*3]*tmpi[k];
        }
      }
}
