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

/** @file aos.cc Blocked version of orbtial evaluation routines.
 Written by Pawel Salek.
 ===================================================================
 RETURNS: GSO: evaluated orbitals for a batch of grid points.
     GSO(:,:,1) contains orbital values.
     GSO(:,:,2:4) contains first geom. derivatives. - if requested.
     GSO(:,:,5:10) contains second derivatives - if requested.
*/

#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "aos.h"


/* get_bf_vals:
   get column(s) corresponding to just values of one shell at set of
   grid points.
*/
static void
get_bf_vals(int nvclen, ergo_real *gao,
	    const ergo_real (*pa)[3], const ergo_real *pa2,
            const BasisInfoStruct& bis,
	    const ShellSpecStruct& currShell)
{
  /*
      DIMENSION GAO(NVCLEN,KCKTA)
      DIMENSION PA(3,NVCLEN), PA2(NVCLEN)
      DIMENSION CSP(KHKTA,KCKTA)
      DIMENSION GA(NVCLEN),CINT(NVCLEN)
  */
  int i, j, k, nTerms, m;
  ergo_real sum, factor;

  ergo_real ga[DFT_MAX_BLLEN];
  memset(ga, 0, nvclen*sizeof(ergo_real));

  /* loop over primitives for selected shell to get the radial part */
  for(i=0; i<currShell.noOfContr; i++) {
    ergo_real fac = currShell.coeffList[i];
    ergo_real malpha = -currShell.exponentList[i];
    for(k=0; k < nvclen; k++)
      ga[k] += fac * template_blas_exp(malpha*pa2[k]);
  }

  /* multiply by the angular part */
  switch (currShell.shellType) {
  case 0: /* s function */
    memcpy(gao, ga, nvclen*sizeof(ga[1]));
    break;
  case 1:
    /* p function */
    for(k=0; k<nvclen; k++) {
      gao[k]            = pa[k][1]*ga[k];
      gao[k + nvclen]   = pa[k][2]*ga[k];
      gao[k + nvclen*2] = pa[k][0]*ga[k];
    }
    break;
  default:
    for(i = 0; i < currShell.noOfBasisFuncs; i++)
      {
	nTerms = bis.basisFuncList[currShell.startIndexInMatrix+i].
	  noOfTermsInPolynomial;
	for(k=0; k<nvclen; k++) {
	  //gao[k]            = pa[k][0]*ga[k];
	  sum = 0;
	  for(j = 0; j < nTerms; j++)
	    {
	      basis_func_term_struct* currTerm = 
		&bis.basisFuncList[currShell.startIndexInMatrix+i].poly[j];
	      factor = currTerm->coeff;
	      for(m = 0; m < currTerm->monomialInts[0]; m++)
		factor *= pa[k][0];
	      for(m = 0; m < currTerm->monomialInts[1]; m++)
		factor *= pa[k][1];
	      for(m = 0; m < currTerm->monomialInts[2]; m++)
		factor *= pa[k][2];
	      sum += factor;
	    } // END FOR j
	  gao[k+nvclen*i] = sum * ga[k];
	} // END FOR k
      } // END FOR i
  }
}

static void
get_bf_vals_derivs(int nvclen, ergo_real *gao,
                   const ergo_real (*pa)[3], const ergo_real *pa2,
                   const BasisInfoStruct& bis,
                   const ShellSpecStruct& currShell)
{
  /*
      DIMENSION GAO(NVCLEN,KCKTA)
      DIMENSION PA(3,NVCLEN), PA2(NVCLEN)
      DIMENSION CSP(KHKTA,KCKTA)
      DIMENSION GA(NVCLEN),CINT(NVCLEN)
  */
  int i, j, k, nTerms, m, nbast=bis.noOfBasisFuncs;

  ergo_real ga[DFT_MAX_BLLEN], gu[DFT_MAX_BLLEN];
  memset(ga, 0, nvclen*sizeof(ergo_real));
  memset(gu, 0, nvclen*sizeof(ergo_real));

  /* loop over primitives for selected shell to get the radial part */
  for(i=0; i<currShell.noOfContr; i++) {
      for(k=0; k < nvclen; k++) {
          ergo_real fac = currShell.coeffList[i] * 
            template_blas_exp(-currShell.exponentList[i]*pa2[k]);
          ga[k] += fac;
          gu[k] -= 2*currShell.exponentList[i]*fac;
      }
  }
  
  /* multiply by the angular part */
  switch (currShell.shellType) {
  case 0: /* s function */
    memcpy(gao, ga, nvclen*sizeof(ergo_real));
    for(k=0; k < nvclen; k++) {
        gao[k + nvclen*nbast]   = pa[k][0]*gu[k]; /* d/dx */
        gao[k + nvclen*nbast*2] = pa[k][1]*gu[k]; /* d/dy */
        gao[k + nvclen*nbast*3] = pa[k][2]*gu[k]; /* d/dz */
    }
    break;
  default: /* p and higher functions */
    for(i = 0; i < currShell.noOfBasisFuncs; i++)
      {
	nTerms = bis.basisFuncList[currShell.startIndexInMatrix+i].
	  noOfTermsInPolynomial;
	for(k=0; k<nvclen; k++) {
          ergo_real x = pa[k][0], y = pa[k][1], z = pa[k][2];
          ergo_real factor, tx, ty, tz, fr = 0, fx = 0, fy = 0, fz = 0;
	  for(j = 0; j < nTerms; j++)
	    {
	      basis_func_term_struct* currTerm = 
		&bis.basisFuncList[currShell.startIndexInMatrix+i].poly[j];
              int ix = currTerm->monomialInts[0];
              int iy = currTerm->monomialInts[1];
              int iz = currTerm->monomialInts[2];
	      factor = currTerm->coeff;
	      for(m = ix-2; m >=0; m--) factor *= x;
	      for(m = iy-2; m >=0; m--) factor *= y;
	      for(m = iz-2; m >=0; m--) factor *= z;
              tx = currTerm->monomialInts[0]>0 ? x : 1;
              ty = currTerm->monomialInts[1]>0 ? y : 1;
              tz = currTerm->monomialInts[2]>0 ? z : 1;
	      fr += factor*tx*ty*tz;
              fx += factor*ix*ty*tz;
              fy += factor*iy*tx*tz;
              fz += factor*iz*tx*ty;
	    } // END FOR j
	  gao[k+nvclen*i]           = fr * ga[k]; /* value */
	  gao[k+nvclen*(i+nbast)]   = fx * ga[k] + fr*gu[k]*x; /* d/dx */
	  gao[k+nvclen*(i+nbast*2)] = fy * ga[k] + fr*gu[k]*y; /* d/dy */
	  gao[k+nvclen*(i+nbast*3)] = fz * ga[k] + fr*gu[k]*z; /* d/dz */
	} // END FOR k
      } // END FOR i
  }
}


void
dft_get_orbs(int nvclen,
             ergo_real *gao,
             const ergo_real (*coor)[3],
             int nblcnt, int (*iblcks)[2],
             int nder,
             const BasisInfoStruct& bis)
{
  /*
      DIMENSION GAO(NVCLEN,NBAST,NTYPSO),WORK(LWORK)
      DIMENSION COOR(3,NVCLEN)
      DIMENSION IBLCKS(2,NBLCNT)
C     PA2 contains distance from the basis function center to respective
c     grid point.
      DIMENSION PA(3,NVCLEN), PA2(NVCLEN)
  */
  int ibl, ishela, i; 
  ergo_real pa[DFT_MAX_BLLEN][3], pa2[DFT_MAX_BLLEN];
  void (*eval)(int nvclen, ergo_real *gao,
               const ergo_real (*pa)[3], const ergo_real *pa2,
               const BasisInfoStruct& bis, const ShellSpecStruct& currShell);

  switch(nder) {
  case 0: eval = get_bf_vals;        break; /* orbital values only */
  case 1: eval = get_bf_vals_derivs; break; 
  default: abort(); /* not reached */
  }

  for(ibl=0; ibl<nblcnt; ibl++) {
    for(ishela = iblcks[ibl][0]; ishela<iblcks[ibl][1]; ishela++) {
      ShellSpecStruct& currShell = bis.shellList[ishela];
      ergo_real cenx = currShell.centerCoords[0];
      ergo_real ceny = currShell.centerCoords[1];
      ergo_real cenz = currShell.centerCoords[2];
      for(i=0; i<nvclen; i++) {
	pa[i][0] = coor[i][0]-cenx;
	pa[i][1] = coor[i][1]-ceny;
	pa[i][2] = coor[i][2]-cenz;
	pa2[i] = pa[i][0]*pa[i][0] + pa[i][1]*pa[i][1] + pa[i][2]*pa[i][2];
      }
      eval(nvclen, &gao[currShell.startIndexInMatrix*nvclen],
           pa, pa2, bis, currShell);
      
    }
  }
  /* done! */
}
