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

/** @file integrals_2el_explicit.cc

    @brief Code for explicit computation of 4-index 2-electron
    integrals.

    @author: Elias Rudberg <em>responsible</em>
*/

/* Written by Elias Rudberg, KTH, Stockholm */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include "integrals_2el_explicit.h"
#include "memorymanag.h"
#include "pi.h"
#include "output.h"
#include "utilities.h"
#include "boysfunction.h"
#include "integral_info.h"
#include "integrals_general.h"
#include "integrals_2el_single.h"



typedef struct
{
  int a, b, c, d;
  int poly_ab_index;
  int poly_cd_index;
} abcd_struct;

#define set_abcd_list_item_macro(i,A,B,C,D) \
list[i].a = A; list[i].b = B; list[i].c = C; list[i].d = D;




static int globalCount = 0;


ergo_real
do_2e_integral(int mu, 
	       int nu, 
	       int la, 
	       int si, 
	       const BasisInfoStruct & basisInfo, 
	       const IntegralInfo & integralInfo) {
  return do_2e_integral_general(mu, 
				nu, 
				la, 
				si, 
				basisInfo, 
				basisInfo, 
				basisInfo, 
				basisInfo, 
				integralInfo);
}

ergo_real 
do_2e_integral_general(int mu, 
		       int nu, 
		       int la, 
		       int si, 
		       const BasisInfoStruct & basisInfo_mu, 
		       const BasisInfoStruct & basisInfo_nu, 
		       const BasisInfoStruct & basisInfo_la, 
		       const BasisInfoStruct & basisInfo_si, 
		       const IntegralInfo & integralInfo)
{
  int n_psi2, i, j;
  ergo_real sum, currIntegral;
  const int maxCount = 1000;
  DistributionSpecStruct list_psi1[maxCount];
  DistributionSpecStruct list_psi2[maxCount];

  /* form product of basisfuncs mu and nu, store product in psi1 */
  int n_psi1 = get_product_simple_primitives(basisInfo_mu, mu,
					     basisInfo_nu, nu,
					     list_psi1,
					     maxCount,
					     0);
  if(n_psi1 <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives\n");
      exit(EXIT_FAILURE);
      return 0;
    }

  /* form product of basisfuncs la and si, store product in psi2 */
  n_psi2 = get_product_simple_primitives(basisInfo_la, la,
					 basisInfo_si, si,
					 list_psi2,
					 maxCount,
					 0);
  if(n_psi2 <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives\n");
      exit(EXIT_FAILURE);
      return 0;
    }

  const JK::ExchWeights CAM_params_not_used;

  sum = 0;
  for(i = 0; i < n_psi1; i++)
    {
      DistributionSpecStruct* prim_psi1 = &list_psi1[i];
      for(j = 0; j < n_psi2; j++)
	{
	  DistributionSpecStruct* prim_psi2 = &list_psi2[j];
	  globalCount++;
	  currIntegral = do_2e_integral_using_symb_info(CAM_params_not_used, prim_psi1, prim_psi2, integralInfo);
	  sum += currIntegral;
	} /* END FOR j */
    } /* END FOR i */
  
  return sum;  
}


/** compute_2e_matrix_simple computes the 2el matrix in the simplest
    possible way. It assumes that the matrix is computed for closed
    shell. The weight of the HF exchange is controlled by @param hf_weight
    which is equal 1 for ordinary Hartree-Fock calculation. No
    assumption are made regarding symmetry of the density matrix
    @param dens . The computed two-electron part of the Fock matrix is
    returned in @param result .
    @param basisInfo info about the used basis set.
    @param integralInfo info needed for evaluation of integrals of Gaussian functions.
*/
int 
compute_2e_matrix_simple(const BasisInfoStruct & basisInfo,
			 const IntegralInfo & integralInfo,
                         ergo_real hf_weight,
			 ergo_real* result,
			 const ergo_real* dens)
{
  int mu, nu, sigma, lambda;
  int nbast = basisInfo.noOfBasisFuncs;
  ergo_real munusila, mulasinu, sum;
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "entering compute_2e_matrix HF_WEIGHT=%f", (double)hf_weight);
  
  for(mu = 0; mu < nbast; mu++)
    {
      for(nu = 0; nu < nbast; nu++)
	{
	  sum = 0;
	  for(lambda = 0; lambda < nbast; lambda++)
	    {
	      for(sigma = 0; sigma < nbast; sigma++)
		{
		  munusila = do_2e_integral(mu, nu, sigma, lambda, basisInfo, integralInfo);
		  mulasinu = do_2e_integral(mu, lambda, sigma, nu, basisInfo, integralInfo);
		  sum += 
		    dens[lambda*nbast+sigma] * 
		    (munusila - 0.5 * hf_weight * mulasinu);
		} /* END FOR sigma */
	    } /* END FOR lambda */
	  result[mu*nbast+nu] = sum;
	} /* END FOR nu */
    } /* END FOR mu */

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_2e_matrix ending OK\n");
  return 0;
}




static int 
compute_J_and_K_integraldriven(const BasisInfoStruct & basisInfo,
			       const IntegralInfo & integralInfo,
			       ergo_real* J,
			       ergo_real* K,
			       ergo_real* dens)
{
  int n, nBytes, i, j, count, a, b, c, d;

  Util::TimeMeter timeMeter;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "entering compute_J_and_K_integraldriven");
  n = basisInfo.noOfBasisFuncs;
  nBytes = n * n * sizeof(ergo_real);
  memset(J, 0, nBytes);
  memset(K, 0, nBytes);

  count = 0;
  a = 0;
  b = 0;
  c = 0;
  d = 0;
  while(d < n)
    {
      abcd_struct list[8];

      /* compute integral */
      ergo_real integralValue = do_2e_integral(a, b, c, d, basisInfo, integralInfo);

      count++;
      
      /* determine unique configurations */
      set_abcd_list_item_macro(0, a, b, c, d);
      set_abcd_list_item_macro(1, a, b, d, c);
      set_abcd_list_item_macro(2, b, a, c, d);
      set_abcd_list_item_macro(3, b, a, d, c);
      set_abcd_list_item_macro(4, c, d, a, b);
      set_abcd_list_item_macro(5, d, c, a, b);
      set_abcd_list_item_macro(6, c, d, b, a);
      set_abcd_list_item_macro(7, d, c, b, a);

      for(i = 0; i < 8; i++)
	{
	  abcd_struct* abcd = &list[i];
	  int aa, bb, cc, dd;

	  /* check if this is a new unique configuration */
	  int unique = 1;
	  for(j = 0; j < i; j++)
	    {
	      if(abcd->a == list[j].a && 
		 abcd->b == list[j].b && 
		 abcd->c == list[j].c && 
		 abcd->d == list[j].d)
		unique = 0;
	    }
	  if(unique == 0)
	    continue;
	  /* now we know that this configuration is unique. */
	  aa = abcd->a;
	  bb = abcd->b;
	  cc = abcd->c;
	  dd = abcd->d;

#if 1
	  /* add contribution to coulomb matrix */
	  J[aa*n+bb] += dens[cc*n+dd] * integralValue;

	  /* add contribution to exchange matrix */
	  K[aa*n+dd] += -0.5 * dens[bb*n+cc] * integralValue;
#endif
	} /* END FOR i go through 8 configurations */
      
      /* now get numbers for next unique integral */
      d++;
      if(d < n)
	continue;

      /* d has hit the roof */
      c++;
      if(c < n)
	{
	  d = c;
	  continue;
	}

      /* c has hit roof */
      b++;
      if(b < n)
	{
	  c = a;
	  d = b;
	  continue;
	}

      /* b has hit roof */
      a++;
      if(a < n)
	{
	  b = a;
	  c = a;
	  d = a;
	  continue;
	}

      /* a has hit roof. This means that we are done. */
      break;

    } /* END WHILE more unique integrals */
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_J_and_K_integraldriven ending OK");
  timeMeter.print(LOG_AREA_INTEGRALS, "compute_J_and_K_integraldriven");
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "number of unique integrals computed: %i", count);
  
  return 0;
}



int
compute_2e_matrix_list_explicit(const BasisInfoStruct & basisInfo,
				const IntegralInfo & integralInfo,
				ergo_real** resultList,
				ergo_real** densList,
				int noOfMatrices,
				ergo_real threshold)
{
  ergo_real* J;
  int n, i, j;

  if(noOfMatrices != 1)
    do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_2e_matrix_list_explicit: (noOfMatrices != 1), will take some time");
  
  n = basisInfo.noOfBasisFuncs;
  J = (ergo_real*)ergo_malloc(n*n*sizeof(ergo_real));
  for(j = 0; j < noOfMatrices; j++)
    {
      if(compute_J_and_K_integraldriven(basisInfo,
					integralInfo,
					J,
					resultList[j],
					densList[j]) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_and_K_integraldriven");
	  return -1;
	}
      for(i = 0; i < n*n; i++)
	{
	  resultList[j][i] += J[i];
	} // END FOR i
    } // END FOR j
  ergo_free(J);
  return 0;
}


