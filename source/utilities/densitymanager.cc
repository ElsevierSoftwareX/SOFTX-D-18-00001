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

/** @file densitymanager.cc

    @brief Functionality for working with the electron density as a
    function of space, for a given basis set and density matrix.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "memorymanag.h"
#include "output.h"
#include "densitymanager.h"
#include "pi.h"
#include "integrals_general.h"
#include "template_blas_common.h"
 

#define EXPONENT_DIFF_LIMIT 1e-22
#define DISTR_CENTER_DIST_LIMIT 1e-22


static ergo_real 
compute_1d_gaussian_integral_recursive(ergo_real a, ergo_real b, int n, ergo_real alpha)
{
  ergo_real result, sqrtalpha, term1, term2;
  ergo_real aToPowerNminus1, bToPowerNminus1;
  if(n == 0)
    {
      sqrtalpha = template_blas_sqrt(alpha);
      result = template_blas_sqrt(pi/(4*alpha)) * (template_blas_erf(sqrtalpha*b) - template_blas_erf(sqrtalpha*a));
      return result;
    }
  if(n == 1)
    {
      result = -(1 / (2*alpha)) * (template_blas_exp(-alpha*b*b) - template_blas_exp(-alpha*a*a));
      return result;
    }
  if(n < 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in 1dintegral: n < 0");
      exit(0);
    }
  /* now we know that n >= 2 */
  term1 = (n - 1) * compute_1d_gaussian_integral_recursive(a, b, n-2, alpha);
  aToPowerNminus1 = template_blas_pow(a, (ergo_real)n-1);
  bToPowerNminus1 = template_blas_pow(b, (ergo_real)n-1);
  term2  = 
    bToPowerNminus1 * template_blas_exp(-alpha*b*b) - 
    aToPowerNminus1 * template_blas_exp(-alpha*a*a);
  result = (term1 - term2) / (2 * alpha);
  /*  return 0; */
  return result;
} /* END compute_1d_gaussian_integral_recursive */



static ergo_real 
compute_integral_over_box(DistributionSpecStruct* distr, 
			  ergo_real* minVect, ergo_real* maxVect,
			  std::vector<int> monomialIntsAdd = std::vector<int>(3, 0))
{
  ergo_real result, a, b, alpha;
  int i, n;

  result = distr->coeff;
  alpha = distr->exponent;
  for(i = 0; i < 3; i++) // run over the coordinates
    {
      n = distr->monomialInts[i];
      a = minVect[i] - distr->centerCoords[i];
      b = maxVect[i] - distr->centerCoords[i];
      result *= compute_1d_gaussian_integral_recursive(a, b, n, alpha); 
    } /* END FOR i */
  return result;
} /* END compute_integral_over_box */







ergo_real
integrate_density_in_box(int nPrims,
			 DistributionSpecStruct* rho,
			 ergo_real mid_x,
			 ergo_real mid_y,
			 ergo_real mid_z,
			 ergo_real box_width)
{
  ergo_real minVect[3];
  ergo_real maxVect[3];
  minVect[0] = mid_x - 0.5 * box_width;
  maxVect[0] = mid_x + 0.5 * box_width;
  minVect[1] = mid_y - 0.5 * box_width;
  maxVect[1] = mid_y + 0.5 * box_width;
  minVect[2] = mid_z - 0.5 * box_width;
  maxVect[2] = mid_z + 0.5 * box_width;
  ergo_real sum = 0;
  int i;
  for(i = 0; i < nPrims; i++)
    sum += compute_integral_over_box(&rho[i], minVect, maxVect);
  return sum;
}



ergo_real
integrate_density_in_box_2(int nPrims,
			   DistributionSpecStruct* rho,
			   ergo_real* minVect, 
			   ergo_real* maxVect,
			   std::vector<int> monomialIntsAdd)
{
  ergo_real sum = 0;
  int i;
  for(i = 0; i < nPrims; i++)
    sum += compute_integral_over_box(&rho[i], minVect, maxVect, monomialIntsAdd);
  return sum;
}




int
get_no_of_primitives_for_density(ergo_real cutoff,
				 const ergo_real *dmat,
				 const BasisInfoStruct & basisInfo)
{
#define MAX_DISTR_IN_TEMP_LIST 888

  int i, j;
  int symmetryFactor;
  int nBasisFuncs, nn;
  
  nBasisFuncs = basisInfo.noOfBasisFuncs;
  nn = 0;
  for(i = 0; i < nBasisFuncs; i++)
    {
      for(j = 0; j < nBasisFuncs; j++)
	{
	  DistributionSpecStruct tempList[MAX_DISTR_IN_TEMP_LIST];
	  int nPrimitives, k;
	  /* the matrix M is symmetric: include diagonal terms once, */
	  /* and include upper off-diagonal terms multiplied by 2 */
	  if(i == j)
              symmetryFactor = 1;
	  else
	    symmetryFactor = 2;
	  if(i > j)
	    continue;
          nPrimitives = 
	    get_product_simple_primitives(basisInfo, i,
					  basisInfo, j,
					  tempList,
					  MAX_DISTR_IN_TEMP_LIST,
					  0);
	  if(nPrimitives <= 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in get_product_simple_primitives");
	      return -1;
	    }
	  for(k = 0; k < nPrimitives; k++)
	    {
	      DistributionSpecStruct* currDistr = &tempList[k];
	      ergo_real Mij = dmat[i*nBasisFuncs+j];
	      ergo_real newCoeff = currDistr->coeff * Mij * symmetryFactor;
	      if(template_blas_fabs(newCoeff) > cutoff)
		nn++;
	    }
	}
    }
  return nn;
}










static int 
do_merge_sort_distrs(int n, 
		     DistributionSpecStruct* list, 
		     DistributionSpecStruct* workList)
{
    /* merge sort:  */
    /* first sort the first half, */
    /* then sort the second half, */
    /* then merge results to form final sorted list. */
  int n1, n2, nn, decision, i1, i2, i;
  DistributionSpecStruct* d1;
  DistributionSpecStruct* d2;

  if(n == 0)
    return 0;

  if(n < 1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "(n < 1)");
      return -1;
    }
  if(n == 1)
    return 0;
  
  n1 = n / 2;
  n2 = n - n1;

  /* sort first half */
  if(do_merge_sort_distrs(n1, list, workList) != 0)
    return -1;

  /* sort second half */
  if(do_merge_sort_distrs(n2, &list[n1], workList) != 0)
    return -1;

  /* merge results */
  nn = 0;
  i1 = 0;
  i2 = 0;
  while(nn < n)
    {
      if((i1 < n1) && (i2 < n2))
	{
            /* compare */
	  d1 = &list[i1];
	  d2 = &list[n1+i2];
	  decision = 0;
	  for(i = 0; i < 3; i++)
	    {
	      if(decision == 0)
		{
		  if(d1->monomialInts[i] != d2->monomialInts[i])
		    {
		      if(d1->monomialInts[i] > d2->monomialInts[i])
			decision = 1;
		      else
			decision = 2;
		    }
		} /* END IF (decision == 0) */
	    } /* END FOR i */
	  if(decision == 0)
	    {
                /* check exponents */
	      if(d1->exponent > d2->exponent)
		decision = 1;
	      else
		decision = 2;
	    }
	}
      else
	{
	  if(i1 == n1)
	      decision = 2;
	  else
	      decision = 1;
	}
      if(decision <= 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "(decision <= 0)");
	  return -1;
	}
      if(decision == 1)
	{
	  memcpy(&workList[nn], &list[i1], sizeof(DistributionSpecStruct));
	  i1++;
	}
      else
	{
	  memcpy(&workList[nn], &list[n1+i2], sizeof(DistributionSpecStruct));
	  i2++;
	}
      nn++;
    } /* END WHILE (nn < n) */
  if(i1 != n1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "(i1 != n1)");
      return -1;
    }
  if(i2 != n2)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "(i2 != n2)");
      return -1;
    }
  if(nn != n)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "(nn != n)");
      return -1;
    }
  memcpy(list, workList, n * sizeof(DistributionSpecStruct));
  return 0;
} /* END do_merge_sort_distrs */








int get_density(const BasisInfoStruct & basisInfo,
		const ergo_real* dmat,
		ergo_real cutoff,
		int maxCountRho,
		DistributionSpecStruct* resultRho)
{
#define MAX_DISTR_IN_TEMP_LIST 888

  int i, j, k, kk;
  DistributionSpecStruct* workList;
  DistributionSpecStruct* rhoSaved;
  ergo_real absvalue;
  ergo_real absdiff;
  ergo_real sqrtValue;
  int sameYesNo, firstIndex, count, withinLimit, resultCount;
  ergo_real coeffSum;
  int* markList;
  int symmetryFactor;
  int nBasisFuncs, nn, nNeededForRho;
  DistributionSpecStruct* rho;

  nNeededForRho = maxCountRho;

  /* allocate rho */
  //rho = (DistributionSpecStruct*)ergo_malloc(nNeededForRho * sizeof(DistributionSpecStruct));
  rho = resultRho;
  
  nBasisFuncs = basisInfo.noOfBasisFuncs;
  nn = 0;
  for(i = 0; i < nBasisFuncs; i++)
    {
      for(j = 0; j < nBasisFuncs; j++)
	{
	  DistributionSpecStruct tempList[MAX_DISTR_IN_TEMP_LIST];
	  int nPrimitives, k;
	  /* the matrix M is symmetric: include diagonal terms once, */
	  /* and include upper off-diagonal terms multiplied by 2 */
	  if(i == j)
              symmetryFactor = 1;
	  else
	    symmetryFactor = 2;
	  if(i > j)
	    continue;
          nPrimitives = 
	    get_product_simple_primitives(basisInfo, i,
					  basisInfo, j,
					  tempList,
					  MAX_DISTR_IN_TEMP_LIST,
					  0);
	  if(nPrimitives <= 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in get_product_simple_primitives");
	      return -1;
	    }
	  for(k = 0; k < nPrimitives; k++)
	    {
	      DistributionSpecStruct* currDistr = &tempList[k];
	      ergo_real Mij = dmat[i*nBasisFuncs+j];
	      ergo_real newCoeff = currDistr->coeff * Mij * symmetryFactor;
	      if(template_blas_fabs(newCoeff) > cutoff)
		{
		  /* add to final list */
		  if(nn > nNeededForRho)
		    {
		      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: (nn > nNeededForRho)");
		      return -1;
		    }
		  memcpy(&rho[nn], currDistr, 
			 sizeof(DistributionSpecStruct));
		  rho[nn].coeff = newCoeff;
		  nn++;  
		}
	    }
	}
    }



  /* Now all distributions are stored in the list 'rho'. */
  /* The number of entries in the list is nn. */
  /* It could happen that all entries are not unique. */
  /* We want to join distributions that have the same center  */
  /* and the same exponent. */
  /* To do this, start with sorting the list by nx, ny, nz, exponent. */
  workList = (DistributionSpecStruct*)ergo_malloc(nn * sizeof(DistributionSpecStruct));
  rhoSaved = (DistributionSpecStruct*)ergo_malloc(nn * sizeof(DistributionSpecStruct));
  memcpy(rhoSaved, rho, nn * sizeof(DistributionSpecStruct));
  
  if(do_merge_sort_distrs(nn, rho, workList) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in do_merge_sort_distrs");
      return -1;
    }
  

  /* check that list is sorted */
  for(i = 0; i < (nn-1); i++)
    {
      if(rho[i].exponent < rho[i+1].exponent)
	{
	  sameYesNo = 1;
	  for(j = 0; j < 3; j++)
	    {
	      if(rho[i].monomialInts[j] != rho[i+1].monomialInts[j])
		sameYesNo = 0;
	    } /* END FOR j */
	  if(sameYesNo == 1)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: distr list NOT properly sorted");
	      return -1;
	    }
	}
    } /* END FOR i */


  markList = (int*)ergo_malloc(nn * sizeof(int));
  for(i = 0; i < nn; i++)
    markList[i] = 0;

  /* now go through sorted list, joining distributions where possible */
  i = 0;
  count = 0;
  firstIndex = 0;
  while(i < nn)
    {
        /* check if this entry has the same nx ny nz as current 'firstIndex' */
      sameYesNo = 1;
      for(j = 0; j < 3; j++)
	{
	  if(rho[i].monomialInts[j] != rho[firstIndex].monomialInts[j])
	    sameYesNo = 0;
	} /* END FOR j */
      /* check exponent */
      absdiff = template_blas_fabs(rho[i].exponent - rho[firstIndex].exponent);
      if(absdiff > EXPONENT_DIFF_LIMIT)
	sameYesNo = 0;
      if(sameYesNo == 0)
	{
	  for(j = firstIndex; j < i; j++)
	    {
	      if(markList[j] == 0)
		{
		  markList[j] = 1;
		  /* join distrs that have centers within  */
		  /* DISTR_CENTER_DIST_LIMIT of this one */
		  coeffSum = rho[j].coeff;
		  for(k = j+1; k < i; k++)
		    {
		      withinLimit = 1;
		      for(kk = 0; kk < 3; kk++)
			{
			  absdiff = template_blas_fabs(rho[j].centerCoords[kk] - 
					 rho[k].centerCoords[kk]);
			  if(absdiff > DISTR_CENTER_DIST_LIMIT)
			    withinLimit = 0;
			} /* END FOR kk */
		      if(withinLimit == 1)
			{
			  coeffSum += rho[k].coeff;
			  markList[k] = 1;
			}
		    } /* END FOR k */
		  memcpy(&workList[count], 
			 &rho[j], 
			 sizeof(DistributionSpecStruct));
		  workList[count].coeff = coeffSum;
		  count++;
		} /* END IF (markList[j] == 0) */
	    } /* END FOR j */
	  firstIndex = i;
	}
      else
	{

	}
      i++;
    } /* END WHILE (i < nn) */
  /* take care of last part */
  for(j = firstIndex; j < nn; j++)
    {
      if(markList[j] == 0)
	{
	  markList[j] = 1;
	  /* join distrs that have centers within  */
	  /* DISTR_CENTER_DIST_LIMIT of this one */
	  coeffSum = rho[j].coeff;
	  for(k = j+1; k < nn; k++)
	    {
	      withinLimit = 1;
	      for(kk = 0; kk < 3; kk++)
		{
		  absdiff = template_blas_fabs(rho[j].centerCoords[kk] - 
				 rho[k].centerCoords[kk]);
		  if(absdiff > DISTR_CENTER_DIST_LIMIT)
		    withinLimit = 0;
		} /* END FOR kk */
	      if(withinLimit == 1)
		{
		  coeffSum += rho[k].coeff;
		  markList[k] = 1;
		}
	    } /* END FOR k */
	  memcpy(&workList[count], &rho[j], sizeof(DistributionSpecStruct));
	  workList[count].coeff = coeffSum;
	  count++;
	} /* END IF (markList[j] == 0) */
    } /* END FOR j */

  for(j = 0; j < nn; j++)
    {
      if(markList[j] != 1)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error: (markList[%i] != 1)", j);
	  return -1;
	}
    } /* END FOR j */


  /* now move results back to list 'rho',  */
  /* skipping those that have too small coeff */
  resultCount = 0;
  for(i = 0; i < count; i++)
    {
      sqrtValue = template_blas_sqrt(pi / workList[i].exponent);
      absvalue = workList[i].coeff * sqrtValue * sqrtValue * sqrtValue;
      if(absvalue < 0) absvalue *= -1;      
      if(absvalue > cutoff)
	{
	  memcpy(&rho[resultCount], 
		 &workList[i], 
		 sizeof(DistributionSpecStruct));
	  resultCount++;
	}
    } /* END FOR i */

  ergo_free(workList);
  ergo_free(markList);
  ergo_free(rhoSaved);

  return resultCount;

}


