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

/** @file integrals_general.cc

    @brief General functionality related to computation of integrals
    involving Gaussian basis functions.

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
#include "memorymanag.h"
#include "pi.h"
#include "output.h"
#include "utilities.h"
#include "boysfunction.h"
#include "integral_info.h"
#include "integrals_general.h"



#define K_MAX_DIM 44


int 
multiply_polynomials(ergo_real result[], 
		     polydeg1struct* polydeg1, 
		     int dim, 
		     ergo_real a[])
{
  int i;
  ergo_real p1[K_MAX_DIM + 1];
  ergo_real p2[K_MAX_DIM + 1];
  if(dim >= (K_MAX_DIM-1))
    return -1;
  for(i = 0; i <= dim; i++)
    p1[i] = a[i]*polydeg1->a0;
  p1[dim+1] = 0;
  p2[0] = 0;
  for(i = 0; i <= dim; i++)
    p2[i+1] = a[i]*polydeg1->a1;
  for(i = 0; i <= (dim+1); i++)
    result[i] = p1[i] + p2[i];
  return 0;
} /* END multiply_polynomials */





/*
get_product_simple_prims
This function calculates the product of two simple primitives.
The result is a list of simple primitives.
*/
int
get_product_simple_prims(const DistributionSpecStruct& primA_in,
			 const DistributionSpecStruct& primB_in,
			 DistributionSpecStruct resultList[],
			 int maxCount,
			 ergo_real threshold)
{
  // Use a coordinate system with primA at the origin.
  // This solves the problem with extreme positions of the primitives.
  DistributionSpecStruct primA_mod = primA_in;
  DistributionSpecStruct primB_mod = primB_in;
  int kk;
  for(kk = 0; kk < 3; kk++)
    {
      primA_mod.centerCoords[kk] -= primA_in.centerCoords[kk];
      primB_mod.centerCoords[kk] -= primA_in.centerCoords[kk];
    }
  DistributionSpecStruct* primA = &primA_mod;
  DistributionSpecStruct* primB = &primB_mod;

  ergo_real CxCyCz, AiAj, alphaNew;
  ergo_real newCenter[3];
  ergo_real poly0[K_MAX_DIM];
  ergo_real poly1[K_MAX_DIM];
  ergo_real poly2[K_MAX_DIM];
  ergo_real tempPoly[K_MAX_DIM];
  ergo_real tempPoly2[K_MAX_DIM];
  ergo_real tempPoly3[K_MAX_DIM];
  int tempPolyDegree, tempPoly2Degree;
  int poly0degree, poly1degree, poly2degree, l, m, nn;
  polydeg1struct polyDeg1;
  ergo_real* poly;
  int* degreePtr;
  /* use the Gaussian product rule */
  ergo_real sum = 0;
  int k;
  for(k = 0; k < 3; k++)
    {
      ergo_real temp = primA->centerCoords[k] - primB->centerCoords[k];
      sum += temp * temp;
    } /* END FOR k */
  CxCyCz = template_blas_exp(-primA->exponent * primB->exponent * 
	       sum / (primA->exponent + primB->exponent));

  // FIXME: do this screening properly!
  if(template_blas_fabs(CxCyCz) < threshold)
    return 0;

  AiAj = primA->coeff * primB->coeff;
  alphaNew = primA->exponent + primB->exponent;
  for(k = 0; k < 3; k++)
    {
      newCenter[k] = 
	(primA->exponent * primA->centerCoords[k] +
	 primB->exponent * primB->centerCoords[k]) /
	(primA->exponent + primB->exponent);
    } /* END FOR k */

  /* do product of polynomials */
  /* one coordinate at a time */
  for(k = 0; k < 3; k++)
    {
      switch(k)
	{
	case 0: poly = poly0; degreePtr = &poly0degree; break;
	case 1: poly = poly1; degreePtr = &poly1degree; break;
	case 2: poly = poly2; degreePtr = &poly2degree; break;
	default: return -1;
	} /* END SWITCH k */
      tempPoly[0] = 1;
      tempPolyDegree = 0;
      for(m = 0; m < primA->monomialInts[k]; m++)
	{
	  polyDeg1.a0 = -primA->centerCoords[k];
	  polyDeg1.a1 = 1;
	  if(multiply_polynomials(tempPoly2, &polyDeg1, 
				  tempPolyDegree, tempPoly) != 0)
	    return -1;
	  tempPolyDegree++;
	  memcpy(tempPoly, 
		 tempPoly2, 
		 (tempPolyDegree+1)*sizeof(ergo_real));
	} /* END FOR m */
      for(m = 0; m < primB->monomialInts[k]; m++)
	{
	  polyDeg1.a0 = -primB->centerCoords[k];
	  polyDeg1.a1 = 1;
	  if(multiply_polynomials(tempPoly2, &polyDeg1, 
				  tempPolyDegree, tempPoly) != 0)
	    return -1;
	  tempPolyDegree++;
	  memcpy(tempPoly, 
		 tempPoly2, 
		 (tempPolyDegree+1)*sizeof(ergo_real));
	} /* END FOR m */

      /* now do variable change */
      for(m = 0; m < K_MAX_DIM; m++)
	poly[m] = 0;
      tempPoly2Degree = 0;
      for(m = 0; m <= tempPolyDegree; m++)
	{
	  tempPoly2[0] = tempPoly[m];
	  tempPoly2Degree = 0;
	  for(l = 0; l < m; l++)
	    {
	      polyDeg1.a0 = newCenter[k];
	      polyDeg1.a1 = 1;
	      if(multiply_polynomials(tempPoly3, 
				      &polyDeg1, 
				      tempPoly2Degree, 
				      tempPoly2) != 0)
		return -1;
	      tempPoly2Degree++;
	      memcpy(tempPoly2, 
		     tempPoly3, 
		     (tempPoly2Degree+1)*sizeof(ergo_real));
	    } /* END FOR l */
	  for(l = 0; l <= tempPoly2Degree; l++)
	    {
	      poly[l] += tempPoly2[l];
	    } /* END FOR l */
	} /* END FOR m */
      *degreePtr = tempPoly2Degree;
    } /* END FOR k */

  nn = 0;
  for(k = 0; k <= poly0degree; k++)
    {
      int l;
      for(l = 0; l <= poly1degree; l++)
	{
	  int m;
	  for(m = 0; m <= poly2degree; m++)
	    {
	      ergo_real newCoeff = AiAj * CxCyCz * poly0[k] * poly1[l] * poly2[m];

	      ergo_real sqrtValue = template_blas_sqrt(pi / alphaNew);
	      ergo_real absvalue = newCoeff * sqrtValue * sqrtValue * sqrtValue;
	      if(absvalue < 0) absvalue *= -1;

	      /* add one function to final list */
	      resultList[nn].coeff = newCoeff;
	      resultList[nn].exponent = alphaNew;

	      memcpy(resultList[nn].centerCoords, 
		     newCenter, 
		     3 * sizeof(ergo_real));
	      resultList[nn].monomialInts[0] = k;
	      resultList[nn].monomialInts[1] = l;
	      resultList[nn].monomialInts[2] = m;

	      // Translate this term of result back to original coordinate system
	      for(kk = 0; kk < 3; kk++)
		resultList[nn].centerCoords[kk] += primA_in.centerCoords[kk];
	      
	      nn++;
	      if(nn >= maxCount)
		{
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_prims: "
			    "maxCount exceeded");
		  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "nn = %i, maxCount = %i", 
			    nn, maxCount);
		  return -1;
		}
	    } /* END FOR m */
	} /* END FOR l */
    } /* END FOR k */

  return nn;
}






int
get_product_simple_primitives(const BasisInfoStruct & basisInfoA, int iA,
			      const BasisInfoStruct & basisInfoB, int iB,
			      DistributionSpecStruct resultList[],
			      int maxCount,
			      ergo_real threshold)
{
  BasisFuncStruct* basisFuncA = &basisInfoA.basisFuncList[iA];
  int nPrimsA = basisFuncA->noOfSimplePrimitives;
  int Astart = basisFuncA->simplePrimitiveIndex;
  BasisFuncStruct* basisFuncB = &basisInfoB.basisFuncList[iB];
  int nPrimsB = basisFuncB->noOfSimplePrimitives;
  int Bstart = basisFuncB->simplePrimitiveIndex;
  int n = 0;
  int i;
  if((nPrimsA <= 0) || (nPrimsB <= 0))
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives: "
		"((nPrimsA <= 0) || (nPrimsB <= 0))\n");
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "nPrimsA = %i, nPrimsB = %i\n", nPrimsA, nPrimsB);
      return -1;
    }
  for(i = 0; i < nPrimsA; i++)
    {
      const DistributionSpecStruct& primA = 
	basisInfoA.simplePrimitiveList[Astart + i];
      int j;
      for(j = 0; j < nPrimsB; j++)
	{
	  const DistributionSpecStruct& primB = 
	    basisInfoB.simplePrimitiveList[Bstart + j];

	  int nNewPrims = get_product_simple_prims(primA, 
						   primB, 
						   &resultList[n],
						   maxCount - n,
						   threshold);
	  if(nNewPrims < 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_prims");
	      return -1;
	    }
	  
	  n += nNewPrims;
	}
    }
  return n;
}


ergo_real
compute_integral_of_simple_prim(const DistributionSpecStruct & distr) {
  ergo_real result, alpha;
  int i, n;
  if( ((distr.monomialInts[0]|
	distr.monomialInts[1]|
	distr.monomialInts[2]) & 1) == 1) /* odd integrals disappear */
    return 0;
  alpha = distr.exponent;
  result = distr.coeff * template_blas_pow((ergo_real)pi/alpha, (ergo_real)1.5);
  ergo_real twoA = 2*alpha;
  for(i = 0; i < 3; i++) {
    n = distr.monomialInts[i];
    for(int j=0; j<n; j+=2)
      result *= (j+1)/twoA;
  } /* END FOR i */
  return result;
}


/**
  Computes the largest integral of any primitive in the basis set,
  when any x y z factors are ignored. This is useful for getting
  rough estimates of basis function extents.
*/
ergo_real
get_largest_simple_integral(const BasisInfoStruct & basisInfo)
{
  int n = basisInfo.noOfBasisFuncs;
  ergo_real A = 0;
  int i;
  for(i = 0; i < n; i++)
    {
      BasisFuncStruct* basisFunc = &basisInfo.basisFuncList[i];
      int nPrims = basisFunc->noOfSimplePrimitives;
      int start = basisFunc->simplePrimitiveIndex;
      int j;
      for(j = 0; j < nPrims; j++)
	{
	  DistributionSpecStruct* prim = &basisInfo.simplePrimitiveList[start + j];
	  DistributionSpecStruct distr;
	  distr = *prim;
	  // Set monomialInts to zero to simplify things
	  distr.monomialInts[0] = 0;
	  distr.monomialInts[1] = 0;
	  distr.monomialInts[2] = 0;
	  ergo_real a = compute_integral_of_simple_prim(distr);
	  if(a > A)
	    A = a;
	} // END FOR j
    } // END FOR i
  return A;
}



/**
   Computes an estimate for the largest absolute value that any basis
   function takes. Useful as "worst case" when you want to find out
   the largest contribution to the density that a basis function can
   be part of.  */
ergo_real get_max_basis_func_abs_value(const BasisInfoStruct & basisInfo) {
  int n = basisInfo.noOfBasisFuncs;
  ergo_real maxValue = 0;
  for(int i = 0; i < n; i++) {
    BasisFuncStruct* basisFunc = &basisInfo.basisFuncList[i];
    int nPrims = basisFunc->noOfSimplePrimitives;
    int start = basisFunc->simplePrimitiveIndex;
    for(int j = 0; j < nPrims; j++) {
      DistributionSpecStruct* prim = &basisInfo.simplePrimitiveList[start + j];
      ergo_real valueAtCenter = template_blas_fabs(prim->coeff); // exp(0) = 1
      if(valueAtCenter > maxValue)
	maxValue = valueAtCenter;
    } // END FOR j
  } // END FOR i
  return maxValue;
}


/**
  Computes an "extent" for each basis function in the basis set.
  The "extent" is such that the value of the function is smaller
  than maxAbsValue at distances beyond the "extent".
*/
int
get_basis_func_extent_list(const BasisInfoStruct & basisInfo, ergo_real* basisFuncExtentList, ergo_real maxAbsValue)
{
  int n = basisInfo.noOfBasisFuncs;
  for(int i = 0; i < n; i++)
    {
      BasisFuncStruct* basisFunc = &basisInfo.basisFuncList[i];
      int nPrims = basisFunc->noOfSimplePrimitives;
      int start = basisFunc->simplePrimitiveIndex;
      ergo_real maxExtent = 0;
      for(int j = 0; j < nPrims; j++)
	{
	  DistributionSpecStruct* prim = &basisInfo.simplePrimitiveList[start + j];
	  ergo_real currExtent = template_blas_sqrt((1.0 / prim->exponent) * template_blas_log(template_blas_fabs(prim->coeff) / maxAbsValue));
	  if(currExtent > maxExtent)
	    maxExtent = currExtent;
	} // END FOR j
      basisFuncExtentList[i] = maxExtent;
    } // END FOR i
  return 0;
}
