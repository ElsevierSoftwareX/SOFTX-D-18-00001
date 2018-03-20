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

/** @file operator_matrix.cc

    @brief Functions for computing the matrix of a
    dipole/quadrupole/etc operator. Full and sparse versions.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include <vector>
#include "pi.h"
#include "output.h"
#include "utilities.h"
#include "boysfunction.h"
#include "integral_info.h"
#include "integrals_general.h"
#include "box_system.h"
#include "operator_matrix.h"


static const ergo_real MATRIX_ELEMENT_THRESHOLD_VALUE = 1e-12; // FIXME: SHOULD NOT USE HARD-CODED VALUE HERE!


int
compute_operator_matrix_full(const BasisInfoStruct & basisInfoA, 
			     const BasisInfoStruct & basisInfoB,
			     int pow_x,
			     int pow_y,
			     int pow_z,
			     ergo_real* result)
{
  int n_A = basisInfoA.noOfBasisFuncs;
  int n_B = basisInfoB.noOfBasisFuncs;
  std::vector<int> nvaluesList(n_A);
  std::vector< std::vector<int> > colindList(n_A);
  std::vector< std::vector<ergo_real> > valuesList(n_A);
  
  if(compute_operator_matrix_sparse(basisInfoA, 
				    basisInfoB,
				    pow_x,
				    pow_y,
				    pow_z,
				    n_A,
				    n_B,
				    nvaluesList,
				    colindList,
				    valuesList) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_operator_matrix_sparse");
      return -1;
    }
  
  // Now populate full result matrix
  memset(result, 0, n_A*n_B*sizeof(ergo_real));
  for(int i = 0; i < n_A; i++)
    {
      int count = nvaluesList[i];
      const std::vector<int> & colind = colindList[i];
      const std::vector<ergo_real> & values = valuesList[i];
      for(int j = 0; j < count; j++)
	result[i*n_B+colind[j]] = values[j];
    } // END FOR i

  return 0;
}


/** computes the matrix of a dipole/quadrupole/etc operator. The
    columns and rows enumerate @param basisInfoA and @param
    basisInfoB respectively. The operator is in the form: X =
    (x^pow_x*y^pow_y*z^pow_z). The resulting matrix (possibly
    rectangular) is returned in nvaluesList, colindList, valuesList. 
    Overlap matrix is
    associated with triple (0,0,0), X component of the dipole moment
    with (1,0,0), etc.
    The parameters @param pow_x @param pow_y @param pow_z determine the operator.
    The parameters @param n_A @param n_B give the number of basis functions in each of the two basis sets.
    The result is stored using the lists @param nvaluesList @param colindList @param valuesList each having length n_A.
*/
int 
compute_operator_matrix_sparse(const BasisInfoStruct & basisInfoA, 
			       const BasisInfoStruct & basisInfoB,
			       int pow_x,
			       int pow_y,
			       int pow_z,
			       int n_A,
			       int n_B,
			       std::vector<int> & nvaluesList,       // length n_A
			       std::vector< std::vector<int> > & colindList,       // length n_A, each element will be allocated
			       std::vector< std::vector<ergo_real> > & valuesList  // length n_A, each element will be allocated
			       )
{
  int internal_error = 0;
  int nBastA = basisInfoA.noOfBasisFuncs;
  int nBastB = basisInfoB.noOfBasisFuncs;

  if(n_A != nBastA || n_B != nBastB)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_operator_matrix_sparse: (n_A != nBastA || n_B != nBastB)");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "compute_operator_matrix_sparse, nBastA = %6i, nBastB = %6i, pows (x y z) = (%i %i %i)", 
	    nBastA, nBastB, pow_x, pow_y, pow_z);

  Util::TimeMeter timeMeter;

  // To reduce scaling we want some kind of "extent" for each basis function.
  // Start by getting largest simple integral for each of the two basis sets.
  ergo_real A_A = get_largest_simple_integral(basisInfoA);
  ergo_real A_B = get_largest_simple_integral(basisInfoB);
  std::vector<ergo_real> basisFuncExtentList_A(nBastA);
  std::vector<ergo_real> basisFuncExtentList_B(nBastB);
  get_basis_func_extent_list(basisInfoA, &basisFuncExtentList_A[0], MATRIX_ELEMENT_THRESHOLD_VALUE / A_A);
  get_basis_func_extent_list(basisInfoB, &basisFuncExtentList_B[0], MATRIX_ELEMENT_THRESHOLD_VALUE / A_B);
  ergo_real maxExtentB = 0;
  for(int i = 0; i < nBastB; i++)
    {
      ergo_real currExtent = basisFuncExtentList_B[i];
      if(currExtent > maxExtentB)
	maxExtentB = currExtent;
    }
  
  // Create box system for basisInfoB.
  std::vector<box_item_struct> itemList(nBastB);
  for(int i = 0; i < nBastB; i++)
    {
      for(int j = 0; j < 3; j++)
	itemList[i].centerCoords[j] = basisInfoB.basisFuncList[i].centerCoords[j];
      itemList[i].originalIndex = i;
    }
  ergo_real toplevelBoxSize = 7.0;
  BoxSystem boxSystem;
  if(boxSystem.create_box_system(&itemList[0],
				 nBastB,
				 toplevelBoxSize) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_operator_matrix_sparse: error creating box system.");
      return -1;
    }
  
  static const int maxDistrsInTempList = 40000;
  static const int maxDistrsInTempList2 = 400;

  int operatorMonomialInts[3];
  operatorMonomialInts[0] = pow_x;
  operatorMonomialInts[1] = pow_y;
  operatorMonomialInts[2] = pow_z;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Allocate vector for results for one row.
    std::vector<ergo_real> rowValueList(nBastB);
    std::vector<DistributionSpecStruct> tempList(maxDistrsInTempList);
    std::vector<int> orgIndexList(nBastB);
#ifdef _OPENMP
#pragma omp for
#endif
  for(int i = 0; i < nBastA; i++)
    {
      int count = 0;
      // Now, instead of looping over all nBastB basis functions, we use box system to find relevant ones.
      ergo_real maxDistance = basisFuncExtentList_A[i] + maxExtentB;
      ergo_real coords[3];
      for(int coordNo = 0; coordNo < 3; coordNo++)
	coords[coordNo] = basisInfoA.basisFuncList[i].centerCoords[coordNo];
      int nRelevant = boxSystem.get_items_near_point(&itemList[0], coords, maxDistance, &orgIndexList[0]);
      for(int jRelevant = 0; jRelevant < nRelevant; jRelevant++)
	{
	  int j = orgIndexList[jRelevant];
	  // Compute distance between basis function centers
	  ergo_real dx = basisInfoA.basisFuncList[i].centerCoords[0] - basisInfoB.basisFuncList[j].centerCoords[0];
	  ergo_real dy = basisInfoA.basisFuncList[i].centerCoords[1] - basisInfoB.basisFuncList[j].centerCoords[1];
	  ergo_real dz = basisInfoA.basisFuncList[i].centerCoords[2] - basisInfoB.basisFuncList[j].centerCoords[2];
	  ergo_real distance = template_blas_sqrt(dx*dx + dy*dy + dz*dz);
	  // We can skip if distance is greater than sum of extents.
	  if(distance > basisFuncExtentList_A[i] + basisFuncExtentList_B[j]) {
            rowValueList[jRelevant] = 0.0;
	    continue;
          }

	  int nPrimitives = 
	    get_product_simple_primitives(basisInfoA, i,
					  basisInfoB, j,
					  &tempList[0],
					  maxDistrsInTempList,
					  0);
	  if(nPrimitives <= 0)
	    {
	      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in get_product_simple_primitives\n");
              internal_error++;
              goto internal_error_occured;
	    }

	  ergo_real sum = 0;
	  for(int k = 0; k < nPrimitives; k++)
	    {
	      DistributionSpecStruct* currDistr = &tempList[k];
	      // now we need to multiply the current distribution by the operator monomial
	      // the result will be a list of new distributions
	      DistributionSpecStruct tempList2[maxDistrsInTempList2];
	      // first put the distribution as the only entry in tempList2
	      // then loop over operator monomial, each time creating a new list in tempList3
	      // and move list back to tempList2 each time
	      memcpy(&tempList2[0], currDistr, sizeof(DistributionSpecStruct));
	      int tempList2_count = 1;
	      for(int coordNo = 0; coordNo < 3; coordNo++)
		{
		  for(int ii = 0; ii < operatorMonomialInts[coordNo]; ii++)
		    {
		      DistributionSpecStruct tempList3[maxDistrsInTempList2];
		      int tempList3_count = 0;
		      // now go through tempList2, and for each entry create two new entries in tempList3
		      for(int jj = 0; jj < tempList2_count; jj++)
			{
			  // multiply this distribution by a single coordinate, x y or z according to coordNo
			  // this gives two new distributions
			  // the first one is the same as the original one multiplied by a constant
			  memcpy(&tempList3[tempList3_count], &tempList2[jj], sizeof(DistributionSpecStruct));
			  tempList3[tempList3_count].coeff *= tempList3[tempList3_count].centerCoords[coordNo];
			  tempList3_count++;
			  if(tempList3_count >= maxDistrsInTempList2)
			    {
			      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_operator_matrix: (tempList3_count >= maxDistrsInTempList2)");
			      internal_error++;
                              goto internal_error_occured;
			    }
			  // the second one is the same as the original one with increased pow of current coordinate
			  memcpy(&tempList3[tempList3_count], &tempList2[jj], sizeof(DistributionSpecStruct));
			  tempList3[tempList3_count].monomialInts[coordNo]++;
			  tempList3_count++;
			  if(tempList3_count >= maxDistrsInTempList2)
			    {
			      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_operator_matrix: (tempList3_count >= maxDistrsInTempList2)");
                              internal_error++;
                              goto internal_error_occured;
			    }
			} // END FOR jj
		      // now tempList3 is complete. copy it back to tempList2
		      memcpy(tempList2, tempList3, tempList3_count * sizeof(DistributionSpecStruct));
		      tempList2_count = tempList3_count;
		    } // END FOR ii
		} // END FOR coordNo
	      // now tempList2 contains all the final distributions
	      for(int ii = 0; ii < tempList2_count; ii++)
		sum += compute_integral_of_simple_prim(tempList2[ii]);
	    } /* END FOR k */
	  rowValueList[jRelevant] = sum;
	  if(template_blas_fabs(sum) > MATRIX_ELEMENT_THRESHOLD_VALUE)
	    count++;
	} /* END FOR jRelevant */
      // OK, this row done.

      nvaluesList[i] = count;
      // Now allocate result vectors for this row.
      colindList[i].resize(count);
      valuesList[i].resize(count);
      count = 0;
      for(int jRelevant = 0; jRelevant < nRelevant; jRelevant++)
	{
	  int j = orgIndexList[jRelevant];
	  ergo_real absVal = template_blas_fabs(rowValueList[jRelevant]);
	  if(absVal > MATRIX_ELEMENT_THRESHOLD_VALUE)
	    {
	      if(count >= nvaluesList[i]) {
                internal_error++;
                goto internal_error_occured;
              }
	      colindList[i][count] = j;
	      valuesList[i][count] = rowValueList[jRelevant];
	      count++;
	    }
	}
    internal_error_occured:;
    } /* END FOR i */
  }

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "compute_operator_matrix_sparse finished.");
  timeMeter.print(LOG_AREA_INTEGRALS, "compute_operator_matrix_sparse");
  
  return -internal_error;
}


int
compute_overlap_matrix(const BasisInfoStruct & basisInfoA, 
		       const BasisInfoStruct & basisInfoB,
		       ergo_real* result)
{
  return compute_operator_matrix_full(basisInfoA, basisInfoB, 0, 0, 0, result);
}


