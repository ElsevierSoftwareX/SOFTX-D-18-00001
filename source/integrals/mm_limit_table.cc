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

/** @file mm_limit_table.cc

    @brief MMLimitTable class used to predict the magnitude of
    contributions when using truncated multipole expansions.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdio.h>
#include <vector>
#include <limits>
#include <stdexcept>

#include "mm_limit_table.h"

#include "matrix_norm.h"
#include "template_blas_basicmath.h"
#include "multipole.h"


const ergo_real HUGE_REAL_NUMBER = template_blas_sqrt(template_blas_get_num_limit_max<ergo_real>()) / 100;
const ergo_real INITIAL_STEP = 0.5;
const ergo_real RANGE_STEP_DIFF_FACTOR = 0.25;
const ergo_real MAX_ALLOWED_DISTANCE = 1e5;


MMLimitTable::MMLimitTable()
{
  noOfRangesUsed = 0;
}

MMLimitTable::~MMLimitTable()
{
}

void MMLimitTable::inittt(const MultipolePrepManager & multipolePrep)
{
  ergo_real maxDistance_input = MAX_ALLOWED_DISTANCE;
  ergo_real maxDistance = maxDistance_input + 0.1; // ELIAS NOTE 2013-11-29: add small value here to avoid problems at exactly the maxDistance_input distance.
  ergo_real r = 0;
  ergo_real currStep = INITIAL_STEP;
  int rangeCount = 0;

  const int NO_OF_SAMPLE_POINTS = 7;
  ergo_real dxlist[NO_OF_SAMPLE_POINTS][3];
  dxlist[0][0] = 1;
  dxlist[0][1] = 0;
  dxlist[0][2] = 0;
  dxlist[1][0] = 0;
  dxlist[1][1] = 1;
  dxlist[1][2] = 0;
  dxlist[2][0] = 0;
  dxlist[2][1] = 0;
  dxlist[2][2] = 1;
  dxlist[3][0] = 1;
  dxlist[3][1] = 1;
  dxlist[3][2] = 0;
  dxlist[4][0] = 1;
  dxlist[4][1] = 0;
  dxlist[4][2] = 1;
  dxlist[5][0] = 0;
  dxlist[5][1] = 1;
  dxlist[5][2] = 1;
  dxlist[6][0] = 1;
  dxlist[6][1] = 1;
  dxlist[6][2] = 1;
  MMInteractor interactor(multipolePrep);

  while(r < maxDistance)
    {
      interaction_matrix_limit_range_struct & range = rangeList[rangeCount];
      range.startDistance = r;
      range.step = currStep;
      for(int i = 0; i < NO_OF_STEPS_PER_RANGE; i++)
	{
	  r = range.startDistance + i*range.step;
	  range.maxDistance = r + range.step;

	  for(int l_large = 0; l_large <= MAX_MULTIPOLE_DEGREE; l_large++)
	    for(int l_small = 0; l_small <= MAX_MULTIPOLE_DEGREE_BASIC; l_small++)
	      range.list[i].x[l_large][l_small] = 0;
      
	  for(int randloop = 0; randloop < NO_OF_SAMPLE_POINTS; randloop++)
	    {
	      ergo_real dx = dxlist[randloop][0];
	      ergo_real dy = dxlist[randloop][1];
	      ergo_real dz = dxlist[randloop][2];
	      ergo_real norm = template_blas_sqrt(dx*dx+dy*dy+dz*dz);
	      dx *= r / norm;
	      dy *= r / norm;
	      dz *= r / norm;
	      
	      ergo_real T[MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC*MAX_NO_OF_MOMENTS_PER_MULTIPOLE];
	      if(r > 0) {
		interactor.getInteractionMatrix(dx, dy, dz,
						MAX_MULTIPOLE_DEGREE_BASIC,
						MAX_MULTIPOLE_DEGREE, T);
	      }
	      else
		{
		  // For r=0 use huge values for all elements in T.
		  for(int k = 0; k < MAX_NO_OF_MOMENTS_PER_MULTIPOLE_BASIC*MAX_NO_OF_MOMENTS_PER_MULTIPOLE; k++)
		    T[k] = HUGE_REAL_NUMBER;
		}
	  
	      // Compute norms for submatrices of T
	      for(int l_large = 0; l_large <= MAX_MULTIPOLE_DEGREE; l_large++)
		{
		  int startIndex_large = l_large*l_large;
		  int endIndex_large = (l_large+1)*(l_large+1);
		  int n_large = endIndex_large - startIndex_large;
		  for(int l_small = 0; l_small <= MAX_MULTIPOLE_DEGREE_BASIC; l_small++)
		    {
		      int startIndex_small = l_small*l_small;
		      int endIndex_small = (l_small+1)*(l_small+1);
		      int n_small = endIndex_small - startIndex_small;

		      ergo_real T_sub[n_small*n_large];
		      for(int ii = 0; ii < n_large; ii++)
			for(int jj = 0; jj < n_small; jj++)
			  T_sub[ii*n_small+jj] = T[(startIndex_small+jj)*MAX_NO_OF_MOMENTS_PER_MULTIPOLE+startIndex_large+ii];

		      ergo_real matrixNorm = get_euclidean_norm(n_large, n_small, T_sub);

		      if(matrixNorm > range.list[i].x[l_large][l_small])
			range.list[i].x[l_large][l_small] = matrixNorm;
		  
		    } // END FOR l_small
		} // END FOR l_large
	    } // END FOR randloop
	  
	} // END FOR i
      currStep = r * RANGE_STEP_DIFF_FACTOR;
      rangeCount++;
      if(rangeCount >= NO_OF_RANGES)
	throw "error in MMLimitTable::Init: (rangeCount >= NO_OF_RANGES)";
    } // END WHILE
  noOfRangesUsed = rangeCount;
}


const MMLimitTable::interaction_matrix_limit_struct & MMLimitTable::get_x_from_distance(ergo_real distance) const
{
  if(distance < 0)
    throw "Error in MMLimitTable::get_x_from_distance(): (distance < 0).";
  int rangeIndex = 0;
  while(rangeList[rangeIndex].maxDistance <= distance) // ELIAS NOTE 2013-11-29: changed from "<" to "<=" here to fix problem found by Anastasia, for case when (distanceLeft / range.step) divides evenly.
    {
      rangeIndex++;
      if(rangeIndex >= noOfRangesUsed)
	throw "error in MMLimitTable::get_x_from_distance: (rangeIndex >= noOfRangesUsed)";
    }
  const interaction_matrix_limit_range_struct & range = rangeList[rangeIndex];
  ergo_real distanceLeft = distance - range.startDistance;
  int i = (int)(distanceLeft / range.step);
  if(i < 0 || i >= NO_OF_STEPS_PER_RANGE)
    {
      if(i < 0)
	throw "error in MMLimitTable::get_x_from_distance: i < 0";
      throw "error in MMLimitTable::get_x_from_distance: i >= NO_OF_STEPS_PER_RANGE";
    }
  const interaction_matrix_limit_struct & x = range.list[i];
  return x;
}


ergo_real MMLimitTable::get_max_abs_mm_contrib(int degree1,
					       const ergo_real* maxMomentVectorNormList1,
					       int degree2,
					       const ergo_real* maxMomentVectorNormList2,
					       ergo_real distance) const
{
  ergo_real maxAbsContributionFromMultipole = 0;
  // Get worst-case interaction matrix limits
  const interaction_matrix_limit_struct & x = get_x_from_distance(distance);
  for(int l_large = degree1; l_large >= 0; l_large--)
    {
      ergo_real contribThisDegree = 0;
      for(int l_small = 0; l_small <= degree2; l_small++)
	{
	  contribThisDegree += 
	    maxMomentVectorNormList1[l_small] *
	    maxMomentVectorNormList2[l_large] *
	    x.x[l_large][l_small];
	} // END FOR l_small
      maxAbsContributionFromMultipole += contribThisDegree;
    } // END FOR l_large
  return maxAbsContributionFromMultipole;
}

int MMLimitTable::get_minimum_multipole_degree_needed(ergo_real distance,
						      const multipole_struct_large* boxMultipole, 
						      int maxDegreeForDistrs, 
						      const ergo_real* maxMomentVectorNormForDistrsList, 
						      ergo_real threshold) const
{
  // Get worst-case interaction matrix limits
  const interaction_matrix_limit_struct & x = get_x_from_distance(distance);
  ergo_real maxAbsContribution = 0;
  int degreeNeeded = boxMultipole->degree;
  for(int l_large = boxMultipole->degree; l_large >= 0; l_large--)
    {
      degreeNeeded = l_large;
      ergo_real contribThisDegree = 0;
      for(int l_small = 0; l_small <= maxDegreeForDistrs; l_small++)
	{
	  contribThisDegree += 
	    maxMomentVectorNormForDistrsList[l_small] * 
	    boxMultipole->euclideanNormList[l_large] * 
	    x.x[l_large][l_small];
	} // END FOR l_small
      maxAbsContribution += contribThisDegree;
      if(maxAbsContribution > threshold)
	break;
    } // END FOR l_large
  return degreeNeeded;
}

void MMLimitTable::write_to_buffer ( char * dataBuffer, size_t const bufferSize ) const {
  char* p = dataBuffer;
  if(bufferSize < get_size())
    throw std::runtime_error("Error in MMLimitTable::write_to_buffer: bufferSize too small.");
  memcpy(p, this, sizeof(MMLimitTable));
}

size_t MMLimitTable::get_size() const {
  return sizeof(MMLimitTable);
}

void MMLimitTable::assign_from_buffer ( char const * dataBuffer, size_t const bufferSize) {
  const char* p = dataBuffer;
  if(bufferSize < sizeof(MMLimitTable))
    throw std::runtime_error("Error in MMLimitTable::assign_from_buffer: bufferSize too small.");
  memcpy(this, p, sizeof(MMLimitTable));
}

