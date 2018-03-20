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

/** @file diis_general.cc

    @brief General DIISManager class that can be used as a base class
    in implementations of direct inversion in the iterative subspace
    (DIIS) for restricted and unrestricted cases.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include <unistd.h>

#include "diis_general.h"

#include "utilities.h"
#include "output.h"
#include "memorymanag.h"


DIISManager::DIISManager()
{
  MaxNoOfIters = 0;
  IterCount = 0;
  B = NULL;
  for(int i = 0; i < 2; i++)
    {
      F_list[i] = NULL;
      E_list[i] = NULL;
    }
}


DIISManager::~DIISManager()
{
  for(int i = 0; i < 2; i++) {
    if(F_list[i])
      delete []F_list[i];
    if(E_list[i])
      delete []E_list[i];
  }
  if(B)
    delete []B;
}


int DIISManager::GetNoOfIters()
{
  return IterCount;
}


typedef symmMatrix* symmMatrixPtr;
typedef normalMatrix* normalMatrixPtr;

int DIISManager::Initialize(int noOfIters)
{
  if(noOfIters <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in DIISManager::Initialize: noOfIters = %i", noOfIters);
      return -1;
    }

  MaxNoOfIters = noOfIters;
  if(noOfIters <= 0)
    return -1;
  B = new ergo_real[(MaxNoOfIters+1)*(MaxNoOfIters+1)];
  // Initially, B is a 1*1 matrix, whose only matrix element should be zero.
  B[0] = 0;

  // Set up lists of matrix pointers
  for(int i = 0; i < 2; i++)
    {
      F_list[i] = new symmMatrixPtr[noOfIters];
      E_list[i] = new normalMatrixPtr[noOfIters];
      int j;
      for(j = 0; j < noOfIters; j++)
	{
	  F_list[i][j] = NULL;
	  E_list[i][j] = NULL;
	}
    }

    return 0;
}



ergo_real DIISManager::DoScalarProductOfErrorMatrices(const normalMatrix & E1, 
						      const normalMatrix & E2)
{
  return normalMatrix::trace_ab(E1, E2);
}

int DIISManager::RemoveOldestIteration()
{
  int dimB    = IterCount + 1;
  int dimBnew = IterCount;
  ergo_real* Bnew = new ergo_real[dimBnew*dimBnew];
  int i, j;
  for(i = 0; i < dimBnew; i++)
    for(j = 0; j < dimBnew; j++)
      {
	if(i == 0 || j == 0)
	  Bnew[i*dimBnew+j] = B[i*dimB+j];
	else
	  Bnew[i*dimBnew+j] = B[(i+1)*dimB+(j+1)];
      }
  memcpy(B, Bnew, dimBnew*dimBnew*sizeof(ergo_real));
  delete [] Bnew;
  IterCount--;
  return 0; /* success */
}


