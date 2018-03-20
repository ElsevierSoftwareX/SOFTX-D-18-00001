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

/** @file diis_restricted.cc

    @brief DIISManagerRestricted class implementing direct inversion
    in the iterative subspace (DIIS) for restricted SCF calculations.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include <string.h>

#include "diis_restricted.h"

#include "output.h"
#include "solve_lin_eq_syst.h"
#include "utilities.h"


DIISManagerRestricted::DIISManagerRestricted()
  : DIISManager()
{
}

DIISManagerRestricted::~DIISManagerRestricted()
{
  ClearList();
}

int DIISManagerRestricted::AddIterationToList(symmMatrix & F, normalMatrix & E)
{
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "entering DIISManagerRestricted::AddIterationToList, IterCount = %2i", IterCount);

  Util::TimeMeter timeMeter;
  
  if(IterCount > MaxNoOfIters)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error:(IterCount > MaxNoOfIters)");
      return -1;
    }

  // check if the list is full
  if(IterCount == MaxNoOfIters)
    {
      // remove oldest iteration
      delete F_list[0][0];
      F_list[0][0] = NULL;
      delete E_list[0][0];
      E_list[0][0] = NULL;
      int i;
      for(i = 0; i < IterCount-1; i++)
	{
	  F_list[0][i] = F_list[0][i+1];
	  F_list[0][i+1] = NULL;
	  E_list[0][i] = E_list[0][i+1];
	  E_list[0][i+1] = NULL;
	}
      RemoveOldestIteration();  /* note that this changes the value of IterCount */
    }

  F_list[0][IterCount] = new symmMatrix(F);
  F_list[0][IterCount]->writeToFile();
  E_list[0][IterCount] = new normalMatrix(E);
  E_list[0][IterCount]->writeToFile();

  // Create new B matrix
  int dimB    = IterCount + 1;
  int dimBnew = IterCount + 2;
  ergo_real* Bnew = new ergo_real[dimBnew*dimBnew];
  memset(Bnew, 0, dimBnew*dimBnew*sizeof(ergo_real));
  int i, j;
  for(i = 0; i < dimB; i++)
    for(j = 0; j < dimB; j++)
      Bnew[i*dimBnew+j] = B[i*dimB+j];
  // Set two matrix elements to -1
  Bnew[0*dimBnew+dimBnew-1] = -1;
  Bnew[(dimBnew-1)*dimBnew+0] = -1;
  // Now it remains to complete B with scalar products of error matrices
  for(i = 0; i < IterCount; i++)
    {
      // compute dot product of error matrix i and E
      E_list[0][i]->readFromFile();
      ergo_real scalarProd = DoScalarProductOfErrorMatrices(E, *E_list[0][i]);
      E_list[0][i]->writeToFile();
      Bnew[(dimBnew-1)*dimBnew+(1+i)] = scalarProd;
      Bnew[(1+i)*dimBnew+(dimBnew-1)] = scalarProd;
    }
  // Do scalar product of the new E with itself
  Bnew[(dimBnew-1)*dimBnew+(dimBnew-1)] = DoScalarProductOfErrorMatrices(E, E);

  // Copy Bnew to B
  memcpy(B, Bnew, dimBnew*dimBnew*sizeof(ergo_real));
  
  delete [] Bnew;
  
  IterCount++;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "DIISManagerRestricted::AddIterationToList ending OK.");
  timeMeter.print(LOG_AREA_SCF, "DIISManagerRestricted::AddIterationToList");
  
  return 0;
}

int DIISManagerRestricted::ClearList()
{
  int i;
  for(i = 0; i < IterCount; i++)
    {
      delete F_list[0][i];
      F_list[0][i] = NULL;
      delete E_list[0][i];
      E_list[0][i] = NULL;
    }
  IterCount = 0;
  return 0;
}

int DIISManagerRestricted::GetCombinedFockMatrix(symmMatrix & result)
{
  if(IterCount <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in DIISManagerRestricted::GetCombinedFockMatrix: (IterCount <= 0)");
      return -1;
    }

  int dimB = IterCount + 1;
  ergo_real* RHS = new ergo_real[dimB];
  ergo_real* cVector = new ergo_real[dimB];
  
  // Construct vector RHS
  RHS[0] = -1;
  int i;
  for(i = 0; i < IterCount; i++)
    RHS[i+1] = 0;

  // Solve equation system B*x = HL
  if(solve_linear_equation_system(dimB, B, RHS, cVector) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in solve_linear_equation_system");
      return -1;
    }
  
  // Create linear combination of Fock matrices using coefficients from cVector.
  F_list[0][0]->readFromFile();
  result = *F_list[0][0];
  F_list[0][0]->writeToFile();
  result *= cVector[1];
  for(i = 1; i < IterCount; i++)
    {
      F_list[0][i]->readFromFile();
      result += cVector[1+i] * (*F_list[0][i]);
      F_list[0][i]->writeToFile();
    } // END FOR i
  

  delete [] RHS;
  delete [] cVector;
  
  return 0;
}



