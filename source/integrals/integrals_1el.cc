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

/** @file integrals_1el.cc

    @brief Code for 1-electron integrals; dense matrix variant of
    H_core matrix computation.

    @author: Elias Rudberg <em>responsible</em>
*/

/* Written by Elias Rudberg */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <memory.h>
#include <time.h>
#include <stdarg.h>
#include <vector>
#include "integrals_1el.h"
#include "integrals_1el_kinetic.h"
#include "integrals_1el_potential.h"
#include "memorymanag.h"
#include "pi.h"
#include "output.h"
#include "utilities.h"
#include "integral_info.h"


 
int 
compute_h_core_matrix_full(const IntegralInfo& integralInfo,
			   const BasisInfoStruct& basisInfo,
			   int nAtoms,
			   const Atom* atomList,
			   ergo_real* result,
			   ergo_real threshold)
{
  Util::TimeMeter timeMeterTot;
  int n = basisInfo.noOfBasisFuncs;

  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "entering compute_h_core_matrix_full, nAtoms = %i, n = %i, threshold = %g", nAtoms, n, (double)threshold);
  
  /* compute T */
  std::vector<ergo_real> T(n*n);
  for(int i = 0 ; i < n*n; i++)
    T[i] = 0;
  Util::TimeMeter timeMeterT;
  if(compute_T_matrix_full(basisInfo, threshold, &T[0]) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_T_matrix\n");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "compute_T_matrix ending OK.");
  timeMeterT.print(LOG_AREA_INTEGRALS, "compute_T_matrix_full");
  
  /* compute V */
  std::vector<ergo_real> V(n*n);
  for(int i = 0 ; i < n*n; i++)
    V[i] = 0;
  if(compute_V_matrix_full(basisInfo, integralInfo, nAtoms, atomList, threshold, &V[0]) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_V_matrix\n");
      return -1;
    }

  for(int i = 0 ; i < n*n; i++)
    result[i] = T[i] + V[i];
  
  do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, 
	    "compute_h_core_matrix_full ending OK.");
  timeMeterTot.print(LOG_AREA_INTEGRALS, "compute_h_core_matrix_full");
  
  return 0;
}


