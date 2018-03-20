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

/** @file integrals_2el_layer.cc

    @brief Functions for computing dense Coulomb and HF exchange
    matrices.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <stdlib.h>
#include <math.h>

#include "integrals_2el_layer.h"
#include "integrals_2el.h"
#include "integrals_2el_boxed.h"
#include "integrals_2el_J.h"
#include "integrals_2el_K.h"
#include "utilities.h"
#include "output.h"
#include "memorymanag.h"


int 
compute_2e_matrix_exchange(const BasisInfoStruct & basisInfo,
			   const IntegralInfo & integralInfo,
			   const JK::ExchWeights & CAM_params,
			   ergo_real* K,
			   ergo_real* dens,
			   ergo_real threshold)
{
  int symmetryFlag = 0;
  JK::Params J_K_params;
  J_K_params.threshold_K = threshold;
  return compute_K_by_boxes_dense(basisInfo,
				  integralInfo,
				  CAM_params,
				  J_K_params,
				  K,
				  dens,
				  symmetryFlag);
}


int 
compute_2e_matrix_coulomb(const BasisInfoStruct & basisInfo,
			  const IntegralInfo & integralInfo,
			  ergo_real* J,
			  ergo_real* dens,
			  const JK::Params& J_K_params)
{
      do_output(LOG_CAT_INFO, LOG_AREA_INTEGRALS, "calling compute_J_by_boxes");
      if(compute_J_by_boxes(basisInfo,
			    integralInfo,
			    J_K_params,
			    J,
			    dens) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_INTEGRALS, "error in compute_J_by_boxes");
	  return -1;
	}
  return 0;
}

