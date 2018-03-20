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

/** @file scf.cc

    @brief Code for classes containing various options related to
    self-consistent field (SCF) calculations.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "scf.h"
#include "basisinfo.h"
#include "molecule.h"
#include "output.h"
#include "utilities.h"
#include "matrix_utilities.h"
#include "AllocatorManager.h"



void SCF::MatOptions::prepare(const BasisInfoStruct& basisInfo)
{
  if(use_allocator_manager)
  {
    int bufferSize = sparse_matrix_block_size * sparse_matrix_block_size;
    mat::AllocatorManager<ergo_real>::instance().init(bufferSize, no_of_buffers_per_allocator);
  }
  
  size_block_info = prepareMatrixSizesAndBlocks(basisInfo.noOfBasisFuncs,
						sparse_matrix_block_size,
						sparse_matrix_block_factor_1,
						sparse_matrix_block_factor_2,
						sparse_matrix_block_factor_3);
  getMatrixPermutation(basisInfo,
		       sparse_matrix_block_size,
		       sparse_matrix_block_factor_1,
		       sparse_matrix_block_factor_2,
		       sparse_matrix_block_factor_3,
		       permutationHML,
		       inversePermutationHML);

  
  static const char *matOptionsInvalidInput =
    "MatOptions: invalid number of OpenMP threads";
  /* Set number of threads to use in matrix library. */
#ifdef _OPENMP
  if (threads < 1) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
	      "Error: Number of threads for matrix library less than one.");
    throw matOptionsInvalidInput;
  }
  if (parallelLevel < 1) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
	      "Error: Parallel level for matrix library less than one.");
    throw matOptionsInvalidInput;
  }
  
  mat::Params::setNProcs((unsigned int)threads);
  mat::Params::setMatrixParallelLevel((unsigned int)parallelLevel);
  
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
            "OpenMP used in matrix library. Number of threads set to %i.",
            mat::Params::getNProcs());
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
            "OpenMP used in matrix library. Parallel level set to    %i.",
            mat::Params::getMatrixParallelLevel());

#else
  if (threads != 1) {
    do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, 
	      "Error: Number of threads for matrix library is set to %i. "
              "Only one thread can be used when openmp is not used.",
              threads);
    throw matOptionsInvalidInput;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, 
	    "OpenMP not used in matrix library. Number of threads set to 1.");
#endif
}

