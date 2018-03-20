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

/** @file density_projection.cc

    @brief Functionality for preparing a starting guess density matrix
    given a previous density matrix. The old density is read from
    file, and a projection between the basis sets is performed.

    @author: Elias Rudberg <em>responsible</em>
*/

/* 
   Projection of electron density from one basis set to another.
*/

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>

#include "density_description_file.h"
#include "density_projection.h"
#include "densfromf_full.h"
#include "integrals_general.h"
#include "operator_matrix.h"
#include "matrix_algebra.h"
#include "memorymanag.h"
#include "output.h"
#include "utilities.h"
#include "matrix_utilities.h"
#include "integral_matrix_wrappers.h"


int
load_density_and_project_full(const char *densityFileName,
			      int noOfDensityMatrices,
			      const IntegralInfo* integralInfo,
			      const BasisInfoStruct & basisInfo,
			      ergo_real** densityMatrixList,
			      int do_purification,
			      const int* noOfElectronsList,
			      ergo_real electronic_temperature)
{
  int n = basisInfo.noOfBasisFuncs;
  ergo_real* densityMatrixListForStartingGuess[2];
  BasisInfoStruct* basisInfoStartingGuess = NULL;
    
  if(ddf_load_density(densityFileName, noOfDensityMatrices,
                      *integralInfo, &basisInfoStartingGuess,
                      densityMatrixListForStartingGuess) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in ddf_load_density");
      return -1;
    }

  /* Projection part, could be a separate routine, too. */

  // get matrix R : overlap of main basis and startingguess basis
  int n_sg = basisInfoStartingGuess->noOfBasisFuncs;
  ergo_real* R = ergo_new(n*n_sg, ergo_real);
  do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "calling compute_overlap_matrix for R");
  if(compute_overlap_matrix(*basisInfoStartingGuess, basisInfo, R) != 0)
  {
      do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in compute_overlap_matrix for matrix R");
      return -1;
  }

  // get matrix S : main overlap matrix
  ergo_real* S = ergo_new(n*n, ergo_real);
  if(compute_overlap_matrix(basisInfo, basisInfo, S) != 0)
  {
      do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in compute_overlap_matrix for matrix S");
      return -1;
  }

  // get matrix Sinv : inverse of main overlap matrix
  //  ergo_real* Sinv = ergo_new(n*n, ergo_real);

  
  do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error: get_inverse_of_posdef_symm_matrix not implemented.");
  return -1;
  /*
  if(get_inverse_of_posdef_symm_matrix(n, S, Sinv) != 0)
  {
      do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in get_inverse_of_posdef_symm_matrix");
      return -1;
  }
  */
  /* ELIAS NOTE 2011-09-15: 
     Remainder of this routine removed since it could anyway not be used since get_inverse_of_posdef_symm_matrix was not implemented. */
}




int
load_density_and_project_sparse(GetDensFromFock &DensFromFock,
				const char *densityFileName,
				int noOfDensityMatrices,
				const IntegralInfo* integralInfo,
				const BasisInfoStruct & basisInfo,
				symmMatrix & S_symm,
				symmMatrix** densityMatrixList,
				const int* noOfElectronsList,
				mat::SizesAndBlocks matrix_size_block_info,
				std::vector<int> const & matrixPermutationVec,
				ergo_real sparse_threshold)
{
  Util::TimeMeter timeMeter;
  BasisInfoStruct* basisInfoStartingGuess = NULL;
  
  long nvaluesList[2];
  int* rowindList[2];
  int* colindList[2];
  ergo_real* valuesList[2];
  // Call ddf_load_density_sparse. It will allocate rowindList, colindList, valuesList,
  // which we will have to delete afterwards.
  int noOfDensitiesRead = 0;
  if(ddf_load_density_sparse(densityFileName,
			     *integralInfo, &basisInfoStartingGuess,
			     &noOfDensitiesRead,
			     rowindList,
			     colindList,
			     valuesList,
			     nvaluesList) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in ddf_load_density_sparse");
      return -1;
    }

  int n_sg = basisInfoStartingGuess->noOfBasisFuncs;

  // Get permutation object for starting guess basis set.
  int sparse_block_size = 20;
  int blockSizeFactor = 8;
  mat::SizesAndBlocks matrix_size_block_info_sg 
    = prepareMatrixSizesAndBlocks(n_sg,
                                  sparse_block_size,
                                  blockSizeFactor, 
                                  blockSizeFactor, 
                                  blockSizeFactor);
  std::vector<int> matrixPermutationVec_sg;
  getMatrixPermutation(*basisInfoStartingGuess,
                       sparse_block_size,
                       blockSizeFactor, blockSizeFactor, blockSizeFactor,
                       matrixPermutationVec_sg);
  
  // Now we have one or two matrices in memory stored as vectors. Convert to symmMatrix form.
  symmMatrix dens_list_1_sparse[2];
  int maxDensities = noOfDensityMatrices > noOfDensitiesRead
    ? noOfDensityMatrices : noOfDensitiesRead;
  
  for(int i = 0; i < maxDensities; i++)
    dens_list_1_sparse[i].resetSizesAndBlocks
      (matrix_size_block_info_sg, matrix_size_block_info_sg);

  for(int i = 0; i < noOfDensitiesRead; i++)
    {
      std::vector<int> rowIndTmp(rowindList[i], rowindList[i] + nvaluesList[i]);
      std::vector<int> colIndTmp(colindList[i], colindList[i] + nvaluesList[i]);
      std::vector<ergo_real> valuesTmp(valuesList[i], valuesList[i] + nvaluesList[i]);
      dens_list_1_sparse[i].assign_from_sparse(rowIndTmp, colIndTmp, valuesTmp,
                                               matrixPermutationVec_sg,
                                               matrixPermutationVec_sg);
      dens_list_1_sparse[i].eucl_thresh(sparse_threshold);
      delete []rowindList[i];
      delete []colindList[i];
      delete []valuesList[i];
    }  

  /* Handle conversions between restricted and unrestricted. */
  if(noOfDensitiesRead == 1 && noOfDensityMatrices == 2) {
    do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "converting restricted density into unrestricted");
    dens_list_1_sparse[1] = dens_list_1_sparse[0];
  } else if(noOfDensitiesRead == 2 && noOfDensityMatrices == 1) {
    do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "converting unrestricted density into restricted");
    dens_list_1_sparse[0] += dens_list_1_sparse[1];
  }

  // Reduce memory fragmentation by writing to file and reading in again.
  for(int i = 0; i < noOfDensityMatrices; i++)
    dens_list_1_sparse[i].writeToFile();
  for(int i = 0; i < noOfDensityMatrices; i++)
    dens_list_1_sparse[i].readFromFile();


  /* Projection part, could be a separate routine, too. */

  // get matrix R : overlap of main basis and startingguess basis
  output_current_memory_usage(LOG_AREA_UNDEFINED, "Before getting compute_R_matrix_sparse");

  normalMatrix R_sparse;
  R_sparse.resetSizesAndBlocks
    (matrix_size_block_info_sg, matrix_size_block_info);
  
  if(compute_R_matrix_sparse(basisInfo,
			     *basisInfoStartingGuess,
			     R_sparse,
			     sparse_threshold,
                             matrixPermutationVec,
                             matrixPermutationVec_sg) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in compute_R_matrix_sparse");
      return -1;
    }

  output_current_memory_usage(LOG_AREA_UNDEFINED, "After compute_R_matrix_sparse");

  // Now calculate RT * P * R for each density matrix
  // FIXME: implement BT * P * B in matrix library. B is normalMatrix, P is symmMatrix. Result is symm.
  for(int i = 0; i < noOfDensityMatrices; i++)
    {
      normalMatrix PR_sparse;
      PR_sparse.resetSizesAndBlocks
        (matrix_size_block_info_sg, matrix_size_block_info);

      do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "trying to create PR_sparse");

      Util::TimeMeter timeMeterPR;
      PR_sparse = (ergo_real)1.0 * dens_list_1_sparse[i] * R_sparse;
      timeMeterPR.print(LOG_AREA_UNDEFINED, "PR_sparse = P * R multiplication");

      PR_sparse.eucl_thresh(sparse_threshold);

      output_current_memory_usage(LOG_AREA_UNDEFINED, "After creating PR_sparse");

      dens_list_1_sparse[i].clear();

      normalMatrix RT;
      RT.resetSizesAndBlocks
        (matrix_size_block_info, matrix_size_block_info_sg);
      RT = transpose(R_sparse);
      
      normalMatrix RT_P_R;
      RT_P_R.resetSizesAndBlocks
        (matrix_size_block_info, matrix_size_block_info);
      
      RT_P_R = (ergo_real)1.0 * RT * PR_sparse;

      output_current_memory_usage(LOG_AREA_UNDEFINED, "After creating RT_P_R");
      RT.clear();
      PR_sparse.clear();

      RT_P_R.eucl_thresh(sparse_threshold);
      *densityMatrixList[i] = RT_P_R;
      densityMatrixList[i]->eucl_thresh(sparse_threshold);
    }

  // Now we do not need R_sparse any more.
  R_sparse.clear();

  // Projection done. Now do purification to force idempotency and correct trace.
  // densityMatrixList now contains RT * P * R for each density matrix.

  output_current_memory_usage(2, "While doing projection of starting guess (2).");

  for(int i = 0; i < noOfDensityMatrices; i++)
    {
      symmMatrix SDS;
      SDS.resetSizesAndBlocks
        (matrix_size_block_info, matrix_size_block_info);
      SDS = *densityMatrixList[i];
      SDS *= -1.0;

      output_current_memory_usage(LOG_AREA_UNDEFINED, "After creating matrix -SDS");
	  
      if(noOfElectronsList == NULL)
	return -1;

      symmMatrix F_ort_prev_dummy;
      F_ort_prev_dummy.resetSizesAndBlocks
        (matrix_size_block_info, matrix_size_block_info);

        int noOfOccupiedOrbs = noOfElectronsList[i];
      if(noOfDensityMatrices == 1)
	noOfOccupiedOrbs /= 2;
      do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, 
		"calling get_dens_from_fock_general for SDS to force idempotency and correct trace of starting guess");

      // The result of purification will be placed in densityMatrixList[i], but it is not used as input,
      // so we can clear it now to free up some memory.
      densityMatrixList[i]->clear();
      densityMatrixList[i]->writeToFile();

      F_ort_prev_dummy.writeToFile();
      SDS.writeToFile();

      DensFromFock.set_no_occupied_orbs(noOfOccupiedOrbs);

      // save old values
      int use_diag_on_error  = DensFromFock.get_use_diag_on_error();

      
      /*
	Needed for the old stopping criterion !!!

	The threshold for the eigenvalue error is multiplied by some
	value (usually making it lower).  Used just for the initial
	guess density matrix, since we want eigenvalues be closer to 0
	and 1.
       */
      ergo_real puri_eig_acc_factor_for_guess = DensFromFock.get_puri_eig_acc_factor_for_guess();
      ergo_real eigvalueErrorLimit_old = DensFromFock.get_eigvalueErrorLimit();
      DensFromFock.set_eigvalueErrorLimit(eigvalueErrorLimit_old*puri_eig_acc_factor_for_guess);


      DensFromFock.set_SCF_step(0);
      if(noOfDensityMatrices == 1)
	DensFromFock.set_generate_figures("");
      else
	{
	  if(i == 0)
	    DensFromFock.set_generate_figures("_alpha");
	  else if (i == 1)
	    DensFromFock.set_generate_figures("_beta");
	  else
	    DensFromFock.set_generate_figures("_i"); // just in case
	}


    // unset parameters which we do not need here
      DensFromFock.unset_use_diag_on_error();

      DensFromFock.clean_eigs_intervals();
      
      if(DensFromFock.get_dens_from_fock(SDS, 
					 *densityMatrixList[i],
					 F_ort_prev_dummy) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_UNDEFINED, "error in get_dens_from_fock.");
	  return -1;
	}
      
      // return old values
      if( use_diag_on_error == 1 )
	DensFromFock.set_use_diag_on_error();
      DensFromFock.unset_generate_figures();
      DensFromFock.set_eigvalueErrorLimit(eigvalueErrorLimit_old); // needed for the old stopping criterion
      
      
      output_current_memory_usage(LOG_AREA_UNDEFINED, "After get_dens_from_fock_general");
      do_output(LOG_CAT_INFO, LOG_AREA_UNDEFINED, "get_dens_from_fock_general finished OK.");
    } // END FOR i
  delete basisInfoStartingGuess;
  
  timeMeter.print(LOG_AREA_UNDEFINED, "load_density_and_project_sparse");
  
  return 0;
}


