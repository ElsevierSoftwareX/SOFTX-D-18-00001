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

/** @file SCF_general.cc

    @brief Class for self-consistent field (SCF) procedure; base class
    that can be used for both restricted and unrestricted cases.

    @author: Elias Rudberg <em>responsible</em>.
*/

#include <sstream>
#include <sys/types.h>
#include <unistd.h>
#include "SCF_general.h"
#include "output.h"
#include "scf_utils.h"
#include "matrix_utilities.h"
#include "utilities.h"
#include "integral_matrix_wrappers.h"
#include "machine_epsilon.h"
#include "units.h"
#include "SCF_statistics.h"
#include "AllocatorManager.h"


static ergo_real get_eucl_norm_try_different_acc(const symmMatrix & A, ergo_real & acc) {
  // First try with acc = sqrt(epsilon) and if that is too difficult try with lower accuracy.
  acc = template_blas_sqrt(mat::getMachineEpsilon<ergo_real>());
  int maxIter = 1000;
  ergo_real result_eucl_norm = -1;
  try {
    result_eucl_norm = A.eucl(acc, maxIter);
    return result_eucl_norm;
  }
  catch(...) {
    ergo_real accNew = template_blas_sqrt(acc);
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_eucl_norm_try_different_acc(): first attempt failed, now reducing accuracy from %g to %g", acc, accNew);
    acc = accNew;
  }
  try {
    result_eucl_norm = A.eucl(acc, maxIter);
    return result_eucl_norm;
  }
  catch(...) {
    ergo_real accNew = template_blas_sqrt(acc);
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_eucl_norm_try_different_acc(): second attempt failed, now reducing accuracy from %g to %g", acc, accNew);
    acc = accNew;
  }
  // Now call eucl without maxIter, let it use as many iterations as it needs to reach the reduced accuracy acc.
  return A.eucl(acc);
}


SCF_general::SCF_general(const Molecule& molecule_,
			 const Molecule& extraCharges_,
			 const BasisInfoStruct & basisInfo_, 
			 const IntegralInfo & integralInfo_,
			 const char* guessDmatFileName_,
			 const JK::Params& J_K_params_,
			 const Dft::GridParams& gridParams_,
			 const SCF::Options& scfoptsPtr,
			 const SCF::MatOptions& matOpts_,
			 ergo_real threshold_integrals_1el_input)
  : 
  molecule(molecule_),
  extraCharges(extraCharges_),
  basisInfo(basisInfo_),
  integralInfo(integralInfo_), 
  guessDmatFileName(guessDmatFileName_), // FIXME: copy this object properly
  J_K_params(J_K_params_), 
  gridParams(gridParams_), 
  scfopts(scfoptsPtr), 
  matOpts(matOpts_),
  threshold_integrals_1el(threshold_integrals_1el_input),
  DIIS(NULL),
  curr_cycle_stats(NULL)
{
  int n = basisInfo.noOfBasisFuncs;
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_general constructor, number of basis functions: %i", n);

  output_current_memory_usage(LOG_AREA_SCF, "beginning of SCF_general constructor");

  // Report info about host name, process ID etc.
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "************** Some general info here **********************");
  do_output_time(LOG_CAT_INFO, LOG_AREA_SCF, "VERSION: " VERSION "  time : ");
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "machine_epsilon = %9.3g", (double)get_machine_epsilon());
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "sizeof(ergo_real) = %i", sizeof(ergo_real));
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "sizeof(size_t)    = %i", sizeof(size_t));
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "sizeof(int)       = %i", sizeof(int));
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "sizeof(long)      = %i", sizeof(long));
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "sizeof(char*)     = %i", sizeof(char*));
  host_name_struct hostName;
  get_host_name(&hostName);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Host name:         '%s'", hostName.s);
  working_directory_struct workingDirectory;
  get_working_directory(&workingDirectory);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Working directory: '%s'", workingDirectory.s);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Process ID (PID): %10i", getpid());
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "************************************************************");
  ergo_real minDist, maxDist;
  molecule.getExtremeInternuclearDistancesQuadratic(minDist, maxDist);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Min internuclear distance: %12.5f a.u.  = %12.5f Angstrom", (double)minDist, (double)(minDist/UNIT_one_Angstrom));
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Max internuclear distance: %12.5f a.u.  = %12.5f Angstrom", (double)maxDist, (double)(maxDist/UNIT_one_Angstrom));

  S_symm.resetSizesAndBlocks(matOpts.size_block_info,
				    matOpts.size_block_info);
  if(compute_overlap_matrix_sparse(basisInfo, S_symm,
				   matOpts.permutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_overlap_matrix_sparse");
      throw "error in compute_overlap_matrix_sparse";
    }

  output_sparsity_symm(n, S_symm, "S_symm before trunc");  
  output_current_memory_usage(LOG_AREA_SCF, "after getting overlap matrix");
  {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "truncating S using threshold value %6.2g", (double)scfopts.sparse_threshold_for_S);
    double nnz_before_trunc_pc  = (double)S_symm.nnz()  * 100 / ((double)n*n); 
    ergo_real truncError = S_symm.eucl_thresh(scfopts.sparse_threshold_for_S);
    double nnz_after_trunc_pc  = (double)S_symm.nnz()  * 100 / ((double)n*n); 
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "Truncated S (eucl), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
	      scfopts.sparse_threshold_for_S, (double)truncError, nnz_before_trunc_pc, nnz_after_trunc_pc);    
  }


  if ( scfopts.create_mtx_file_S == 1 ) {
    // Write overlap matrix in matrix market format
    std::stringstream ss_id;
    ss_id << scfopts.calculation_identifier << " - overlap matrix";
    write_matrix_in_matrix_market_format( S_symm, matOpts.inversePermutationHML, "S_matrix", 
					  ss_id.str(), scfopts.method_and_basis_set );
  }
  if( scfopts.create_mtx_files_S_and_quit == 1 ) {
    // Write overlap matrix in matrix market format using different basis function orderings, and then quit.
    create_mtx_files_with_different_orderings(S_symm,
					      scfopts.calculation_identifier,
					      scfopts.method_and_basis_set,
					      matOpts.inversePermutationHML,
					      basisInfo);
    throw "Breaking because create_mtx_files_S_and_quit was set.";
  }
  if ( scfopts.create_basis_func_coord_file == 1 ) {
    write_basis_func_coord_file(basisInfo);
  }
  if ( scfopts.create_2el_integral_m_file == 1 ) 
    write_2el_integral_m_file(basisInfo, integralInfo);

  {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Computing Euclidean norm of overlap matrix (just for fun)...");
    Util::TimeMeter timeMeterEuclS;
    ergo_real acc = -1; // value set by call below
    ergo_real S_symm_euclnorm = get_eucl_norm_try_different_acc(S_symm, acc);
    timeMeterEuclS.print(LOG_AREA_SCF, "S_symm.eucl()");
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Euclidean norm of overlap matrix = %22.11f  (acc %7.3g)",
	      (double)S_symm_euclnorm, (double)acc);
  }

  invCholFactor.resetSizesAndBlocks(matOpts.size_block_info,
					   matOpts.size_block_info);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Calling invCholFactor.inch with matOpts.threshold_inch = %g", 
	    (double)matOpts.threshold_inch);
  Util::TimeMeter timeMeterInch;
  invCholFactor.inch(S_symm, matOpts.threshold_inch, mat::right);
  timeMeterInch.print(LOG_AREA_SCF, "invCholFactor.inch");

  output_sparsity_triang(n, invCholFactor, "invCholFactor before truncation");
  {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "truncating Z using threshold value %6.2g", (double)scfopts.sparse_threshold_for_Z);
    double nnz_before_trunc_pc  = (double)invCholFactor.nnz()  * 100 / ((double)n*n); 
    ergo_real truncError = invCholFactor.eucl_thresh(scfopts.sparse_threshold_for_Z);
    double nnz_after_trunc_pc  = (double)invCholFactor.nnz()  * 100 / ((double)n*n); 
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "Truncated Z (eucl), selected threshold = %10.6g, returned error = %10.6g, nnz before = %3.4f %%, nnz after = %3.4f %%",
	      scfopts.sparse_threshold_for_Z, (double)truncError, nnz_before_trunc_pc, nnz_after_trunc_pc);    
  }
  output_sparsity_triang(n, invCholFactor, "invCholFactor after truncation");
  {
    ergo_real invCholFactor_frobnorm = invCholFactor.frob();
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Frobenius norm of invCholFactor after truncation = %22.11f", 
	      (double)invCholFactor_frobnorm);
    Util::TimeMeter timeMeterEucl;
    invCholFactor_euclnorm = invCholFactor.eucl( template_blas_sqrt(mat::getMachineEpsilon<ergo_real>()) );
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Euclidean norm of invCholFactor after truncation = %22.11f", 
	      (double)invCholFactor_euclnorm);
    timeMeterEucl.print(LOG_AREA_SCF, "invCholFactor.eucl()");
  }
  output_current_memory_usage(LOG_AREA_SCF, "after getting invCholFactor");

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "testing invCholFactor by computing ZT*S*Z");
  {
    Util::TimeMeter timeMeterZZTS;
    normalMatrix noOver(S_symm);
    normalMatrix inchcopy(invCholFactor);
    normalMatrix tmpSZ, ID;
    normalMatrix the_real_identity;
    tmpSZ.resetSizesAndBlocks(matOpts.size_block_info,
			      matOpts.size_block_info);
    ID.resetSizesAndBlocks(matOpts.size_block_info,
                           matOpts.size_block_info);
    the_real_identity.resetSizesAndBlocks(matOpts.size_block_info,
                                          matOpts.size_block_info);
    tmpSZ = noOver * inchcopy;
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "Truncating tmp matrix S*Z using eucl_thresh() with threshold value %6.3g", 
	      (double)matOpts.sparse_threshold);
    tmpSZ.eucl_thresh(matOpts.sparse_threshold);
    ID =  transpose(inchcopy) * tmpSZ;
    the_real_identity = 1;
    ergo_real frobenius_error = normalMatrix::frob_diff(ID, the_real_identity);
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "frobenius_error for ZT*S*Z is %10g", (double)frobenius_error);
    timeMeterZZTS.print(LOG_AREA_SCF, "testing invCholFactor by computing ZT*S*Z");
  }
  

  S_symm.writeToFile();
  if(scfopts.write_overlap_matrix) {
    if(save_symmetric_matrix(S_symm, basisInfo, "overlap.bin",
			     matOpts.inversePermutationHML) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF,
                "error in ddf_writeShellListAndDensityMatricesToFile");
      throw "error in ddf_writeShellListAndDensityMatricesToFile";
    }
  }
  
  invCholFactor.writeToFile();
  output_current_memory_usage(LOG_AREA_SCF, "after writing invCholFactor and S_symm to file");

  H_core_Matrix.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  ergo_real nuclearRepulsionEnergyTmp = 0;
  if(scfopts.skip_H_core == 1)
    {
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "NOTE: skip_H_core parameter set, will skip construction of H_core matrix!");
      do_output(LOG_CAT_WARNING, LOG_AREA_SCF, "NOTE: skip_H_core parameter set, results will be bogus!");      
    }
  else
    {
      if(scfopts.use_simple_dense_H_core == 1) {
	if(extraCharges.getNoOfAtoms() != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error: (extraCharges.noOfAtoms != 0) not implemented for use_simple_dense_H_core case.");
	  throw "error: (extraCharges.noOfAtoms != 0) not implemented for use_simple_dense_H_core case.";
	}
	if(scfopts.electric_field.v[0] != 0 || scfopts.electric_field.v[1] != 0 || scfopts.electric_field.v[2] != 0) {
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error: electric field != 0 not implemented for use_simple_dense_H_core case.");
	  throw "error: electric field != 0 not implemented for use_simple_dense_H_core case.";
	}
	if(compute_h_core_matrix_simple_dense(integralInfo,
					      molecule,
					      basisInfo,
					      H_core_Matrix,
					      threshold_integrals_1el,
					      scfopts.no_of_threads_for_V,
					      matOpts.size_block_info,
					      matOpts.permutationHML,
					      nuclearRepulsionEnergyTmp) != 0)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_h_core_matrix_simple_dense");
	    throw "error in compute_h_core_matrix_simple_dense";
	  }
      }
      else {
	if(compute_h_core_matrix_sparse(integralInfo, 
					molecule, 
					extraCharges,
					scfopts.electric_field.v[0],
					scfopts.electric_field.v[1],
					scfopts.electric_field.v[2],
					basisInfo, 
					H_core_Matrix,
					threshold_integrals_1el,
					scfopts.no_of_threads_for_V,
					scfopts.box_size_for_V_and_T,
					nuclearRepulsionEnergyTmp,
					matOpts.size_block_info,
					matOpts.permutationHML,
					scfopts.create_mtx_files_dipole,
					&matOpts.inversePermutationHML,
					&scfopts.calculation_identifier,
					&scfopts.method_and_basis_set) != 0)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_h_core_matrix_sparse");
	    throw "error in compute_h_core_matrix_sparse";
	  }
      }
    }
  output_sparsity_symm(n, H_core_Matrix, "H_core_Matrix before trunc");
  {
    ergo_real subspaceThr = 0.001 * scfopts.purification_subspace_err_limit;    
    ergo_real threshold_Hcore = subspaceThr * scfopts.gap_expected_lower_bound / (1+subspaceThr);
    double nnz_before_trunc_pc  = (double)H_core_Matrix.nnz()  * 100 / ((double)n*n); 
    invCholFactor.readFromFile();
    ergo_real truncError = H_core_Matrix.eucl_thresh( threshold_Hcore, &invCholFactor );
    invCholFactor.writeToFile();
    double nnz_after_trunc_pc  = (double)H_core_Matrix.nnz()  * 100 / ((double)n*n); 
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	      "Truncated H_core_Matrix (eucl with Z), selected threshold = %10.6g, returned error = %10.6g, "
	      "nnz before = %3.4f %%, nnz after = %3.4f %%",
	      threshold_Hcore, truncError, nnz_before_trunc_pc, nnz_after_trunc_pc);    
  }
  output_sparsity_symm(n, H_core_Matrix, "H_core_Matrix after trunc");
  
  output_current_memory_usage(LOG_AREA_SCF, "after getting H_core_Matrix");

  ergo_real nuclearElectricFieldEnergyTmp = molecule.getNuclearElectricFieldEnergy(scfopts.electric_field);
  nuclearEnergy = nuclearRepulsionEnergyTmp + nuclearElectricFieldEnergyTmp;

  if ( scfopts.create_mtx_file_H_core == 1 ) {
    // Write H_core matrix in matrix market format
    std::stringstream ss_id;
    ss_id << scfopts.calculation_identifier << " - H_core matrix";
    write_matrix_in_matrix_market_format( H_core_Matrix, matOpts.inversePermutationHML, "H_core_matrix", 
					  ss_id.str(), scfopts.method_and_basis_set );
  }

  H_core_Matrix.writeToFile();
  output_current_memory_usage(LOG_AREA_SCF, "after writing H_core_Matrix to file");

  noOfElectrons = molecule.getNumberOfElectrons();  

  get_hf_weight_and_cam_params(scfopts.use_dft, &CAM_params.alpha,
                               &CAM_params.beta, &CAM_params.mu);

  CAM_params.computeRangeSeparatedExchange =
    CAM_params.beta != ergo_real(0.0);

  energy_2el = 0;
  energy = 0;
  energy_2el_core = 0; // only used when "core density matrix" is used
  energy_2el_valence = 0; // only used when "core density matrix" is used
  energy_of_valence = 0; // only used when "core density matrix" is used
  energy_reference = 0; // only used when "core density matrix" is used
  electronicEntropyTerm = 0;


  /*
    Set parameters in the GetDensFromFock class.
    All parameters must be set up, otherwise you will get exception.
   */
  
  // common parameters 
  DensFromFock.set_general_params(n,
	             		  matOpts.size_block_info);

#ifdef USE_CHUNKS_AND_TASKS
  DensFromFock.set_cht_matrix_params(scfopts.cht_leavesSizeMax, scfopts.cht_blocksize);
#endif

  DensFromFock.set_truncationNormPurification(scfopts.purification_truncation_norm);
  DensFromFock.set_stopCriterionNormPurification(scfopts.purification_stop_crit_norm);
  DensFromFock.set_invCholFactor(invCholFactor, invCholFactor_euclnorm);
  DensFromFock.set_gap_expected_lower_bound(scfopts.gap_expected_lower_bound);
  DensFromFock.set_purification_maxmul(scfopts.purification_maxmul);


  if( scfopts.purification_create_m_files > 0 )
    DensFromFock.set_purification_create_m_files();
  else
    DensFromFock.unset_purification_create_m_files();



  if( scfopts.output_homo_and_lumo_eigenvectors > 0 )
    DensFromFock.set_output_homo_and_lumo_eigenvectors();
  else
    DensFromFock.unset_output_homo_and_lumo_eigenvectors();



  /*
    do not plot figures from the purification by default
    set flag before calling function for purification
    unset after calling function for purification
    otherwise, it will be output from all its calls
  */
  DensFromFock.unset_generate_figures();


  if( scfopts.purification_ignore_failure > 0 )
    DensFromFock.set_purification_ignore_failure();
  else
    DensFromFock.unset_purification_ignore_failure();



  if( scfopts.purification_use_rand_perturbation_for_alleigsint > 0 )
    DensFromFock.set_purification_use_rand_perturbation_for_alleigsint();
  else
    DensFromFock.unset_purification_use_rand_perturbation_for_alleigsint();


  // purification or something else?
  if(scfopts.use_diagonalization == 0) // then use purification
    DensFromFock.set_use_purification();
  else
    DensFromFock.unset_use_purification();

  // purification settings
  if(DensFromFock.get_use_purification() > 0)
    {
      DensFromFock.set_purification_limits(scfopts.purification_subspace_err_limit,
					   scfopts.purification_eigvalue_err_limit, 
					   scfopts.puri_eig_acc_factor_for_guess);

      if(scfopts.purification_with_acceleration > 0)
	DensFromFock.set_use_acceleration();
      else
	DensFromFock.unset_use_acceleration();

      if(scfopts.use_new_stopping_criterion > 0)
	DensFromFock.set_use_new_stopping_criterion();
      else
	DensFromFock.unset_use_new_stopping_criterion();
    }


  if(scfopts.use_diag_on_error > 0)
    DensFromFock.set_use_diag_on_error();
  else
    DensFromFock.unset_use_diag_on_error();

  if(scfopts.use_diag_on_error_guess > 0)
    DensFromFock.set_use_diag_on_error_guess();
  else
    DensFromFock.unset_use_diag_on_error_guess();




  // diagonalization settings
  if(scfopts.use_diagonalization > 0)
    DensFromFock.set_use_diagonalization();
  else
    DensFromFock.unset_use_diagonalization();
  
  // may be needed if scfopts.use_diag_on_error is set
  DensFromFock.set_diagonalization_params(scfopts.electronic_temperature,
					  S_symm);

  if(scfopts.use_diagonalization > 0 && scfopts.store_all_eigenvalues_to_file > 0)
    DensFromFock.set_store_all_eigenvalues_to_file();
  else
    DensFromFock.unset_store_all_eigenvalues_to_file();


  if(scfopts.save_permuted_F_matrix_in_bin > 0)
    DensFromFock.set_save_permuted_F_matrix_in_bin();
  else
    DensFromFock.unset_save_permuted_F_matrix_in_bin();
}


SCF_general::~SCF_general()
{
  delete curr_cycle_stats;
}


ergo_real SCF_general::GetEuclideanNormOfMatrix(const symmMatrix & A)
{
  ergo_real acc = template_blas_sqrt(mat::getMachineEpsilon<ergo_real>());
  return A.eucl(acc);
}


void SCF_general::get_overlap_matrix(symmMatrix & S)
{
  S_symm.readFromFile();
  S = S_symm;
  S_symm.writeToFile();
}


void SCF_general::get_invCholFactor_matrix(triangMatrix & invCholFactor_)
{
  invCholFactor.readFromFile();
  invCholFactor_ = invCholFactor;
  invCholFactor.writeToFile();
}


void SCF_general::get_H_core_matrix(symmMatrix & H_core)
{
  H_core_Matrix.readFromFile();
  H_core = H_core_Matrix;
  H_core_Matrix.writeToFile();
}


void SCF_general::get_energy(ergo_real & E, ergo_real & E_nuclear)
{
  E = energy;
  E_nuclear = nuclearEnergy;
}


void SCF_general::do_SCF_iterations()
{
  Util::TimeMeter timeMeterTot;
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_general::do_SCF_iterations");

  // Initialize DIIS
  if(DIIS == NULL)
    throw "ERROR: (DIIS == NULL)";
  if(DIIS->Initialize(scfopts.max_no_of_diis_matrices) != 0)
    throw "Error in DIIS->Initialize";

  initialize_matrices();

  if(J_K_params.threshold_J <= 0 || J_K_params.threshold_K <= 0)
    throw "Error in SCF_general::do_SCF_iterations: (J_K_params.threshold_J <= 0 || J_K_params.threshold_K <= 0).";

  // Check that parameters are reasonable, even number of electrons for restricted etc.
  check_params();

  // set up starting guess
  get_starting_guess_density();

  if(scfopts.spin_flip_atom_count > 0)
    do_spin_flip(scfopts.spin_flip_atom_count);

  if(scfopts.starting_guess_disturbance > 0)
    add_random_disturbance_to_starting_guess();

  if(scfopts.write_guess_density_only == 1) {
    write_density_to_file();
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "throwing exception after write_density_to_file.");
    throw "exiting after scfopts.write_guess_density_only";
  }

  if(scfopts.output_density_images_only == 1) {
    output_density_images();
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "throwing exception after output_density_images.");
    throw "exiting after scfopts.output_density_images_only";
  }


  // Use at least two iterations if disimple_starting_guess used.
  int min_number_of_iterations = scfopts.min_number_of_iterations;
  if(scfopts.use_simple_starting_guess && min_number_of_iterations < 2)
    min_number_of_iterations = 2;

  ergo_real stepLength = scfopts.step_length_start;
  ergo_real best_energy_so_far = 0;
  int step = 0;
  SCF_step = 0;
  int noOfFailuresInARow = 0;
  int restartCount = 0;
  int noOfImprovementsInARow = 0;
  int diisUsedInLastIteration = 0;

  curr_subspace_diff = 0;

  initialize_homo_lumo_limits();

  // The different Fockmatrix objects are empty at this point,
  // But we write them to file because they are supposed to be on file 
  // in the beginning of each SCF cycle.
  write_matrices_to_file();


  output_current_memory_usage(LOG_AREA_SCF, "before main SCF loop");

  Util::TimeMeter timeMeterScfMainLoop;
  
  // main SCF loop
  while(1)
    {
      curr_cycle_stats = new SCF_statistics;
      curr_cycle_stats->start_timer("scf_cycle_time");
      curr_cycle_stats->add_value("no_of_basis_func", basisInfo.noOfBasisFuncs);
      curr_cycle_stats->add_value("sparse_matrix_block_size", matOpts.sparse_matrix_block_size);
      step++;
      SCF_step++;
      curr_cycle_stats->add_value("scf_cycle", step);
      char infoString[888];
      sprintf(infoString, "Beginning of SCF cycle %i: ", step);
      do_output_time(LOG_CAT_INFO, LOG_AREA_SCF, infoString);

      output_current_memory_usage(LOG_AREA_SCF, infoString);

      Util::TimeMeter timeMeterStep;

      save_density_as_prevdens();

      get_2e_part_and_energy();

      //if(scfopts.use_artificial_subspace_disturbances == 1)
      //disturb_fock_matrix(curr_subspace_diff * scfopts.subspace_factor_fock);

      // Now we have created Fock matrix (or matrices), and written to file. Memory usage should be the same as before.
      output_current_memory_usage(LOG_AREA_SCF, "After get_2e_part_and_energy");

      output_sparsity_S_F_D(*curr_cycle_stats);

      double virtMem = 0, resMem = 0, virtPeakMem = 0;
      get_memory_usage_by_procfile(&virtMem, &resMem, &virtPeakMem);
      curr_cycle_stats->add_value("peak_virt_mem_usage_GB", virtPeakMem);
      

      calculate_energy();

      
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "nuclearEnergy, energy_2el, energy = %f, %f, %f", 
		(double)nuclearEnergy, (double)energy_2el, (double)energy);
      // ELIAS NOTE 2014-01-01: FIXME: consider adding another output message here showing the "electronic energy" so that the difference between "electronic energy" and "total energy" becomes more clear.
      if( step > 2)
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Energy %3i = %22.11f  ( diff %22.15f )", step, (double)energy, (double)(energy-best_energy_so_far));
      else
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Energy %3i = %22.11f", step, (double)energy);

      if(scfopts.compute_core_density == 1 && step > 1) {
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Valence-only energy %3i  = %22.11f", step, (double)energy_of_valence);
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Reference energy    %3i  = %22.11f", step, (double)energy_reference);
      }

      output_current_memory_usage(LOG_AREA_SCF, "After computing energy");

      curr_cycle_stats->start_timer("get_FDSminusSDF");
      Util::TimeMeter timeMeterFDSminusSDF;
      get_FDSminusSDF();
      timeMeterFDSminusSDF.print(LOG_AREA_SCF, "get_FDSminusSDF");
      output_current_memory_usage(LOG_AREA_SCF, "After computing FDS-SDF");
      curr_cycle_stats->stop_timer("get_FDSminusSDF");

      Util::TimeMeter timeMeterGerErrorMeasure;
      get_error_measure();
      timeMeterGerErrorMeasure.print(LOG_AREA_SCF, "get_error_measure");
      output_current_memory_usage(LOG_AREA_SCF, "After computing error measure");

      // Check if converged
      if(scfopts.use_artificial_subspace_disturbances == 1 && step > 1 && curr_subspace_diff < 1e-6)
	{
	  do_output(LOG_CAT_RESULTS, LOG_AREA_SCF, "CONVERGED due to curr_subspace_diff after %3i iterations.", step);
	  do_output(LOG_CAT_RESULTS, LOG_AREA_SCF, "FINAL ENERGY: %22.11f", (double)energy);
	  report_final_results();
	  break;
	}
      if(errorMeasure < scfopts.convergence_threshold && step >= min_number_of_iterations)
	{
	  // Write Fock matrix in matrix market format in the last SCF cycle
	  if ( scfopts.create_mtx_files_F == 1 )
	    create_mtx_files_F( step );

	  do_output(LOG_CAT_RESULTS, LOG_AREA_SCF, "CONVERGED after %3i iterations.", step);
	  do_output(LOG_CAT_RESULTS, LOG_AREA_SCF, "FINAL ENERGY: %22.11f", (double)energy);
	  report_final_results();
	  break;
	}

      // Check if max number of iterations reached
      if(scfopts.max_number_of_iterations > 0 && step >= scfopts.max_number_of_iterations)
	{
	  // Write Fock matrix in matrix market format in the last SCF cycle
	  if ( scfopts.create_mtx_files_F == 1 )
	    create_mtx_files_F( step );

	  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Maximum number of SCF iterations reached. Breaking SCF procedure.");
	  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Gave up after %i iterations.", step);
	  break;
	}

      

      if(step == 1)
        {
          // first time
	  add_to_DIIS_list();
        }
      else
        {
          // not first time
          if((energy < best_energy_so_far) && 
	     (step > noOfImprovementsInARow+3) && 
	     (noOfImprovementsInARow < scfopts.no_of_impr_req_for_diis || 
	      errorMeasure > scfopts.error_maxabs_for_diis) && 
	     !(diisUsedInLastIteration == 1 && noOfImprovementsInARow > 1) &&
	     !(scfopts.use_diis_always == 1))
            {
              noOfImprovementsInARow++;
              do_output(LOG_CAT_INFO, LOG_AREA_SCF, "the energy got better, "
			"but we do not dare to try DIIS yet, noOfImprovementsInARow = %i", 
			noOfImprovementsInARow);
              best_energy_so_far = energy;
              noOfFailuresInARow = 0;

	      update_best_fock_so_far();
	      
              do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
			"mixing best Fock matrix so far with next one, stepLength = %10.6f", 
			(double)stepLength);
	      combine_old_fock_matrices(stepLength);
	      	      
              diisUsedInLastIteration = 0;
              stepLength *= 1.1;
              do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
			"increased steplength by factor 1.1, stepLength = %10.6f", (double)stepLength);
            }
          else if(step <= scfopts.no_of_careful_first_scf_steps)
	    {
	      // "careful" option chosen: in this case we do not use
	      // DIIS for the early steps. Can be useful when the
	      // starting guess is very bad.
              do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Not considering DIIS now because 'careful' option "
			"chosen, no_of_careful_first_scf_steps = %2d", scfopts.no_of_careful_first_scf_steps);
	      /* Elias note 2010-05-12: added || step == 2 in this if
		 statement condition to handle the case when
		 no_of_careful_first_scf_steps is used and the energy
		 increased in step 2, which heppened for the Umeda
		 protein molecule.  */
	      if(energy < best_energy_so_far || step == 2) {
		best_energy_so_far = energy;
		update_best_fock_so_far();
		stepLength *= 1.1;
		do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
			  "increased steplength by factor 1.1, stepLength = %10.6f", (double)stepLength);
	      }
	      else {
		// Energy got worse.
		noOfImprovementsInARow = 0;
		noOfFailuresInARow++;
		do_output(LOG_CAT_INFO, LOG_AREA_SCF, "restarting DIIS because energy did not improve.");
		clear_diis_list();
		ergo_real newStepLength = stepLength * 0.5;
		do_output(LOG_CAT_INFO, LOG_AREA_SCF, "changing stepLength from %10.6f to %10.6f", 
			  (double)stepLength, (double)newStepLength);
		stepLength = newStepLength;
	      }
              do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
			"mixing best Fock matrix so far with next one, stepLength = %10.6f", 
			(double)stepLength);
	      combine_old_fock_matrices(stepLength);
              diisUsedInLastIteration = 0;
	    }
          else if(step == 2 || energy < best_energy_so_far || scfopts.use_diis_always == 1)
            {
              if(step > 2 && energy < best_energy_so_far)
                noOfImprovementsInARow++;
              // energy got better
              noOfFailuresInARow = 0;
              best_energy_so_far = energy;
      
	      update_best_fock_so_far();

              diisUsedInLastIteration = 1;

	      add_to_DIIS_list();
	      
	      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "using DIIS to get combined Fock matrix, "
			"number of iters used for DIIS: %2i", DIIS->GetNoOfIters());
	      
	      use_diis_to_get_new_fock_matrix();
            }
          else
            {
              // energy got worse
	      if(scfopts.break_on_energy_increase == 1)
		{
		  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Energy increased. Breaking SCF.");
		  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Stopped SCF after %i iterations.", step);
		  break;
		}
              if(stepLength < scfopts.step_length_giveup)
                {
                  if(restartCount > scfopts.max_restart_count)
                    {
                      do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
				"The energy does not seem to get any lower. We give up!");
                      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Gave up after %i iterations.", step);
                      break;
                    }
                  restartCount++;
                  stepLength = scfopts.step_length_start;
                  noOfFailuresInARow = 0;
                  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "reset stepLength to %10.6f, restartCount = %i", 
			    (double)scfopts.step_length_start, restartCount);
                }
              noOfImprovementsInARow = 0;
              noOfFailuresInARow++;
              do_output(LOG_CAT_INFO, LOG_AREA_SCF, "restarting DIIS because energy did not improve.");
	      clear_diis_list();
	      
              if(diisUsedInLastIteration == 0)
		{
		  ergo_real newStepLength = stepLength * 0.5;
		  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "changing stepLength from %10.6f to %10.6f", 
			    (double)stepLength, (double)newStepLength);
		  stepLength = newStepLength;
		}
              do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
			"mixing best Fock matrix so far with next one, stepLength = %10.6f", 
			(double)stepLength);
	      combine_old_fock_matrices(stepLength);

              diisUsedInLastIteration = 0;
            }
        }

      output_current_memory_usage(LOG_AREA_SCF, "After creating lin comb F");

      // Free memory used by Err_sparse
      Util::TimeMeter timeMeterClearErrorMatrices;
      clear_error_matrices();
      timeMeterClearErrorMatrices.print(LOG_AREA_SCF, "clear_error_matrices");
      output_current_memory_usage(LOG_AREA_SCF, "After clear_error_matrices");

      // FIXME: Compare F and Fprev here to get info about gap of
      // F? Such info will be needed as input to get_dens_from_fock?

      Util::TimeMeter timeMeterSaveFockAsFprev;
      save_current_fock_as_fprev();
      timeMeterSaveFockAsFprev.print(LOG_AREA_SCF, "save_current_fock_as_fprev");

      curr_cycle_stats->start_timer("get_new_density_matrix");
      Util::TimeMeter timeMeterGetNewDensityMatrix;
      get_new_density_matrix();
      timeMeterGetNewDensityMatrix.print(LOG_AREA_SCF, "get_new_density_matrix");
      curr_cycle_stats->stop_timer("get_new_density_matrix");
      
      // At this point a new density matrix has just been computed, so the electronic entropy term has also been computed in the nonzero-temperature case.
      if(scfopts.electronic_temperature > 0) {
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Electronic temperature: %22.11f a.u. = %22.2f Kelvin", 
		  (double)scfopts.electronic_temperature, (double)scfopts.electronic_temperature/UNIT_one_Kelvin);
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Electronic entropy term: %22.11f, energy with entropy term included: %22.11f", 
		  (double)electronicEntropyTerm, (double)(energy + electronicEntropyTerm));
      }

      if(step > 1)
	report_density_difference();
      
      if(scfopts.use_artificial_subspace_disturbances == 1 && step > 2)
	disturb_dens_matrix_exact(curr_subspace_diff * scfopts.subspace_factor_dens);
      //disturb_dens_matrix(curr_subspace_diff * scfopts.subspace_factor_dens);

      if(scfopts.use_artificial_subspace_disturbances == 1)
	update_subspace_diff();

      if(scfopts.output_density_at_every_step == 1) {
	Util::TimeMeter timeMeterWriteDensityToFile;
	output_current_memory_usage(LOG_AREA_SCF, "before write_density_to_file()");
	write_density_to_file();
	output_current_memory_usage(LOG_AREA_SCF, "after  write_density_to_file()");
	timeMeterWriteDensityToFile.print(LOG_AREA_SCF, "write_density_to_file");
      }

      if ( scfopts.create_mtx_files_F == 1 )
	// Write Fock matrix in matrix market format
	create_mtx_files_F( step );
      if ( scfopts.create_mtx_files_D == 1 )
	// Write Fock matrix in matrix market format
	create_mtx_files_D( step );

      if ( scfopts.output_homo_and_lumo_eigenvectors ) {
      	// Write homo and lumo eigenvectors to file
      	create_eigenvectors_files();
      	create_gabedit_file();
      }
            
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF cycle %3i finished.", step);
      char tempString[88];
      sprintf(tempString, "SCF cycle %3i", step);
      timeMeterStep.print(LOG_AREA_SCF, tempString);
      curr_cycle_stats->stop_timer("scf_cycle_time");
      std::stringstream ss;
      ss << "scf_cycle_" << step;
      if(scfopts.output_statistics_mfiles)
	curr_cycle_stats->output_mfile( ss.str() );
      delete curr_cycle_stats;
      curr_cycle_stats = NULL; // since curr_cycle_stats is deleted in destructor, we need to set to null here to avoid double-delete if we exit loop in some unusual way.
    } // END WHILE main SCF loop

  timeMeterScfMainLoop.print(LOG_AREA_SCF, "Main SCF loop");
  output_current_memory_usage(LOG_AREA_SCF, "after main SCF loop");
  std::string allocStatsStr = mat::AllocatorManager<ergo_real>::instance().getStatistics();
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "After main SCF loop: %s", allocStatsStr.c_str());

  if(scfopts.output_density_at_every_step == 1 && step == 1) {
    Util::TimeMeter timeMeterWriteDensityToFile;
    output_current_memory_usage(LOG_AREA_SCF, "before write_density_to_file()");
    write_density_to_file();
    output_current_memory_usage(LOG_AREA_SCF, "after  write_density_to_file()");
    timeMeterWriteDensityToFile.print(LOG_AREA_SCF, "write_density_to_file");
  }

  compute_dipole_moment();

  if(scfopts.output_mulliken_pop == 1)
    do_mulliken_pop_stuff();

  if(scfopts.compute_gradient_fixeddens == 1)
    compute_gradient_fixeddens();

  if(scfopts.save_full_matrices_for_matlab == 1)
    save_full_matrices_for_matlab();
  
  if(scfopts.save_final_potential == 1)
    save_final_potential();
  
  if(scfopts.output_expected_values_pos_operator == 1)
    output_expected_values_pos_operator();

  if(scfopts.output_density_images == 1)
    output_density_images();
  
  if(scfopts.write_diag_dens_to_file == 1)
    write_diag_dens_to_file();

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "SCF_general::do_SCF_iterations finished.");
  timeMeterTot.print(LOG_AREA_SCF, "SCF_general::do_SCF_iterations");
}
