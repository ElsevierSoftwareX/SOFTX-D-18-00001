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

/** @file scf_utils.h

    @brief Various utilities used by self-consistent field (SCF)
    code. For example, interface routines converting between HML
    matrix format used in main SCF code and elementwise format used by
    integral code.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef SCF_UTILS_HEADER
#define SCF_UTILS_HEADER

#include "molecule.h"
#include "basisinfo.h"
#include "integrals_2el.h"
#include "matrix_typedefs.h"
#include "grid_stream.h"
#include "SCF_statistics.h"


void output_sparsity(int n, const normalMatrix & M, const char* matrixName);
void output_sparsity_symm(int n, const symmMatrix & M, const char* matrixName);
void output_sparsity_triang(int n, const triangMatrix & M, const char* matrixName);

int
compute_h_core_matrix_sparse(const IntegralInfo& integralInfo,
			     const Molecule& molecule,
			     const Molecule& extraCharges,
			     ergo_real electric_field_x,
			     ergo_real electric_field_y,
			     ergo_real electric_field_z,
			     const BasisInfoStruct& basisInfo,
			     symmMatrix & H_core_Matrix_sparse,
			     ergo_real threshold_integrals_1el,
			     int noOfThreadsForV,
			     ergo_real boxSizeForVT,
			     ergo_real & result_nuclearRepulsionEnergy,
			     mat::SizesAndBlocks const & matrix_size_block_info,
			     std::vector<int> const & permutationHML,
			     int const create_dipole_mtx = 0,
			     std::vector<int> const * const inversePermutationHML = 0,
			     std::string const * const calculation_identifier = 0,
			     std::string const * const method_and_basis_set = 0);

int
compute_h_core_matrix_simple_dense(const IntegralInfo& integralInfo,
				   const Molecule& molecule,
				   const BasisInfoStruct& basisInfo,
				   symmMatrix & H_core_Matrix_sparse,
				   ergo_real threshold_integrals_1el,
				   int noOfThreadsForV,
				   mat::SizesAndBlocks const & matrix_size_block_info,
				   std::vector<int> const & permutationHML,
				   ergo_real & result_nuclearRepulsionEnergy);

int
get_gradient_for_given_mol_and_dens(const IntegralInfo& integralInfo,
				    const Molecule& molecule,
				    const BasisInfoStruct& basisInfo,
				    const symmMatrix & D,
				    ergo_real threshold_integrals_1el,
				    mat::SizesAndBlocks const & matrix_size_block_info,
				    std::vector<int> const & permutationHML,
				    ergo_real* result_gradient_list);

int save_symmetric_matrix(symmMatrix& A, 
                          const BasisInfoStruct & basisInfo,
                          const char *name,
			  std::vector<int> const & inversePermutationHML);

int 
add_disturbance_to_matrix(int n, 
			  symmMatrix & A,
			  ergo_real disturbance,
			  int specificElementCount,
			  const int* elementIndexVector,
			  std::vector<int> const & permutationHML);

int
get_simple_starting_guess_sparse(int n, 
				 int noOfElectrons, 
				 symmMatrix & densityMatrix);

int
write_diag_elements_to_file(int n, 
			    const symmMatrix & M, 
			    const char* fileName,
			    std::vector<int> const & permutationHML);

int
get_diag_matrix_from_file(int n, 
			  symmMatrix & M, 
			  const char* fileName,
			  std::vector<int> const & permutationHML);

int 
write_full_matrix(int n, 
		  const symmMatrix & M, 
		  const char* fileName,
		  std::vector<int> const & inversePermutationHML);

int
write_basis_func_coord_file(const BasisInfoStruct & basisInfo);

int
write_2el_integral_m_file(const BasisInfoStruct & basisInfo, const IntegralInfo & integralInfo);

int
get_2e_matrix_and_energy_sparse(const BasisInfoStruct & basisInfo,
				const Molecule& molecule,
				const IntegralInfo& integralInfo, 
				symmMatrix & twoelMatrix_sparse, 
				symmMatrix & densityMatrix_sparse,
				const JK::Params& J_K_params,
				const JK::ExchWeights & CAM_params,
				const Dft::GridParams& gridParams,
				int do_xc,
				ergo_real* energy_2el,
				int noOfElectrons,
				mat::SizesAndBlocks const & matrix_size_block_info,
				std::vector<int> const & permutationHML,
				std::vector<int> const & inversePermutationHML,
				int get_J_K_Fxc_matrices,
				symmMatrix & J_matrix,
				symmMatrix & K_matrix,
				symmMatrix & Fxc_matrix,
				SCF_statistics & stats);

int
get_2e_matrices_and_energy_sparse_unrestricted(const BasisInfoStruct & basisInfo, 
					       const Molecule& molecule,
					       const IntegralInfo& integralInfo, 
					       const JK::ExchWeights & CAM_params,
					       symmMatrix & twoelMatrix_sparse_alpha, 
					       symmMatrix & twoelMatrix_sparse_beta, 
					       symmMatrix & densityMatrix_sparse_alpha,
					       symmMatrix & densityMatrix_sparse_beta,
					       const JK::Params& J_K_params,
					       const Dft::GridParams& gridParams,
					       int do_xc,
					       ergo_real* energy_2el,
					       int noOfElectrons,
					       mat::SizesAndBlocks const & matrix_size_block_info,
					       std::vector<int> const & permutationHML,
					       std::vector<int> const & inversePermutationHML);

int
get_2e_matrices_and_energy_restricted_open(const BasisInfoStruct & basisInfo, 
					   const Molecule& molecule,
					   const IntegralInfo& integralInfo, 
					   const JK::ExchWeights & CAM_params,
					   symmMatrix & twoelMatrix_Fc, 
					   symmMatrix & twoelMatrix_Fo, 
					   symmMatrix & densityMatrix_sparse_alpha,
					   symmMatrix & densityMatrix_sparse_beta,
					   const JK::Params& J_K_params,
				           const Dft::GridParams& gridParams,
					   int do_xc,
					   ergo_real* energy_2el,
					   int noOfElectrons,
					   mat::SizesAndBlocks const & matrix_size_block_info,
					   std::vector<int> const & permutationHML,
					   std::vector<int> const & inversePermutationHML);

int
compute_FDSminusSDF_sparse(int n, 
			   symmMatrix & F_symm, 
			   symmMatrix & D_symm, 
			   symmMatrix & S_symm, 
			   normalMatrix & result, 
			   ergo_real sparse_threshold);

int 
determine_number_of_electrons_unrestricted(int noOfElectrons, 
					   int alpha_beta_diff, 
					   int* noOfElectrons_alpha, 
					   int* noOfElectrons_beta);

void 
get_hf_weight_and_cam_params(int use_dft, 
			     ergo_real* exch_param_alpha, 
			     ergo_real* exch_param_beta, 
			     ergo_real* exch_param_mu);

int 
determine_number_of_electrons_unrestricted(int noOfElectrons, 
					   int alpha_beta_diff, 
					   int* noOfElectrons_alpha, 
					   int* noOfElectrons_beta);

void
do_mulliken_atomic_charges(const symmMatrix & densityMatrix,
			   const symmMatrix & S_symm,
			   const BasisInfoStruct & basisInfo,
			   mat::SizesAndBlocks const & matrix_size_block_info,
			   std::vector<int> const & permutationHML,
			   std::vector<int> const & inversePermutationHML,
			   const Molecule& molecule);

void
do_mulliken_spin_densities(const symmMatrix & spinDensityMatrix,
			   const symmMatrix & S_symm,
			   const BasisInfoStruct & basisInfo,
			   mat::SizesAndBlocks const & matrix_size_block_info,
			   std::vector<int> const & permutationHML,
			   std::vector<int> const & inversePermutationHML,
			   const Molecule& molecule);

void get_exp_value_pos_operator(const BasisInfoStruct & basisInfo,
				const Molecule& molecule,
				const symmMatrix & densityMatrix,
				mat::SizesAndBlocks const & matrix_size_block_info,
				std::vector<int> const & permutationHML,
				std::vector<ergo_real> &mean,
				std::vector<ergo_real> &std);


void
do_density_images(const BasisInfoStruct & basisInfo,
		  const Molecule& molecule,
		  const ergo_real* densityMatrixFull_tot, 
		  const ergo_real* densityMatrixFull_spin,
		  double output_density_images_boxwidth, 
		  const std::string &filename_id = "");

void
do_acc_scan_J(const symmMatrix & D,
	      const IntegralInfo & integralInfo,
	      const BasisInfoStruct & basisInfo,
	      triangMatrix & invCholFactor,
	      bool doInvCholFactorTransformation,
	      const JK::Params & J_K_params,
	      mat::SizesAndBlocks const & matrix_size_block_info,
	      std::vector<int> const & permutationHML,
	      int nSteps,
	      ergo_real startThresh,
	      ergo_real stepFactor);

void
do_acc_scan_K(symmMatrix & D,
	      const IntegralInfo & integralInfo,
	      const BasisInfoStruct & basisInfo,
	      triangMatrix & invCholFactor,
	      bool doInvCholFactorTransformation,
	      const JK::ExchWeights & CAM_params,
	      const JK::Params & J_K_params,
	      mat::SizesAndBlocks const & matrix_size_block_info,
	      std::vector<int> const & permutationHML,
	      std::vector<int> const & inversePermutationHML,
	      int nSteps,
	      ergo_real startThresh,
	      ergo_real stepFactor);

void
do_acc_scan_Vxc(symmMatrix & D,
		const IntegralInfo & integralInfo,
		const BasisInfoStruct & basisInfo,
		const Molecule & molecule,
		const Dft::GridParams & gridParams,
		int noOfElectrons,
		triangMatrix & invCholFactor,
		bool doInvCholFactorTransformation,
		mat::SizesAndBlocks const & matrix_size_block_info,
		std::vector<int> const & permutationHML,
		std::vector<int> const & inversePermutationHML,
		int nSteps,
		ergo_real startThresh,
		ergo_real stepFactor);

void
create_mtx_files_with_different_orderings(const symmMatrix & A,
					  const std::string & calculation_identifier,
					  const std::string & method_and_basis_set,
					  const std::vector<int> & inversePermutationHML,
					  const BasisInfoStruct & basisInfo);

#endif
