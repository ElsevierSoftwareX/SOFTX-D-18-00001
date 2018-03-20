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

/** @file integral_matrix_wrappers.h

    @brief Wrapper routines for different parts of the integral code,
    including conversion of matrices from/to the hierarchical matrix
    library (HML) format.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef INTEGRAL_MATRIX_WRAPPERS_HEADER
#define INTEGRAL_MATRIX_WRAPPERS_HEADER

#include "basisinfo.h"
#include "matrix_typedefs.h"
#include "integrals_2el.h"


int
compute_V_sparse(const BasisInfoStruct& basisInfo,
		 const IntegralInfo& integralInfo,
		 const Molecule& molecule,
		 ergo_real threshold,
		 ergo_real boxSize,
		 symmMatrix & V,
		 std::vector<int> const & permutationHML,
		 ergo_real & result_nuclearRepulsionEnergy);

int
compute_V_sparse_hierarchical(const BasisInfoStruct& basisInfo,
			      const IntegralInfo& integralInfo,
			      const Molecule& molecule,
			      ergo_real threshold,
			      ergo_real boxSize,
			      symmMatrix & V,
			      std::vector<int> const & permutationHML,
			      ergo_real & result_nuclearRepulsionEnergy);


int
compute_gradient_of_nucl_and_trDV(const BasisInfoStruct& basisInfo,
				  const IntegralInfo& integralInfo,
				  const Molecule& molecule,
				  ergo_real threshold,
				  ergo_real boxSize,
				  const symmMatrix & densityMatrix_sparse,
				  std::vector<int> const & permutationHML,
				  ergo_real* result_gradient_list);


ergo_real
get_electron_nuclear_attraction_energy(const IntegralInfo& integralInfo,
				       const Molecule& molecule,
				       const BasisInfoStruct& basisInfo,
				       const symmMatrix & D,
				       ergo_real threshold_integrals_1el,
				       mat::SizesAndBlocks const & matrix_size_block_info,
				       std::vector<int> const & permutationHML);


int
compute_T_sparse_linear(const BasisInfoStruct& basisInfo,
			const IntegralInfo& integralInfo,
			ergo_real threshold,
			ergo_real boxSize,
			symmMatrix & T,
			std::vector<int> const & permutationHML);


int
compute_overlap_matrix_sparse(const BasisInfoStruct& basisInfo,
			      symmMatrix & S_symm,
			      std::vector<int> const & permutationHML);


/** compute_R_matrix_sparse computes the overlap matrix between two
    different basis sets.
*/
int
compute_R_matrix_sparse(const BasisInfoStruct & basisInfo_A,
			const BasisInfoStruct & basisInfo_B,
			normalMatrix & result_R,
			ergo_real sparse_threshold,
                        std::vector<int> const & matrixPermutationVec_A,
                        std::vector<int> const & matrixPermutationVec_B);


int
compute_operator_matrix_sparse_symm(const BasisInfoStruct& basisInfo,
				    int pow_x,
				    int pow_y,
				    int pow_z,
				    symmMatrix & A_symm,
				    std::vector<int> const & permutationHML);


int
compute_J_by_boxes_sparse(const BasisInfoStruct& basisInfo,
			  const IntegralInfo& integralInfo, 
			  const JK::Params& J_K_params,
			  symmMatrix & J,
 			  const symmMatrix & densityMatrix_sparse,
			  std::vector<int> const & permutationHML);


int
compute_K_by_boxes_sparse(const BasisInfoStruct& basisInfo,
			  const IntegralInfo& integralInfo, 
			  const JK::ExchWeights & CAM_params,
			  const JK::Params& J_K_params,
			  symmMatrix & K,
			  symmMatrix & densityMatrix_sparse,
			  std::vector<int> const & permutationHML,
			  std::vector<int> const & inversePermutationHML);

int
compute_K_by_boxes_sparse_nosymm(const BasisInfoStruct& basisInfo,
				 const IntegralInfo& integralInfo, 
				 const JK::ExchWeights & CAM_params,
				 const JK::Params& J_K_params,
				 normalMatrix & K,
				 normalMatrix & densityMatrix_sparse,
				 std::vector<int> const & permutationHML,
				 std::vector<int> const & inversePermutationHML);


#endif
