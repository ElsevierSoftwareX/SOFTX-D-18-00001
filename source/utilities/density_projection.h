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

/** @file density_projection.h

    @brief Functionality for preparing a starting guess density matrix
    given a previous density matrix. The old density is read from
    file, and a projection between the basis sets is performed.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef DENSITY_PROJECTION
#define DENSITY_PROJECTION

#include "basisinfo.h"
#include "matrix_typedefs.h"
#include "GetDensFromFock.h"

/** load_density_and_project_full loads one or two density matrices (depending
    on value of noOfDensityMatrices) from the file specified by
    densityFileName. 
    @param densityFileName Name of file to load density matrices from
    @param noOfDensityMatrices Number of density matrices to load
    @param integralInfo static helper object for
    integral evaluation. 
    @param basisInfo the basis set that the density is
    to be expanded into (if the original density was expressed with
    help of other basis set, an apriopriate projection will be
    performed). 
    @param densityMatrixList must be already allocated and have
    proper dimension. 
    @param do_purification determines whether an additional
    purification is to be run after the projection. 
    @param noOfElectronsList is an one or two element array specyfying 
    the number of total electrons (one element) 
    or alpha and beta electrons (two elements).
    @param electronic_temperature Electronic temperature
*/
int load_density_and_project_full(const char *densityFileName,
				  int noOfDensityMatrices,
				  const IntegralInfo* integralInfo,
				  const BasisInfoStruct & basisInfo,
				  ergo_real** densityMatrixList,
				  int do_purification,
				  const int* noOfElectronsList,
				  ergo_real electronic_temperature);


/** load_density_and_project_sparse loads one or two density matrices (depending
    on value of noOfDensityMatrices) from the file specified by
    densityFileName. 

    The projection is done as follows:
    First, a matrix R is computed. 
    R is the overlap matrix between the two basis sets.
    Then RT * P * R is computed, where P is the starting guess density matrix
    read from file.
    To get a final projected density one could then multiply by S_inv from
    both sides, but to prepare for purification the matrix S*D*S is needed,
    so we skip multiplication by S_inv since it will anyway be cancelled out.

    @param DensFromFock Instance of GetDensFromFock class contatining all data for computing the density matrix.
    @param densityFileName Name of file to load density matrices from
    @param noOfDensityMatrices Number of density matrices to load
    @param integralInfo static helper object for
    integral evaluation. 
    @param basisInfo the basis set that the density is
    to be expanded into (if the original density was expressed with
    help of other basis set, an apriopriate projection will be
    performed). 
    @param S_symm Overlap matrix
    @param densityMatrixList pointers to one or two empty matrices that will 
    contain the result.
    Purification is always run after the projection. 
    @param noOfElectronsList is an one or two element array specyfying 
    the number of total electrons (one element) 
    or alpha and beta electrons (two elements).
    @param matrix_size_block_info Information about HML matrix block sizes etc.
    @param matrixPermutationVec Permutation vector used when calling matrix lib.
    @param sparse_threshold Threshold used when truncating matrices.

*/
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
				ergo_real sparse_threshold);




#endif /*  DENSITY_PROJECTION */
