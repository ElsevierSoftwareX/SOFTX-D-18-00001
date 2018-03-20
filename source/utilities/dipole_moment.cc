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

/** @file dipole_moment.cc

    @brief Functionality for computing the dipole moment of a molecule
    for a given density matrix.

    @author: Elias Rudberg <em>responsible</em>
*/

#include "dipole_moment.h"
#include "integral_matrix_wrappers.h"
#include "output.h"
#include "matrix_utilities.h"

static ergo_real
compute_dipole_moment_onecoord(const symmMatrix & densityMatrix,
			       const BasisInfoStruct & basisInfo,
			       mat::SizesAndBlocks const & matrix_size_block_info,
			       std::vector<int> const & permutationHML,
			       const Molecule& molecule,
			       int coordIdx)
{
  int ix = 0;
  int iy = 0;
  int iz = 0;
  switch(coordIdx) {
  case 0: ix = 1; break;
  case 1: iy = 1; break;
  case 2: iz = 1; break;
  default: throw "Error in compute_dipole_moment_onecoord.";
  }
  symmMatrix opMatrix;
  opMatrix.resetSizesAndBlocks(matrix_size_block_info, matrix_size_block_info);
  if(compute_operator_matrix_sparse_symm(basisInfo, ix, iy, iz, opMatrix, permutationHML) != 0)
    throw "Error in compute_operator_matrix_sparse_symm";
  ergo_real density_contrib = symmMatrix::trace_ab(densityMatrix, opMatrix);
  // Now compute contribution from nuclei.
  ergo_real nuclear_contrib = 0;
  for(int i = 0; i < molecule.getNoOfAtoms(); i++)
    nuclear_contrib += molecule.getAtom(i).charge * molecule.getAtom(i).coords[coordIdx];
  ergo_real dipole_moment_oneCoord = nuclear_contrib - density_contrib;
  return dipole_moment_oneCoord;
}
  
void 
get_dipole_moment(const symmMatrix & densityMatrix,
		  const BasisInfoStruct & basisInfo,
		  mat::SizesAndBlocks const & matrix_size_block_info,
		  std::vector<int> const & permutationHML,
		  const Molecule& molecule,
		  int logArea,
		  const char* label)
{
  ergo_real dipole_moment_x = compute_dipole_moment_onecoord(densityMatrix, basisInfo, matrix_size_block_info, permutationHML, molecule, 0);
  ergo_real dipole_moment_y = compute_dipole_moment_onecoord(densityMatrix, basisInfo, matrix_size_block_info, permutationHML, molecule, 1);
  ergo_real dipole_moment_z = compute_dipole_moment_onecoord(densityMatrix, basisInfo, matrix_size_block_info, permutationHML, molecule, 2);
  do_output(LOG_CAT_INFO, logArea, "%s: dipole_moment_x [atomic units] = %16.8f", label, (double)dipole_moment_x);
  do_output(LOG_CAT_INFO, logArea, "%s: dipole_moment_y [atomic units] = %16.8f", label, (double)dipole_moment_y);
  do_output(LOG_CAT_INFO, logArea, "%s: dipole_moment_z [atomic units] = %16.8f", label, (double)dipole_moment_z);
}

void
get_dipole_moment_fullmat(int n,
			  const ergo_real* densityMatrix,
			  const BasisInfoStruct & basisInfo,
			  const Molecule& molecule,
			  int logArea,
			  const char* label)
{
  mat::SizesAndBlocks size_block_info;
  std::vector<int> permutationHML(n);
  size_block_info = prepareMatrixSizesAndBlocks(basisInfo.noOfBasisFuncs, 1000, 10, 10, 10);
  for(int i = 0; i < n; i++)
    permutationHML[i] = i;
  symmMatrix D;
  D.resetSizesAndBlocks(size_block_info, size_block_info);
  std::vector<ergo_real> vec(n*n);
  for(int i = 0; i < n*n; i++)
    vec[i] = densityMatrix[i];
  D.assignFromFull(vec, permutationHML, permutationHML);
  return get_dipole_moment(D, basisInfo, size_block_info, permutationHML, molecule, logArea, label);
}

