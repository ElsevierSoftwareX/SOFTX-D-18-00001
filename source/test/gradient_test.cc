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

/** @file gradient_test.cc

    @brief Tests the gradient computation.

    @author Elias Rudberg <em>responsible</em>
*/

#include "integrals_1el_potential.h"
#include "integrals_1el_single.h"
#include "integrals_general.h"
#include "matrix_typedefs.h"
#include "integral_matrix_wrappers.h"
#include "matrix_utilities.h"
#include "output.h"

static void
preparePermutations(const BasisInfoStruct& basisInfo,
		    mat::SizesAndBlocks& sizeBlockInfo, 
		    std::vector<int>& permutation,
		    std::vector<int>& inversePermutation)
{
  static const int sparseMatrixBlockSize = 16, sparseMatrixBlockFactor = 4;
  sizeBlockInfo =
    prepareMatrixSizesAndBlocks(basisInfo.noOfBasisFuncs,
				sparseMatrixBlockSize,
				sparseMatrixBlockFactor,
				sparseMatrixBlockFactor,
				sparseMatrixBlockFactor);
  getMatrixPermutation(basisInfo,
		       sparseMatrixBlockSize,
		       sparseMatrixBlockFactor,
		       sparseMatrixBlockFactor,
		       sparseMatrixBlockFactor,
		       permutation,
		       inversePermutation);
}

static ergo_real get_nucl_energy_for_given_mol_and_dens(const IntegralInfo& integralInfo,
							const Molecule& molecule,
							const BasisInfoStruct& basisInfo,
							const symmMatrix & D,
							ergo_real threshold_integrals_1el,
							mat::SizesAndBlocks const & matrix_size_block_info,
							std::vector<int> const & permutationHML) {
  ergo_real nuclearRepulsionEnergy = molecule.getNuclearRepulsionEnergyQuadratic();
  ergo_real elecNuclEnergy = get_electron_nuclear_attraction_energy(integralInfo,
								    molecule,
								    basisInfo,
								    D,
								    threshold_integrals_1el,
								    matrix_size_block_info,
								    permutationHML);
  return nuclearRepulsionEnergy + elecNuclEnergy;
}

static int get_gradient_using_finite_differences(const IntegralInfo& integralInfo,
						 const BasisInfoStruct& basisInfo,
						 const Molecule & molecule, 
						 const symmMatrix & densityMatrix,
						 ergo_real threshold_integrals_1el,
						 const mat::SizesAndBlocks & matrix_size_block_info,
						 const std::vector<int> & permutationHML,
						 ergo_real* resultGradient) {
  int nAtoms = molecule.getNoOfAtoms();
  for(int i = 0; i < nAtoms; i++) {
    for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
      const ergo_real h = 1e-3;
      Molecule moleculeTmp = molecule;
      Atom atomTmp = molecule.getAtom(i);
      atomTmp.coords[coordIdx] += h;
      moleculeTmp.replaceAtom(i, atomTmp);
      ergo_real E1 = get_nucl_energy_for_given_mol_and_dens(integralInfo,
							    moleculeTmp,
							    basisInfo,
							    densityMatrix,
							    threshold_integrals_1el,
							    matrix_size_block_info,
							    permutationHML);
      moleculeTmp = molecule;
      atomTmp = molecule.getAtom(i);
      atomTmp.coords[coordIdx] -= h;
      moleculeTmp.replaceAtom(i, atomTmp);
      ergo_real E2 = get_nucl_energy_for_given_mol_and_dens(integralInfo,
							    moleculeTmp,
							    basisInfo,
							    densityMatrix,
							    threshold_integrals_1el,
							    matrix_size_block_info,
							    permutationHML);
      ergo_real gradientComponent = (E1 - E2) / (2 * h);
      resultGradient[i*3+coordIdx] = gradientComponent;
    } // END FOR coordIdx
  } // END FOR i
  return 0;
}

static int get_gradient_using_explicit_integrals(const IntegralInfo& integralInfo,
						 const BasisInfoStruct& basisInfo,
						 const Molecule & molecule, 
						 const symmMatrix & densityMatrix,
						 ergo_real threshold_integrals_1el,
						 const mat::SizesAndBlocks & matrix_size_block_info,
						 const std::vector<int> & permutationHML,
						 const std::vector<int> & inversePermutationHML,
						 ergo_real* resultGradient) {  
  int nAtoms = molecule.getNoOfAtoms();
  memset(resultGradient, 0, nAtoms*3*sizeof(ergo_real));
  int n = basisInfo.noOfBasisFuncs;
  std::vector<ergo_real> densityMatrix_full(n*n);
  densityMatrix.fullMatrix(densityMatrix_full, inversePermutationHML, inversePermutationHML);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      // Consider gradient contributions related to V matrix element (i,j)
      // Compute list of simpleprimitives for product of basis funcs i and j.
      const int maxCount = 888;
      DistributionSpecStruct list[maxCount];
      int nPrims = get_product_simple_primitives(basisInfo, i, basisInfo, j, list, maxCount, 0);
      if(nPrims < 0)
	return -1;
      for(int m = 0; m < nPrims; m++)
	for(int k = 0; k < nAtoms; k++) {
	  std::vector<ergo_real> integralValues = do_1e_repulsion_integral_derivatives_using_symb_info(&list[m], 
												       molecule.getAtom(k).charge,
												       molecule.getAtom(k).coords,
												       integralInfo);
	  for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
	    ergo_real integralValue = integralValues[coordIdx];
	    resultGradient[k*3+coordIdx] += -1 * integralValue * densityMatrix_full[i*n+j];
	  }
	}
    }
  molecule.getNuclearRepulsionEnergyGradientContribQuadratic(resultGradient);
  return 0;
}

static int
test_gradient_by_explicit_comparison(const IntegralInfo & integralInfo)
{
  Molecule molecule;
  
  // Put some atoms far away to make sure multipoles are used, but
  // still some of them close together so that significant overlaps
  // exist.
  molecule.addAtom(5.0, 0.2, 1.5, 0.5);
  molecule.addAtom(9.0, 3.6, 0.1, 0.2);
  molecule.addAtom(4.0, 1.6, 2.1, 0.9);
  molecule.addAtom(6.0, 10.1, 21.4, 85.3);
  molecule.addAtom(7.0, 13.6, 20.2, 85.1);
  molecule.addAtom(3.0, 11.2, 22.4, 85.6);
  molecule.addAtom(6.0, 10.4, 31.3, 74.7);
  molecule.addAtom(9.0, 13.5, 30.2, 74.4);
#if 0
  molecule.addAtom(7.0, 16.6, 15.2, 85.1);
  molecule.addAtom(3.0, 15.2, 18.4, 85.6);
  molecule.addAtom(6.0, 14.4, 26.3, 74.7);
  molecule.addAtom(9.0, 13.5, 28.2, 74.1);
  molecule.addAtom(6.0, 4.4, 26.3, 34.7);
  molecule.addAtom(9.0, 3.5, 28.2, 34.1);
  molecule.addAtom(6.0, 4.4, 26.3, 36.9);
  molecule.addAtom(9.0, 3.5, 27.2, 36.5);
#endif

  BasisInfoStruct basisInfo;
  if(basisInfo.addBasisfuncsForMolecule(molecule,
					 "6-31Gss",
					 0,
					 NULL,
					 integralInfo,
					 0,
					 1,
					 0) != 0)
    {
      puts("error in basisInfo.addBasisfuncsGorMolecule");
      return -1;
    }
  
  std::vector<int> permutationHML, inversePermutationHML;
  mat::SizesAndBlocks sizeBlockInfo;
  preparePermutations(basisInfo, sizeBlockInfo,
		      permutationHML, inversePermutationHML);
  ergo_real threshold = 1e-12;
  ergo_real boxSize = 3.3; // Use small box size to provoke errors.
  symmMatrix V;
  V.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);

  symmMatrix D;
  D.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  {
    /* Add values to density matrix diagonal and one step 
       next to diagonal. */
    const int nvalues1 = basisInfo.noOfBasisFuncs;
    std::vector<int> idxrow(nvalues1);
    std::vector<int> idxcol(nvalues1);
    std::vector<ergo_real> values(nvalues1);
    for(int i=0; i<nvalues1; i++) {
      idxrow[i] = i;
      idxcol[i] = i;
      values[i] = 1.0;
    }
    D.add_values(idxrow, idxcol, values, permutationHML, permutationHML);
    const int nvalues2 = basisInfo.noOfBasisFuncs-1;
    for(int i=0; i<nvalues2; i++) {
      idxrow[i] = i;
      idxcol[i] = i+1;
      values[i] = 0.3;
    }
    idxrow.resize(nvalues2); 
    idxcol.resize(nvalues2); 
    values.resize(nvalues2); 
    D.add_values(idxrow, idxcol, values, permutationHML, permutationHML);
  }  

  int nAtoms = molecule.getNoOfAtoms();

  // Compute gradient in 3 ways:
  // - 1: Using finite differences
  // - 2: Using explicit integrals
  // - 3: By calling the real linear-scaling gradient computation routine
  std::vector<ergo_real> gradient1(nAtoms*3);
  std::vector<ergo_real> gradient2(nAtoms*3);
  std::vector<ergo_real> gradient3(nAtoms*3);

  // - 1: Compute gradient using finite differences
  if(get_gradient_using_finite_differences(integralInfo,
					   basisInfo,
					   molecule, 
					   D,
					   threshold,
					   sizeBlockInfo,
					   permutationHML,
					   &gradient1[0]) != 0) {
    printf("Error in get_gradient_using_finite_differences.\n");
    throw "Error in get_gradient_using_finite_differences";
  }

  // - 2: Compute gradient using explicit integrals
  if(get_gradient_using_explicit_integrals(integralInfo,
					   basisInfo,
					   molecule, 
					   D,
					   threshold,
					   sizeBlockInfo,
					   permutationHML,
					   inversePermutationHML,
					   &gradient2[0]) != 0) {
    printf("Error in get_gradient_using_explicit_integrals.\n");
    throw "Error in get_gradient_using_explicit_integrals";
  }

  // - 3: Compute gradient using by calling the real linear-scaling gradient computation routine
  if(compute_gradient_of_nucl_and_trDV(basisInfo,
				       integralInfo, 
				       molecule,
				       threshold,
				       boxSize,
				       D,
				       permutationHML,
				       &gradient3[0]) != 0) {
    printf("Error in compute_gradient_of_nucl_and_trDV.\n");
    throw "Error in compute_gradient_of_nucl_and_trDV";
  }

  ergo_real maxabsdiff_1_2 = 0;
  ergo_real maxabsdiff_1_3 = 0;
  ergo_real maxabsdiff_2_3 = 0;
  for(int i = 0; i < nAtoms; i++) 
    for(int coordIdx = 0; coordIdx < 3; coordIdx++) {
      int idx = i*3+coordIdx;
      ergo_real absdiff_1_2 = template_blas_fabs(gradient1[idx] - gradient2[idx]);
      ergo_real absdiff_1_3 = template_blas_fabs(gradient1[idx] - gradient3[idx]);
      ergo_real absdiff_2_3 = template_blas_fabs(gradient2[idx] - gradient3[idx]);
      if(absdiff_1_2 > maxabsdiff_1_2)
	maxabsdiff_1_2 = absdiff_1_2;
      if(absdiff_1_3 > maxabsdiff_1_3)
	maxabsdiff_1_3 = absdiff_1_3;
      if(absdiff_2_3 > maxabsdiff_2_3)
	maxabsdiff_2_3 = absdiff_2_3;
    }  
  printf("maxabsdiff_1_2 = %22.11f\n", (double)maxabsdiff_1_2);
  printf("maxabsdiff_1_3 = %22.11f\n", (double)maxabsdiff_1_3);
  printf("maxabsdiff_2_3 = %22.11f\n", (double)maxabsdiff_2_3);
  if(maxabsdiff_1_2 > 2e-5) {
    printf("Error: maxabsdiff_1_2 too large.\n");
    return -1;
  }
  if(maxabsdiff_2_3 > 5e-9) {
    printf("Error: maxabsdiff_2_3 too large.\n");
    return -1;
  }
  printf("OK!\n");
  puts("test_gradient_by_explicit_comparison OK.");
  return 0;
}



int main(int argc, char *argv[])
{
#ifdef PRECISION_SINGLE
  printf("gradient_test: skipping this test because PRECISION_SINGLE is defined.\n");
  return 0;
#endif

  IntegralInfo integralInfo(true);
  int errorCount = 0;

  // enable_output(); // Do this if you want the ergoscf.out file, to see timings etc.

  if(test_gradient_by_explicit_comparison(integralInfo) != 0)
    {
      puts("error in test_gradient_by_explicit_comparison.");
      errorCount++;
    }

  if(errorCount != 0)
    {
      puts("ERROR in gradient tests.");
      return -1;
    }

  puts("Gradient tests OK.");
  return 0;
}
