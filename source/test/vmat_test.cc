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

/** @file vmat_test.cc Tests the potential energy matrix
    construction. The purpose of the test in its current form is
    mostly to verify compilation correctness. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>

#include "integrals_1el_potential.h"
#include "integrals_1el_single.h"
#include "integrals_general.h"
#include "matrix_typedefs.h"
#include "integral_matrix_wrappers.h"
#include "matrix_utilities.h"

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

static int 
test_S_V_comparison(const IntegralInfo & integralInfo)
{
  Molecule molecule1;
  Molecule molecule2;
  
  const ergo_real R = 22222.0;

  molecule1.addAtom(5.0, 0.0, 1.0, 0.0);
  molecule1.addAtom(9.0, 3.0, 0.0, 0.0);

  molecule2.addAtom(5.0, 0.0, 1.0, R  );
  molecule2.addAtom(9.0, 3.0, 0.0, R  );
  
  BasisInfoStruct basisInfo;
  if(basisInfo.addBasisfuncsForMolecule(molecule1,
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

  // Get overlap matrix
  std::vector<int> permutationHML, inversePermutationHML;
  mat::SizesAndBlocks sizeBlockInfo;
  preparePermutations(basisInfo, sizeBlockInfo,
		      permutationHML, inversePermutationHML);

  symmMatrix S;
  S.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
 
   if(compute_overlap_matrix_sparse(basisInfo, S, 
				   permutationHML) != 0)
    {
      puts("error in compute_overlap_matrix_sparse");
      return -1;
    }

  // Get V matrix
  ergo_real threshold = 1e-11;
  ergo_real boxSize = 3.3; // Use small box size to provoke errors.
  symmMatrix V;
  V.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  ergo_real nuclearRepulsionEnergy = 0;
  if(compute_V_sparse(basisInfo,
		      integralInfo,
		      molecule2,
		      threshold,
		      boxSize,
		      V,
		      permutationHML,
		      nuclearRepulsionEnergy) != 0)
    {
      puts("error in compute_V_sparse");
      return -1;
    }
  ergo_real nuclearRepulsionEnergyDiff = template_blas_fabs(nuclearRepulsionEnergy - molecule2.getNuclearRepulsionEnergyQuadratic());
  if(nuclearRepulsionEnergyDiff > 1e-9) {
    puts("Error: nuclearRepulsionEnergyDiff too large");
    return -1;
  }

  // Now V matrix should be approximately equal to -1 * (1/R) * (sum of charges) * S
  ergo_real sumOfCharges = 14.0;
  symmMatrix X(V);
  ergo_real factor = -1 * (1/R) * sumOfCharges;
  X *= (1/factor);
  
  ergo_real diff = symmMatrix::frob_diff(X, S);
  printf("diff = %22.11f\n", (double)diff);
  if(diff > 0.0003)
    {
      puts("error in V test: too large diff.");
      return -1;
    }

  puts("S vs V comparison test OK.");
  
  return 0;
}


static int
test_V_by_explicit_comparison(const IntegralInfo & integralInfo)
{
  Molecule molecule;
  
  // Put some atoms far away to make sure multipoles are used, but
  // still some of them close together so that significant overlaps
  // exist.
  molecule.addAtom(5.0, 0.0, 1.6, 0.4);
  molecule.addAtom(9.0, 3.7, 0.0, 0.3);
  molecule.addAtom(4.0, 1.7, 2.0, 0.8);
  molecule.addAtom(6.0, 10.0, 21.6, 88.4);
  molecule.addAtom(7.0, 13.7, 20.0, 88.3);
  molecule.addAtom(3.0, 11.7, 22.0, 88.8);
  molecule.addAtom(6.0, 10.0, 31.6, 77.4);
  molecule.addAtom(9.0, 13.7, 30.0, 77.3);
  molecule.addAtom(5.0, 11.7, 32.0, 77.8);
  molecule.addAtom(8.0, 44.7, 44.0, 44.8);  

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
  
  // Get V matrix
  std::vector<int> permutationHML, inversePermutationHML;
  mat::SizesAndBlocks sizeBlockInfo;
  preparePermutations(basisInfo, sizeBlockInfo,
		      permutationHML, inversePermutationHML);
  ergo_real threshold = 1e-12;
  ergo_real boxSize = 3.3; // Use small box size to provoke errors.
  symmMatrix V;
  V.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  ergo_real nuclearRepulsionEnergy_1 = 0;
  if(compute_V_sparse(basisInfo,
		      integralInfo,
		      molecule,
		      threshold,
		      boxSize,
		      V,
		      permutationHML,
		      nuclearRepulsionEnergy_1) != 0)
    {
      puts("error in compute_V_sparse");
      return -1;
    }

  {
    // Test also compute_V_sparse_hierarchical
    symmMatrix V2;
    V2.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
    ergo_real nuclearRepulsionEnergy_2 = 0;
    if(compute_V_sparse_hierarchical(basisInfo,
				     integralInfo,
				     molecule,
				     threshold,
				     boxSize,
				     V2,
				     permutationHML,
				     nuclearRepulsionEnergy_2) != 0)
      {
	puts("error in compute_V_sparse_hierarchical");
	return -1;
      }
    ergo_real diff = symmMatrix::frob_diff(V, V2);
    ergo_real tolerance = 2e-9;
    printf("compute_V_sparse_hierarchical comparison gave diff = %g\n", (double)diff);
    if(diff > tolerance) {
      printf("Error in test_V_by_explicit_comparison, compute_V_sparse_hierarchical comparison gave diff = %g which is greater than tolerance %g\n", (double)diff, (double)tolerance);
      return -1;
    }
    ergo_real nuclRepDiff = template_blas_fabs(nuclearRepulsionEnergy_2 - nuclearRepulsionEnergy_1);
    if(nuclRepDiff > tolerance) {
      printf("Error in test_V_by_explicit_comparison, compute_V_sparse_hierarchical comparison gave nuclRepDiff = %g which is greater than tolerance %g\n", (double)nuclRepDiff, (double)tolerance);
      return -1;
    }
  }

  // Convert V to full matrix format
  int n = basisInfo.noOfBasisFuncs;
  std::vector<ergo_real> V_full(n*n);
  V.fullMatrix(V_full, 
	       inversePermutationHML,
	       inversePermutationHML);

  // Check each element by explicit computation.
  ergo_real maxAbsDiff = 0;
  ergo_real sum_for_maxAbsDiff = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	// Compute list of simpleprimitives for product of basis funcs i and j.
	const int maxCount = 888;
	DistributionSpecStruct list[maxCount];
	int nPrims = get_product_simple_primitives(basisInfo, i, basisInfo, j, list, maxCount, 0);
	if(nPrims < 0)
	  return -1;
	ergo_real sum = 0;
	for(int m = 0; m < nPrims; m++)
	  for(int k = 0; k < molecule.getNoOfAtoms(); k++)
	    sum += -1 * do_1e_repulsion_integral_using_symb_info(list[m], molecule.getAtom(k).charge,
								 molecule.getAtom(k).coords,
								 integralInfo);
	ergo_real absDiff = template_blas_fabs(V_full[i*n+j] - sum);
	if(absDiff > maxAbsDiff) {
	  maxAbsDiff = absDiff;
	  sum_for_maxAbsDiff = sum;
	}
      }
  printf("maxAbsDiff = %22.15f  (sum_for_maxAbsDiff = %22.15f)\n", (double)maxAbsDiff, (double)sum_for_maxAbsDiff);
  ergo_real tolerance1 = 2e-9;
#ifdef PRECISION_SINGLE
  tolerance1 = 1e-2;
#endif
  if(maxAbsDiff > tolerance1)
    {
      puts("error in test_V_by_explicit_comparison: maxAbsDiff too large.");
      return -1;
    }

  // Also test the compute_h_core_matrix_full routine.
  std::vector<ergo_real> V_full_other(n*n);
  if(compute_V_matrix_full(basisInfo, integralInfo, molecule.getNoOfAtoms(),
			   molecule.getAtomListPtr(), threshold, &V_full_other[0]) != 0)
    {
      puts("error in compute_V_matrix_full");
      return -1;
    }
  ergo_real maxAbsDiff2 = 0;
  ergo_real value_for_maxAbsDiff2 = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      ergo_real absDiff = template_blas_fabs(V_full[i*n+j] - V_full_other[i*n+j]);
      if(absDiff > maxAbsDiff2) {
	maxAbsDiff2 = absDiff;
	value_for_maxAbsDiff2 = V_full[i*n+j];
      }
    }
  printf("maxAbsDiff2 = %22.15f  (value_for_maxAbsDiff2 = %22.15f)\n", (double)maxAbsDiff2, (double)value_for_maxAbsDiff2);
  ergo_real tolerance2 = 2e-9;
#ifdef PRECISION_SINGLE
  tolerance2 = 1e-2;
#endif
  if(maxAbsDiff2 > tolerance2)
    {
      puts("error in test_V_by_explicit_comparison: maxAbsDiff2 too large.");
      return -1;
    }

  puts("test_V_by_explicit_comparison OK.");
  return 0;
}




static int
test_V_by_explicit_comparison_tight(const IntegralInfo & integralInfo)
{
  Molecule molecule;

  // Put atoms close together so that basis funcs overlap with most
  // other basis funcs.
  molecule.addAtom(1.0, 0.0, 1.6, 0.4);
  molecule.addAtom(1.0, 3.7, 0.0, 0.3);
  molecule.addAtom(1.0, 1.7, 2.0, 0.8);
  molecule.addAtom(1.0, 0.1, 1.3, 0.4);
  molecule.addAtom(1.0, 3.8, 0.3, 0.3);
  molecule.addAtom(1.0, 1.8, 2.3, 0.8);
  molecule.addAtom(1.0, 0.3, 1.7, 0.4);
  molecule.addAtom(1.0, 3.5, 0.7, 0.3);
  molecule.addAtom(1.0, 1.2, 2.7, 0.8);
  molecule.addAtom(1.0, 0.3, 1.1, 0.4);
  molecule.addAtom(1.0, 3.6, 0.2, 0.3);
  molecule.addAtom(1.0, 1.9, 2.3, 0.8);
  molecule.addAtom(1.0, 0.3, 1.1, 1.1);
  molecule.addAtom(1.0, 3.6, 0.2, 1.2);
  molecule.addAtom(1.0, 1.9, 2.3, 1.4);

  BasisInfoStruct basisInfo;
  if(basisInfo.addBasisfuncsForMolecule(molecule,
					 "STO-3G",
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
  
  // Get V matrix
  std::vector<int> permutationHML, inversePermutationHML;
  mat::SizesAndBlocks sizeBlockInfo;
  preparePermutations(basisInfo, sizeBlockInfo,
		      permutationHML, inversePermutationHML);
  ergo_real threshold = 1e-12;
  ergo_real boxSize = 3.3; // Use small box size to provoke errors.
  symmMatrix V;
  V.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  ergo_real nuclearRepulsionEnergy_dummy = 0;
  if(compute_V_sparse(basisInfo,
		      integralInfo,
		      molecule,
		      threshold,
		      boxSize,
		      V,
		      permutationHML,
		      nuclearRepulsionEnergy_dummy) != 0)
    {
      puts("error in compute_V_sparse");
      return -1;
    }

  // Convert V to full matrix format
  int n = basisInfo.noOfBasisFuncs;
  std::vector<ergo_real> V_full(n*n);
  V.fullMatrix(V_full, inversePermutationHML, inversePermutationHML);

  // Check each element by explicit computation.
  ergo_real maxAbsDiff = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      {
	// Compute list of simpleprimitives for product of basis funcs i and j.
	const int maxCount = 888;
	DistributionSpecStruct list[maxCount];
	int nPrims = get_product_simple_primitives(basisInfo, i, basisInfo, j, list, maxCount, 0);
	if(nPrims < 0)
	  return -1;
	ergo_real sum = 0;
	for(int m = 0; m < nPrims; m++)
	  for(int k = 0; k < molecule.getNoOfAtoms(); k++)
	    sum += -1 * do_1e_repulsion_integral_using_symb_info(list[m], molecule.getAtom(k).charge,
								 molecule.getAtom(k).coords,
								 integralInfo);
	ergo_real absDiff = template_blas_fabs(V_full[i*n+j] - sum);
	if(absDiff > maxAbsDiff)
	  maxAbsDiff = absDiff;
      }
  ergo_real tolerance = 1e-11;
#ifdef PRECISION_SINGLE
  tolerance = 1e-5;
#endif
  if(maxAbsDiff > tolerance)
    {
      printf("error in V matrix test: maxAbsDiff = %22.15f\n", (double)maxAbsDiff);
      return -1;
    }

  return 0;
}



int main(int argc, char *argv[])
{
  IntegralInfo integralInfo(true);
  int errorCount = 0;
  if(test_S_V_comparison(integralInfo) != 0)
    {
      puts("error in test_S_V_comparison.");
      errorCount++;
    }

  if(test_V_by_explicit_comparison(integralInfo) != 0)
    {
      puts("error in test_V_by_explicit_comparison.");
      errorCount++;
    }

  int N = 100;
  if(argc == 2)
    N = atoi(argv[1]);
  for(int i = 0; i < N; i++)
    {
      if(test_V_by_explicit_comparison_tight(integralInfo) != 0)
	errorCount++;
    }
  printf("test_V_by_explicit_comparison_tight test repeated, %i times, errorCount = %i\n", N, errorCount);
  
  if(errorCount != 0)
    {
      puts("ERROR in V matrix tests.");
      return -1;
    }

  puts("V matrix tests OK.");
  return 0;
}
