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

/** @file kmat_nosymm_test.cc Tests the sparse exchange matrix
    construction for non-symmetric density matrices. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include <vector>

#include "matrix_utilities.h"
#include "integral_matrix_wrappers.h"
#include "integrals_2el_explicit.h"


int main(int argc, char *argv[])
{
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  if ( !(env && (defThreads=atoi(env)) > 0) ) {
    defThreads = 1;
  }  
#ifdef _OPENMP
  mat::Params::setNProcs(defThreads);
  mat::Params::setMatrixParallelLevel(2);
  std::cout<<"OpenMP is used, number of threads set to "
	   <<mat::Params::getNProcs()<<". Matrix parallel level: "
	   <<mat::Params::getMatrixParallelLevel()<<"."<<std::endl;
#endif

  IntegralInfo integralInfo(true);
  BasisInfoStruct bis;

  Molecule m;

  int nx, ny, nz;
  if(getenv("RUN_BENCHMARK")) 
    {
      nx = 5;
      ny = 4;
      nz = 4;
    }
  else 
    {
      nx = 2;
      ny = 2;
      nz = 2;
    }

  const ergo_real space = 8.8;
  int atomCount = 0;
  for(int ix = 0; ix < nx; ix++)
    for(int iy = 0; iy < ny; iy++)
      for(int iz = 0; iz < nz; iz++)
	{
	  ergo_real x = ix*space + 0.4*std::cos((ix+iy+iz)*0.2+0.0)*space;
	  ergo_real y = iy*space + 0.4*std::cos((ix+iy+iz)*0.2+0.3)*space;
	  ergo_real z = iz*space + 0.4*std::cos((ix+iy+iz)*0.2+0.6)*space;
	  /* Use a mix of charges: H, C, Zn.
	     It is good to have some Zn there so we check also usage
	     of basis functions of f type. */
	  int charge = 1;
	  if(atomCount%3 == 0)
	    charge = 6;
	  //	  if(atomCount%9 == 0)
	  //	    charge = 30;
	  m.addAtom(charge, x, y, z);
	  atomCount++;
	}
	
  if(bis.addBasisfuncsForMolecule(m, ERGO_SPREFIX "/basis/3-21G",
                                   0, NULL, integralInfo, 0, 0, 0) != 0) {
    printf("bis.addBasisfuncsForMolecule failed.\n");
    return -1;
  }

  mat::SizesAndBlocks matrix_size_block_info =
    prepareMatrixSizesAndBlocks(bis.noOfBasisFuncs,
				20, 8, 8, 8);
  std::vector<int> permutationHML(bis.noOfBasisFuncs);
  std::vector<int> inversePermutationHML(bis.noOfBasisFuncs);
  getMatrixPermutation(bis, 20, 8, 8, 8, 
		       permutationHML,
		       inversePermutationHML);

  symmMatrix D_symm;
  D_symm.resetSizesAndBlocks(matrix_size_block_info,
			     matrix_size_block_info);
  symmMatrix K_symm;
  K_symm.resetSizesAndBlocks(matrix_size_block_info,
			     matrix_size_block_info);


  {
    /* Add values to density matrix. */
    int n = bis.noOfBasisFuncs;
    /* Put values on diagonal and next to diagonal. */
    int nvalues = n + 1 * (n-1);
    std::vector<int> rowind(nvalues);
    std::vector<int> colind(nvalues);
    std::vector<ergo_real> values(nvalues);
    int count = 0;
    // Put ones on diagonal.
    for(int i = 0; i < n; i++) {
      rowind[count] = i;
      colind[count] = i;
      values[count] = 1;
      count++;
    }
    // Put small random values next to diagonal.
    for(int i = 0; i < n-1; i++) {
      int j = i + 1;
      ergo_real value = 0.2 - 0.1 * ((double)rand() / RAND_MAX);
      rowind[count] = i;
      colind[count] = j;
      values[count] = value;
      count++;
    }
    D_symm.assign_from_sparse(rowind,
			      colind,
			      values,
			      permutationHML,
			      permutationHML);
  }

  JK::Params J_K_params;
  J_K_params. noOfThreads_K = defThreads;
  static const ergo_real EPS = mat::getMachineEpsilon<ergo_real>();
  J_K_params.threshold_K = template_blas_sqrt(EPS);

  JK::ExchWeights CAM_params_not_used;


  J_K_params.exchange_box_size = 4;
  if(compute_K_by_boxes_sparse(bis, integralInfo, CAM_params_not_used, J_K_params, K_symm, D_symm,
			       permutationHML,
			       inversePermutationHML) != 0)
    {
      printf("Error in compute_K_by_boxes_sparse\n");
      return -1;
    }

  //  ergo_real acc = template_blas_sqrt(EPS);
  //  ergo_real diffNorm = K_diff.eucl(acc);
  ergo_real requestedAcc;
  if(getenv("RUN_BENCHMARK")) 
    requestedAcc = J_K_params.threshold_K*25*nx*ny*nz;
  else
    requestedAcc = J_K_params.threshold_K*25*nx*ny*nz;

  // Also test result by explicit computation of all integrals.
  ergo_real maxabsdiff_explicitcheck = 0;
  {
    std::vector<int> rowind_D;
    std::vector<int> colind_D;
    std::vector<ergo_real> values_D;
    D_symm.get_all_values(rowind_D,
			  colind_D,
			  values_D,
			  inversePermutationHML,
			  inversePermutationHML);
    int nvalues_D = values_D.size();
    std::vector<int> rowind_K;
    std::vector<int> colind_K;
    std::vector<ergo_real> values_K;
    K_symm.get_all_values(rowind_K,
			  colind_K,
			  values_K,
			  inversePermutationHML,
			  inversePermutationHML);
    int nvalues_K = values_K.size();
    int n = bis.noOfBasisFuncs;
    // Create full matrix version of D.
    std::vector<ergo_real> D_full(n*n);
    for(int ii = 0; ii < n*n; ii++)
      D_full[ii] = 0;
    for(int ii = 0; ii < nvalues_D; ii++) {
      int i = rowind_D[ii];
      int j = colind_D[ii];
      D_full[i*n+j] = values_D[ii];
      D_full[j*n+i] = values_D[ii];
    }
    // Create full matrix version of K.
    std::vector<ergo_real> K_full(n*n);
    for(int ii = 0; ii < n*n; ii++)
      K_full[ii] = 0;
    for(int ii = 0; ii < nvalues_K; ii++) {
      int i = rowind_K[ii];
      int j = colind_K[ii];
      K_full[i*n+j] = values_K[ii];
      K_full[j*n+i] = values_K[ii];
    }
    // Now check K.
    for(int i = 0; i < n; i++)
      for(int j = i; j < n; j++) {
	// Compute element K_ij explicitly.
	ergo_real sum = 0;
	for(int k = 0; k < n; k++)
	  for(int l = 0; l < n; l++) {
	    ergo_real integral_iklj = do_2e_integral(i, k, l, j, bis, integralInfo);
	    ergo_real contrib = (-0.5) * D_full[k*n+l] * integral_iklj;
	    sum += contrib;
	  }
	ergo_real K_value_explicit = sum;
	ergo_real K_value_to_check = K_full[i*n+j];
	ergo_real absdiff = template_blas_fabs(K_value_explicit - K_value_to_check);
	if(absdiff > maxabsdiff_explicitcheck)
	  maxabsdiff_explicitcheck = absdiff;
      }
    if(maxabsdiff_explicitcheck > requestedAcc) {
      printf("Error in K test: diff too large in explicit comparison!\n");
      printf("maxabsdiff_explicitcheck  = %8.4g\n", (double)maxabsdiff_explicitcheck);
      printf("requestedAcc              = %8.4g\n", (double)requestedAcc);
      return -1;
    }
  }

  unlink("ergoscf.out");

  printf("K test OK, maxabsdiff_explicitcheck = %7.4g\n", (double)maxabsdiff_explicitcheck);
  return 0;
}
