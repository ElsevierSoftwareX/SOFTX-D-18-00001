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

/** @file jmat_test.cc Tests the sparse Coulomb matrix construction. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include <vector>

#include "matrix_utilities.h"
#include "integral_matrix_wrappers.h"
#include "output.h"


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

  //  enable_output(); // Do this if you want the ergoscf.out file, to see error messages, timings etc.
  
  IntegralInfo integralInfo(true);
  BasisInfoStruct bis;

  Molecule m;

  int nx, ny, nz;
  if(getenv("RUN_BENCHMARK")) 
    {
      nx = 5;
      ny = 2;
      nz = 4;
    }
  else 
    {
      nx = 2;
      ny = 2;
      nz = 2;
    }

  const ergo_real space = 8.0;
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
	  if(atomCount%9 == 0)
	    charge = 30;
	  m.addAtom(charge, x, y, z);
	  atomCount++;
	}
	
  if(bis.addBasisfuncsForMolecule(m, ERGO_SPREFIX "/basis/6-31Gss",
                                   0, NULL, integralInfo, 0, 0, 0) != 0) {
    printf("bis.addBasisfuncsForMolecule failed.\n");
    return 1;
  }
  
  mat::SizesAndBlocks matrix_size_block_info =
    prepareMatrixSizesAndBlocks(bis.noOfBasisFuncs,
				20, 8, 8, 8);
  std::vector<int> permutationHML(bis.noOfBasisFuncs);
  getMatrixPermutation(bis, 20, 8, 8, 8, permutationHML);

  symmMatrix D;
  D.resetSizesAndBlocks(matrix_size_block_info,
			       matrix_size_block_info);
  symmMatrix J_1;
  J_1.resetSizesAndBlocks(matrix_size_block_info,
				 matrix_size_block_info);
  symmMatrix J_2;
  J_2.resetSizesAndBlocks(matrix_size_block_info,
				 matrix_size_block_info);
  symmMatrix J_diff;
  J_diff.resetSizesAndBlocks(matrix_size_block_info,
				    matrix_size_block_info);
  
  {
    /* Add values to density matrix diagonal and one step 
       next to diagonal. */
    const int nvalues1 = bis.noOfBasisFuncs;
    std::vector<int> idxrow(nvalues1);
    std::vector<int> idxcol(nvalues1);
    std::vector<ergo_real> values(nvalues1);
    for(int i=0; i<nvalues1; i++) {
      idxrow[i] = i;
      idxcol[i] = i;
      values[i] = 1.0;
    }
    D.add_values(idxrow, idxcol, values, 
		 permutationHML,
		 permutationHML);
    const int nvalues2 = bis.noOfBasisFuncs-1;
    for(int i=0; i<nvalues2; i++) {
      idxrow[i] = i;
      idxcol[i] = i+1;
      values[i] = 0.3;
    }
    idxrow.resize(nvalues2); 
    idxcol.resize(nvalues2); 
    values.resize(nvalues2); 
    D.add_values(idxrow, idxcol, values, 
		 permutationHML,
		 permutationHML);
  }

  JK::Params J_K_params;
  J_K_params. noOfThreads_J = defThreads;
  static const ergo_real EPS = mat::getMachineEpsilon<ergo_real>();
  J_K_params.threshold_J = template_blas_sqrt(EPS);

  J_K_params.fmm_box_size = 10;
  if(compute_J_by_boxes_sparse(bis, integralInfo, J_K_params, J_1, D,
			       permutationHML) != 0)
    {
      printf("Error in compute_J_by_boxes_sparse\n");
      return -1;
    }
  
  J_K_params.fmm_box_size = 4;
  if(compute_J_by_boxes_sparse(bis, integralInfo, J_K_params, J_2, D,
			       permutationHML) != 0)
    {
      printf("Error in compute_J_by_boxes_sparse\n");
      return -1;
    }
  
  J_diff = J_1;
  J_diff += (ergo_real)(-1.0) * J_2;

  ergo_real acc = template_blas_sqrt(EPS);
  ergo_real diffNorm = J_diff.eucl(acc);
  ergo_real requestedAcc;
  if(getenv("RUN_BENCHMARK")) 
    requestedAcc = J_K_params.threshold_J*50*nx*ny*nz;
  else
    requestedAcc = J_K_params.threshold_J*25*nx*ny*nz;
  if(diffNorm > requestedAcc)
    {
      printf("Error in J test: diff too large!\n");
      printf("diffNorm     = %8.4g\n", (double)diffNorm);
      printf("requestedAcc = %8.4g\n", (double)requestedAcc);
      return -1;
    }

  printf("J test OK, diffNorm = %7.4g\n", (double)diffNorm);
  return 0;
}
