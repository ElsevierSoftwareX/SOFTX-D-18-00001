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

/** @file jmat_extrapolate_test.cc Tests the error extrapolation scheme
    for the Coulomb matrix. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include <vector>

#include "scf_utils.h"
#include "matrix_utilities.h"
#include "mat_acc_extrapolate.h"
#include "integral_matrix_wrappers.h"


class Jworkertest
{
public:
  Jworkertest(const symmMatrix & D_,
	      const IntegralInfo & integralInfo_,
	      const BasisInfoStruct & basisInfo_,
	      const JK::Params & J_K_params_,
	      std::vector<int> const & permutationHML_);
  void ComputeMatrix(ergo_real param,
		     symmMatrix & result) const;
private:
  const symmMatrix & D;
  const IntegralInfo & integralInfo;
  const BasisInfoStruct & basisInfo;
  JK::Params J_K_params;
  std::vector<int> const & permutationHML;
};

Jworkertest::Jworkertest(const symmMatrix & D_,
			 const IntegralInfo & integralInfo_,
			 const BasisInfoStruct & basisInfo_,
			 const JK::Params & J_K_params_,
			 std::vector<int> const & permutationHML_) : 
  D(D_), 
  integralInfo(integralInfo_), 
  basisInfo(basisInfo_),
  permutationHML(permutationHML_)
{
  J_K_params = J_K_params_;
}

void Jworkertest::ComputeMatrix(ergo_real param,
			    symmMatrix & result) const
{
  JK::Params J_K_params_tmp = J_K_params;
  J_K_params_tmp.threshold_J = param;
  if(compute_J_by_boxes_sparse(basisInfo, 
			       integralInfo, 
			       J_K_params_tmp, 
			       result, 
			       D,
			       permutationHML) != 0)
    throw "Jworkertest::ComputeMatrix: error in compute_J_by_boxes_sparse";
}


int main(int argc, char *argv[])
{
#ifdef _OPENMP
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  if ( !(env && (defThreads=atoi(env)) > 0) ) {
    defThreads = 1;
  }  
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
      nx = 4;
      ny = 4;
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
    idxrow.resize(nvalues2); 
    idxcol.resize(nvalues2); 
    values.resize(nvalues2); 
    for(int i=0; i<nvalues2; i++) {
      idxrow[i] = i;
      idxcol[i] = i+1;
      values[i] = 0.3;
    }
    D.add_values(idxrow, idxcol, values, 
		 permutationHML,
		 permutationHML);
  }

  JK::Params J_K_params;




  Jworkertest worker(D, integralInfo, bis, J_K_params, permutationHML);
  

  MatAccInvestigator<ergo_real, Jworkertest> 
    investigator(matrix_size_block_info);
  int nSteps = 4;
  investigator.Scan(worker, 1e-8, sqrt(sqrt(sqrt(10.0))), nSteps);
  
  ergo_real threshList[nSteps];
  ergo_real errorList_frob[nSteps];
  ergo_real errorList_eucl[nSteps];
  ergo_real errorList_maxe[nSteps];
  ergo_real timeList[nSteps];
  investigator.GetScanResult(threshList, 
			     errorList_frob,
			     errorList_eucl,
			     errorList_maxe,
			     timeList);

  printf("Scan result:\n");
  for(int i = 0; i < nSteps; i++)
    printf("%9.5g %9.5g %9.5g %9.5g\n", 
	   (double)threshList[i], 
	   (double)errorList_frob[i],
	   (double)errorList_eucl[i],
	   (double)errorList_maxe[i]);



#if 0
  static const ergo_real EPS = mat::getMachineEpsilon<ergo_real>();
  J_K_params.threshold = sqrt(EPS);

  J_K_params.fmm_box_size = 10;
  if(compute_J_by_boxes_sparse(bis, integralInfo, J_K_params, J_1, D) != 0)
    {
      printf("Error in compute_J_by_boxes_sparse\n");
      return -1;
    }
  
  J_K_params.fmm_box_size = 4;
  if(compute_J_by_boxes_sparse(bis, integralInfo, J_K_params, J_2, D) != 0)
    {
      printf("Error in compute_J_by_boxes_sparse\n");
      return -1;
    }
  
  J_diff = J_1;
  J_diff += (ergo_real)(-1.0) * J_2;

  generalVector vector;
  vector.build_data_structure(perm);
  ergo_real acc = sqrt(EPS);
  ergo_real diffNorm = J_diff.eucl(vector, acc);
  ergo_real requestedAcc;
  if(getenv("RUN_BENCHMARK")) 
    requestedAcc = J_K_params.threshold*4000;
  else
    requestedAcc = J_K_params.threshold*200;
  if(diffNorm > requestedAcc)
    {
      printf("Error in J test: diff too large!\n");
      printf("diffNorm     = %8.4g\n", diffNorm);
      printf("requestedAcc = %8.4g\n", requestedAcc);
      return -1;
    }
#endif

  unlink("ergoscf.out");

  printf("J extrapolate test OK\n");
  return 0;
}
