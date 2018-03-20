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

/** @file many_h_atoms_test.cc Tests the simple case of many
    well-separated H atoms with just one basis function for each
    atom. */

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include "basisinfo.h"

#include "matrix_utilities.h"
#include "integrals_general.h"
#include "integrals_2el_explicit.h"
#include "integral_matrix_wrappers.h"
#include "utilities.h"
#include "densfromf_full.h"
#include "GetDensFromFock.h"


static void preparePermutationsHML(const BasisInfoStruct& basisInfo,
				   mat::SizesAndBlocks& sizeBlockInfo, 
				   std::vector<int>& permutation,
				   std::vector<int>& inversePermutation,
				   int blockSizeHML) {
  const int sparseMatrixBlockFactor = 4;
  sizeBlockInfo = prepareMatrixSizesAndBlocks(basisInfo.noOfBasisFuncs,
					      blockSizeHML,
					      sparseMatrixBlockFactor,
					      sparseMatrixBlockFactor,
					      sparseMatrixBlockFactor);
  getMatrixPermutation(basisInfo,
                       blockSizeHML,
                       sparseMatrixBlockFactor,
                       sparseMatrixBlockFactor,
                       sparseMatrixBlockFactor,
                       permutation,
                       inversePermutation);
}

static void report_timing(const Util::TimeMeter & tm, const char* s) {
  double secondsTaken = Util::TimeMeter::get_wall_seconds() - tm.get_start_time_wall_seconds();
  printf("'%s' took %12.5f wall seconds.\n", s, secondsTaken);
}

void get_HML_J(const mat::SizesAndBlocks & sizeBlockInfo,
	       const IntegralInfo & integralInfo,
	       const BasisInfoStruct & bis,
	       const std::vector<int> & permutationHML,
	       const symmMatrix & D,
	       symmMatrix & J) {
  J.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  JK::Params J_K_params;
  if(compute_J_by_boxes_sparse(bis,
			       integralInfo,
			       J_K_params,
			       J,
			       D,
			       permutationHML) != 0)
    throw std::runtime_error("Error in compute_J_by_boxes_sparse.");
}

void get_HML_K(const mat::SizesAndBlocks & sizeBlockInfo,
	       const IntegralInfo & integralInfo,
	       const BasisInfoStruct & bis,
	       const std::vector<int> & permutationHML,
	       const std::vector<int> & inversePermutationHML,
	       symmMatrix & D,
	       symmMatrix & K) {
  K.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  JK::Params J_K_params;
  JK::ExchWeights CAM_params;
  if(compute_K_by_boxes_sparse(bis,
			       integralInfo,
			       CAM_params,
			       J_K_params,
			       K,
			       D,
			       permutationHML,
			       inversePermutationHML) != 0)
    throw std::runtime_error("Error in compute_K_by_boxes_sparse.");
}

void get_HML_G_matrix(const mat::SizesAndBlocks & sizeBlockInfo,
		      const IntegralInfo & integralInfo,
		      const BasisInfoStruct & bis,
		      const std::vector<int> & permutationHML,
		      const std::vector<int> & inversePermutationHML,
		      symmMatrix & D,
		      symmMatrix & G) {
  symmMatrix J, K;
  get_HML_J(sizeBlockInfo, integralInfo, bis, permutationHML, D, J);
  get_HML_K(sizeBlockInfo, integralInfo, bis, permutationHML, inversePermutationHML, D, K);
  G = J;
  G += K;
}

static void get_diag_dens_mat(symmMatrix & D,
			      int n,
			      int offset_0_or_1,
			      const mat::SizesAndBlocks & sizeBlockInfo,
			      const std::vector<int> & permutationHML) {
  D.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  std::vector<int> rowind(n);
  std::vector<int> colind(n);
  std::vector<ergo_real> values(n);
  for(int i = 0; i < n; i++) {
    rowind[i] = i;
    colind[i] = i;
    values[i] = 0;
    if(i % 2 == offset_0_or_1)
      values[i] = 1;
  }
  D.assign_from_sparse(rowind, colind, values, permutationHML, permutationHML);
}

static void report_max_abs_element_of_matrix(const symmMatrix & A,
					     const char* name,
					     const std::vector<int> & inversePermutationHML) {
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values;
  A.get_all_values(rowind,
		   colind,
		   values,
		   inversePermutationHML,
		   inversePermutationHML);
  ergo_real maxabsvalue = 0;
  ergo_real savedValue = 0;
  long nvalues = values.size();
  for(long i = 0; i < nvalues; i++) {
    ergo_real value = values[i];
    ergo_real absvalue = template_blas_fabs(value);
    if(absvalue > maxabsvalue) {
      maxabsvalue = absvalue;
      savedValue = value;
    }
  }
  printf("report_max_abs_element_of_matrix() for '%s': maxabsvalue = %g, savedValue = %g\n", name, (double)maxabsvalue, (double)savedValue);
}


static void do_add_atom(Molecule & m, ergo_real spacing, int ix, int iy, int iz) {
  ergo_real x = ix*spacing;
  ergo_real y = iy*spacing;
  ergo_real z = iz*spacing;
  m.addAtom(1, x, y, z);
}

static ergo_real get_energy(int nx, int ny, int nz, int & nAtoms, ergo_real boxSizeForJ, bool specialCase) {
  
  IntegralInfo integralInfo(true);
  BasisInfoStruct basisInfo;

  // Set up molecule
  Util::TimeMeter tmSetUpMolecule;
  const ergo_real spacing = 32.7; // this should be a large enough spacing so that the interaction between atoms becomes negligible
  Molecule m;
  if(specialCase) {
    // Specific case with 29 atoms that turned out to expose a problem with the multipole method implementation.
    do_add_atom(m, spacing, 0, 7, 1);
    do_add_atom(m, spacing, 5, 5, 5);
    do_add_atom(m, spacing, 5, 5, 7);
    do_add_atom(m, spacing, 5, 5, 8);
    do_add_atom(m, spacing, 5, 6, 5);
    do_add_atom(m, spacing, 5, 6, 7);
    do_add_atom(m, spacing, 5, 6, 8);
    do_add_atom(m, spacing, 5, 7, 5);
    do_add_atom(m, spacing, 5, 7, 7);
    do_add_atom(m, spacing, 5, 7, 8);
    do_add_atom(m, spacing, 6, 5, 5);
    do_add_atom(m, spacing, 6, 5, 7);
    do_add_atom(m, spacing, 6, 5, 8);
    do_add_atom(m, spacing, 6, 6, 5);
    do_add_atom(m, spacing, 6, 6, 7);
    do_add_atom(m, spacing, 6, 7, 5);
    do_add_atom(m, spacing, 6, 7, 7);
    do_add_atom(m, spacing, 6, 7, 8);
    do_add_atom(m, spacing, 7, 0, 8);
    do_add_atom(m, spacing, 7, 5, 5);
    do_add_atom(m, spacing, 7, 5, 7);
    do_add_atom(m, spacing, 7, 5, 8);
    do_add_atom(m, spacing, 7, 6, 5);
    do_add_atom(m, spacing, 7, 6, 7);
    do_add_atom(m, spacing, 7, 6, 8);
    do_add_atom(m, spacing, 7, 7, 0);
    do_add_atom(m, spacing, 7, 7, 5);
    do_add_atom(m, spacing, 7, 7, 7);
    do_add_atom(m, spacing, 7, 7, 8);
  }
  else {
    for(int ix = 0; ix < nx; ix++)
      for(int iy = 0; iy < ny; iy++)
	for(int iz = 0; iz < nz; iz++) {
	  ergo_real x = ix*spacing;
	  ergo_real y = iy*spacing;
	  ergo_real z = iz*spacing;
	  m.addAtom(1, x, y, z);
	}
  }
  nAtoms = m.getNoOfAtoms();
  printf("nAtoms = %d\n", nAtoms);
  m.setNetCharge(0);
  int noOfElectrons = m.getNumberOfElectrons();
  if(noOfElectrons != nAtoms)
    throw std::runtime_error("Error: (noOfElectrons != nAtoms).");
  report_timing(tmSetUpMolecule, "Setting up molecule");

  const char* basisSetName = "STO-1G";

  static const char *dirv[] = {
    ".", "basis", "../basis",
    ERGO_DATA_PREFIX "/basis",
    ERGO_DATA_PREFIX,
    ERGO_SPREFIX "/basis",
    ERGO_SPREFIX
  };
  basisset_info basissetDef;
  if(read_basisset_file(basissetDef, basisSetName, 6, dirv, 0) != 0)
    throw std::runtime_error("Error in read_basisset_file().");

  if(basisInfo.addBasisfuncsForMolecule(m, basisSetName,
					0, NULL, integralInfo, 0, 0, 0) != 0)
    throw std::runtime_error("bis.addBasisfuncsForMolecule failed.");
  int nBasisFuncs = basisInfo.noOfBasisFuncs;
  printf("nBasisFuncs = %d\n", nBasisFuncs);
  if(nBasisFuncs != nAtoms)
    throw std::runtime_error("Error: (nBasisFuncs != nAtoms).");

  const int blockSizeHML = 16;

  ergo_real common_threshold = 1e-7;

  ergo_real threshold_for_V = common_threshold;
  ergo_real boxSize_for_V = 78.5;
  ergo_real box_size_for_T = 65.1;
  JK::Params J_K_params;
  J_K_params.threshold_J = common_threshold;
  J_K_params.threshold_K = common_threshold;
  J_K_params.fmm_box_size = boxSizeForJ;
  J_K_params.exchange_box_size = 77.4;

  // Get HML overlap matrix
  std::vector<int> permutationHML, inversePermutationHML;
  mat::SizesAndBlocks sizeBlockInfo;
  preparePermutationsHML(basisInfo, sizeBlockInfo, permutationHML, inversePermutationHML, blockSizeHML);
  symmMatrix S_notrunc;
  S_notrunc.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo); 
  if(compute_overlap_matrix_sparse(basisInfo, S_notrunc, permutationHML) != 0)
    throw std::runtime_error("Error in compute_overlap_matrix_sparse.");
  report_max_abs_element_of_matrix(S_notrunc, "S_notrunc", inversePermutationHML);

  symmMatrix Ssymm(S_notrunc);
  printf("Ssymm.eucl() = %9.5f\n", (double)Ssymm.eucl(1e-6));

  symmMatrix T_notrunc;
  symmMatrix V_notrunc;

  // Get HML T matrix
  Util::TimeMeter tmGetMatrixT;
  T_notrunc.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo); 
  ergo_real threshold_for_T = 1e-12;
  if(compute_T_sparse_linear(basisInfo,
			     integralInfo,
			     threshold_for_T,
			     box_size_for_T,
			     T_notrunc,
			     permutationHML) != 0)
    throw std::runtime_error("Error in compute_T_sparse_linear.");
  report_timing(tmGetMatrixT, "GetMatrixT");
  report_max_abs_element_of_matrix(T_notrunc, "T_notrunc", inversePermutationHML);

  // Get HML V matrix
  Util::TimeMeter tmGetMatrixV;
  V_notrunc.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  ergo_real nuclearRepulsionEnergy = 0;
  if(compute_V_sparse(basisInfo,
		      integralInfo,
		      m,
		      threshold_for_V,
		      boxSize_for_V,
		      V_notrunc,
		      permutationHML,
		      nuclearRepulsionEnergy) != 0)
    throw std::runtime_error("Error in compute_V_sparse.");
  report_timing(tmGetMatrixV, "GetMatrixV");
  report_max_abs_element_of_matrix(V_notrunc, "V_notrunc", inversePermutationHML);
  printf("nuclearRepulsionEnergy (computed by compute_V_sparse)  = %22.11f\n", (double)nuclearRepulsionEnergy);

  symmMatrix D_alpha;
  symmMatrix D_beta;
  get_diag_dens_mat(D_alpha, nBasisFuncs, 0, sizeBlockInfo, permutationHML);
  get_diag_dens_mat(D_beta , nBasisFuncs, 1, sizeBlockInfo, permutationHML);

  symmMatrix D_tot(D_alpha);
  D_tot += D_beta;

  Util::TimeMeter tmGetMatrixJ;
  symmMatrix J;
  J.resetSizesAndBlocks(sizeBlockInfo,sizeBlockInfo);
  if(compute_J_by_boxes_sparse(basisInfo,
			       integralInfo,
			       J_K_params,
			       J,
			       D_tot,
			       permutationHML) != 0)
    throw std::runtime_error("error in compute_J_by_boxes_sparse");
  report_timing(tmGetMatrixJ, "GetMatrixJ");
  report_max_abs_element_of_matrix(J, "J", inversePermutationHML);

  Util::TimeMeter tmGetMatrixK;
  symmMatrix K_alpha;
  K_alpha.resetSizesAndBlocks(sizeBlockInfo,sizeBlockInfo);
  symmMatrix K_beta;
  K_beta.resetSizesAndBlocks(sizeBlockInfo,sizeBlockInfo);
  JK::ExchWeights CAM_params_not_used;
  if(compute_K_by_boxes_sparse(basisInfo, integralInfo, CAM_params_not_used,
			       J_K_params, K_alpha, D_alpha,
			       permutationHML, inversePermutationHML) != 0)
    throw std::runtime_error("error in compute_K_by_boxes_sparse");
  if(compute_K_by_boxes_sparse(basisInfo, integralInfo, CAM_params_not_used,
			       J_K_params, K_beta, D_beta,
			       permutationHML, inversePermutationHML) != 0)
    throw std::runtime_error("error in compute_K_by_boxes_sparse");
  report_timing(tmGetMatrixK, "GetMatrixK");
  report_max_abs_element_of_matrix(K_alpha, "K_alpha", inversePermutationHML);
  report_max_abs_element_of_matrix(K_beta , "K_beta ", inversePermutationHML);

  K_alpha *= 2;
  K_beta *= 2;

  symmMatrix twoelMatrix_sparse_alpha;
  twoelMatrix_sparse_alpha = J + K_alpha;
  symmMatrix twoelMatrix_sparse_beta;
  twoelMatrix_sparse_beta = J + K_beta;

  ergo_real energy_2el = 0;
  energy_2el +=  0.5 * symmMatrix::trace_ab(D_alpha, twoelMatrix_sparse_alpha);
  energy_2el +=  0.5 * symmMatrix::trace_ab(D_beta , twoelMatrix_sparse_beta );

  symmMatrix H_core;
  H_core = T_notrunc + V_notrunc;
  ergo_real energy_1el = symmMatrix::trace_ab(D_tot, H_core);

  printf("energy_1el              = %22.11f\n", (double)energy_1el);
  printf("energy_2el              = %22.11f\n", (double)energy_2el);
  printf("nuclearRepulsionEnergy  = %22.11f\n", (double)nuclearRepulsionEnergy);

  ergo_real energy = energy_1el + energy_2el + nuclearRepulsionEnergy;

  {
    // Other way of computing energy
    ergo_real halfTraceOfDJ = 0.5 * symmMatrix::trace_ab(D_tot, J);
    ergo_real exchange_energy_alpha = 0.5 * symmMatrix::trace_ab(D_alpha, K_alpha);
    ergo_real exchange_energy_beta = 0.5 * symmMatrix::trace_ab(D_beta, K_beta);
    ergo_real kinetic_energy = symmMatrix::trace_ab(D_tot, T_notrunc);
    ergo_real halfTraceDV = 0.5 * symmMatrix::trace_ab(D_tot, V_notrunc);
    ergo_real energy_alternative1 = halfTraceOfDJ + exchange_energy_alpha + exchange_energy_beta + kinetic_energy + halfTraceDV + halfTraceDV + nuclearRepulsionEnergy;

    printf("halfTraceOfDJ = %22.11f\n", (double)halfTraceOfDJ);
    printf("halfTraceDV   = %22.11f\n", (double)halfTraceDV);

    symmMatrix J_plus_V(J);
    J_plus_V += V_notrunc;
    ergo_real halfTraceOfDtimesJplusV = 0.5 * symmMatrix::trace_ab(D_tot, J_plus_V);
    ergo_real energy_alternative2 = halfTraceOfDtimesJplusV + exchange_energy_alpha + exchange_energy_beta + kinetic_energy + halfTraceDV + nuclearRepulsionEnergy;

    printf("energy              = %22.11f\n", (double)energy);
    printf("energy_alternative1 = %22.11f\n", (double)energy_alternative1);
    printf("energy_alternative2 = %22.11f\n", (double)energy_alternative2);

    ergo_real testSum1 = halfTraceOfDJ + exchange_energy_alpha + exchange_energy_beta;
    printf("testSum1    = %22.11f\n", (double)testSum1);
    printf("energy_2el  = %22.11f\n", (double)energy_2el);
  }
  
  return energy;
}

static int do_energy_comparison(int nx, int ny, int nz, ergo_real boxSizeForJ, bool specialCase) {
  int nAtoms1 = 0;
  ergo_real E_one_atom = get_energy(1, 1, 1, nAtoms1, boxSizeForJ, false);
  assert(nAtoms1 == 1);
  int nAtoms = 0;
  ergo_real E_many_atoms = get_energy(nx, ny, nz, nAtoms, boxSizeForJ, specialCase);
  ergo_real expectedEnergy = E_one_atom * nAtoms;
  ergo_real absDiff = template_blas_fabs(E_many_atoms - expectedEnergy);
  ergo_real relDiff = absDiff / template_blas_fabs(expectedEnergy);
  printf("E_one_atom   = %22.11f\n", (double)E_one_atom);
  printf("E_many_atoms = %22.11f\n", (double)E_many_atoms);
  printf("absDiff = %g\n", (double)absDiff);
  printf("relDiff = %g\n", (double)relDiff);
  ergo_real tolerance = 1e-5;
#ifdef PRECISION_SINGLE
  tolerance = 2e-2;
#endif
  if(relDiff > tolerance) {
    printf("Error: relDiff too large.\n");
    return -1;
  }
  return 0;
}

int main(int argc, char *argv[])
{
  Util::TimeMeter tmEverything;
  int nx = 2;
  int ny = 2;
  int nz = 2;
  if(argc == 4) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nz = atoi(argv[3]);
  }
  printf("many_h_atoms_test, nx ny nz = %d %d %d\n", nx, ny, nz);

#if 0
  // Do this if you want the ergoscf.out file, to see timings etc.
  enable_output();
  enable_memory_usage_output();
#endif

#ifdef _OPENMP
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  if ( !(env && (defThreads=atoi(env)) > 0) )
    defThreads = 1;
  mat::Params::setNProcs(defThreads);
  mat::Params::setMatrixParallelLevel(2);
  std::cout<<"OpenMP is used, number of threads set to "
	   <<mat::Params::getNProcs()<<". Matrix parallel level: "
	   <<mat::Params::getMatrixParallelLevel()<<"."<<std::endl;
#endif

  ergo_real boxSize_for_J = 77.9;

  bool try_one_special_case = true;
  if(try_one_special_case) {
    if(do_energy_comparison(nx, ny, nz, boxSize_for_J, true) != 0) {
      printf("FAILED.\n");
      return -1;
    }
  }
  bool try_many_box_sizes = false;
  if(try_many_box_sizes) {
    for(int i = 0; i < 100; i++) {
      ergo_real boxSizeForJ_curr = 7.9 + 0.9*i;
      printf("i = %d, boxSizeForJ_curr = %f\n", i, (double)boxSizeForJ_curr);
      if(do_energy_comparison(nx, ny, nz, boxSizeForJ_curr, true) != 0) {
	printf("FAILED.\n");
	return -1;
      }
    }
  } // end if try_many_box_sizes

  if(do_energy_comparison(nx, ny, nz, boxSize_for_J, false) != 0) {
    printf("FAILED.\n");
    return -1;
  }

  puts("many_h_atoms_test finished OK.");
  report_timing(tmEverything, "Everything");
  return 0;
}

