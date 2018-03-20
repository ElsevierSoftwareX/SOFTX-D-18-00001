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

/** @file naive_response_test.cc

    @brief Tests naive implementation of linear response calculation.

    @author Elias Rudberg <em>responsible</em>
*/

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


const int MAX_AOS = 30;

typedef struct {
  ergo_real x[MAX_AOS][MAX_AOS][MAX_AOS][MAX_AOS];
} four_idx_AO_struct;


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
			       D, // use S as "density matrix" for J matrix test
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

void get_HML_dens_matrix(int noOfBasisFuncs,
			 int noOfOccupiedOrbitals,
			 const mat::SizesAndBlocks & sizeBlockInfo,
			 const std::vector<int> & permutationHML,
			 const std::vector<int> & inversePermutationHML,
			 symmMatrix & S,
			 symmMatrix & F,
			 symmMatrix & Dnew) {
  triangMatrix Z;
  Z.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  Z.inch(S, 0, mat::right);


  S.writeToFile();
  Z.writeToFile();

  // it is enough to set parameters which will be used or checked
  // if parameter is not specified nad checked somewhere during the execution, we should get an exception
  GetDensFromFock DensFromFock;
  DensFromFock.set_general_params(noOfBasisFuncs, sizeBlockInfo);
  DensFromFock.set_no_occupied_orbs(noOfOccupiedOrbitals);
  DensFromFock.set_use_diagonalization();
  DensFromFock.set_use_diag_on_error();
  DensFromFock.unset_use_purification();
  DensFromFock.unset_store_all_eigenvalues_to_file();

  ergo_real electronic_temperature = 0;
  DensFromFock.set_diagonalization_params(electronic_temperature,
					  S);

  DensFromFock.do_restricted_calculations(); // set factor = 2
  ergo_real invCholFactor_euclnorm = 0;
  DensFromFock.set_invCholFactor(Z, invCholFactor_euclnorm);
  ergo_real gap_expected_lower_bound = 0;
  DensFromFock.set_gap_expected_lower_bound(gap_expected_lower_bound);

  // before we set use_diagonalization flag
  // if purification will not be executed, these lines can be deleted
  ergo_real purification_subspace_err_limit = 1e-4;
  ergo_real purification_eigvalue_err_limit= 1e-4;
  ergo_real puri_eig_acc_factor_for_guess = 0;
  DensFromFock.set_purification_limits(purification_subspace_err_limit,
				       purification_eigvalue_err_limit, 
				       puri_eig_acc_factor_for_guess);

  DensFromFock.unset_use_acceleration();
  mat::normType normType = mat::frobNorm;
  DensFromFock.set_truncationNormPurification(normType);
  DensFromFock.set_stopCriterionNormPurification(normType);
  DensFromFock.unset_purification_create_m_files(); // this paramerter can be removed now
  DensFromFock.unset_output_homo_and_lumo_eigenvectors();
  DensFromFock.unset_use_new_stopping_criterion();
  DensFromFock.unset_use_diag_on_error_guess();
  /**********/


  DensFromFock.set_purification_maxmul(88);
  DensFromFock.unset_purification_ignore_failure();
  DensFromFock.unset_purification_use_rand_perturbation_for_alleigsint();
  DensFromFock.set_gap_expected_lower_bound(0);


  DensFromFock.clean_eigs_intervals();
  std::map<std::string, double> puri_stats;
  symmMatrix F_ort_prev_dummy;
  F_ort_prev_dummy.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  Dnew.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  Dnew.writeToFile();
  F.writeToFile();
  F_ort_prev_dummy.writeToFile();

  if(DensFromFock.get_dens_from_fock(F, 
				     Dnew,
				     F_ort_prev_dummy) != 0)
    throw "Error in get_dens_from_fock.";
	
  DensFromFock.get_puri_stats(puri_stats);

  S.readFromFile();
  Z.readFromFile();
  Dnew.readFromFile();
  F.readFromFile();
}

void do_HF_HML(int noOfBasisFuncs,
	       int noOfOccupiedOrbitals,
	       const mat::SizesAndBlocks & sizeBlockInfo,
	       const IntegralInfo & integralInfo,
	       const BasisInfoStruct & bis,
	       const std::vector<int> & permutationHML,
	       const std::vector<int> & inversePermutationHML,
	       const symmMatrix & S,
	       const symmMatrix & T,
	       const symmMatrix & V,
	       double nuclearRepulsionEnergy,
	       int noOfIterations,
	       symmMatrix & finalFockMatrix
	       ) {
  symmMatrix Scopy(S);
  symmMatrix Hcore;
  Hcore = T + V;
  symmMatrix D(S);
  int count = 0;
  while(1) {
    // Get new Fock matrix
    symmMatrix G;
    get_HML_G_matrix(sizeBlockInfo, integralInfo, bis, permutationHML, inversePermutationHML, D, G);
    double E = symmMatrix::trace_ab(D, Hcore) + 0.5 * symmMatrix::trace_ab(D, G) + nuclearRepulsionEnergy;
    printf("Energy %3d: %22.11f  (symmMatrix::trace_ab(D, G) = %12.6f, symmMatrix::trace_ab(D, Hcore) = %12.6f)\n",
	   count, E, (double)symmMatrix::trace_ab(D, G), (double)symmMatrix::trace_ab(D, Hcore));
    symmMatrix F(Hcore);
    F += G;
    symmMatrix Dnew;
    get_HML_dens_matrix(noOfBasisFuncs, noOfOccupiedOrbitals, sizeBlockInfo, permutationHML, inversePermutationHML, Scopy, F, Dnew);
    // Compute energy
    D = Dnew;
    count++;
    if(count == noOfIterations) {
      finalFockMatrix = F;
      break;
    }
  }
}

/* The get_matrices_A_and_B routine computed matrices A and B as given
   in equations (3.35) and (3.36) in the following paper: 
   "Polarization propagator methods in atomic and molecular calculations"
   Jens Oddershede, Poul Jørgensen, and Danny L. Yeager
   Computer Physics Reports
   Volume 2, Issue 2, November–December 1984, Pages 33-92
   http://dx.doi.org/10.1016/0167-7977(84)90003-0
*/
static void get_matrices_A_and_B(int nBasisFuncs,
				 int noOfOccupiedOrbitals,
				 ergo_real* A,
				 ergo_real* B,
				 const ergo_real* eigv,
				 const four_idx_AO_struct* g_MO) {
  int nBasf = nBasisFuncs;
  int nOcc = noOfOccupiedOrbitals;
  int n2 = nBasf*nBasf;
  for(int m = 0; m < nBasf; m++)
    for(int alpha = 0; alpha < nBasf; alpha++)
      for(int n = 0; n < nBasf; n++)
	for(int beta = 0; beta < nBasf; beta++) {
	  int idx1 = m * nBasf + alpha;
	  int idx2 = n * nBasf + beta;
	  if(m < nOcc || n < nOcc || alpha >= nOcc || beta >= nOcc) {
	    A[idx1*n2+idx2] = 0;
	    B[idx1*n2+idx2] = 0;
	    continue;
	  }
	  int delta_m_n = 0;
	  if(m == n)
	    delta_m_n = 1;
	  int delta_alpha_beta = 0;
	  if(alpha == beta)
	    delta_alpha_beta = 1;
	  A[idx1*n2+idx2] = (eigv[m] - eigv[alpha]) * delta_m_n * delta_alpha_beta
	    + 2*g_MO->x[m][alpha][beta][n] - g_MO->x[m][n][beta][alpha];
	  B[idx1*n2+idx2] = g_MO->x[alpha][n][beta][m] - 2 * g_MO->x[alpha][m][beta][n];
	}
}

static void get_all_generalized_eigenvalues(int n, const ergo_real* A_in, const ergo_real* B_in, ergo_real* eigvalList) {
  int lwork = 8*n*n;
  std::vector<ergo_real> work(lwork);
  std::vector<ergo_real> A(n*n);
  std::vector<ergo_real> B(n*n);
  memcpy(&A[0], A_in, n*n*sizeof(ergo_real));
  memcpy(&B[0], B_in, n*n*sizeof(ergo_real));
  std::vector<ergo_real> alpha_r(n);
  std::vector<ergo_real> alpha_i(n);
  std::vector<ergo_real> beta(n);
  ergo_real* vl = NULL;
  ergo_real* vr = NULL;
  int ldvl = 1;
  int ldvr = 1;
  int info = -1;
  mat::ggev("N", "N", &n, &A[0], &n, &B[0], &n, &alpha_r[0], 
	    &alpha_i[0], &beta[0], vl, &ldvl, 
	    vr, &ldvr, &work[0], &lwork, 
	    &info);
  if(info != 0)
    throw std::runtime_error("Error in mat::ggev: (info != 0)");
  // Check that all eigenvalues are real
  for(int i = 0; i < n; i++) {
    if(alpha_i[i] != 0)
      throw std::runtime_error("Error, eigenvalue from mat::ggev with nonzero imaginary part found.");
  }
  for(int i = 0; i < n; i++)
    eigvalList[i] = alpha_r[i] / beta[i];
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
  
  Util::TimeMeter tmEverything;
  std::string moleculeStr;
  if(argc > 1)
    moleculeStr = argv[1];
  else
    moleculeStr = "default";
  int blockSizeHML = 16;

  printf("===== Parameters (set by command-line args in the same order as below) =======\n");
  printf("moleculeStr = '%s'\n", moleculeStr.c_str());
  printf("===== End of parameters ======================================================\n");

  IntegralInfo integralInfo(true);
  BasisInfoStruct bis;
  Molecule m;

  if(moleculeStr == "default") {
    // Use default h2o molecule.
    m.addAtom(8, 0, 0, 0);
    m.addAtom(1, -1.809, 0, 0);
    m.addAtom(1, 0.453549, 1.751221, 0);
    int nAtoms = m.getNoOfAtoms();
    printf("Using default molecule, nAtoms = %d\n", nAtoms);
  }
  else {
    // We expect a filename for a molecule file.
    if(m.setFromMoleculeFile(moleculeStr.c_str(), 
			     0, /* we are guessing the net charge here */
			     NULL) != 0) {
      std::cerr << "Error in setFromMoleculeFile for filename '" << moleculeStr << "'." << std::endl;
      throw std::runtime_error("Error in m.setFromMoleculeFile().");
    }
    // Verify that nAtoms is -1 in this case.
    int nAtoms = m.getNoOfAtoms();
    printf("Molecule file '%s' read OK, nAtoms = %d\n", moleculeStr.c_str(), nAtoms);
  }

  Util::TimeMeter tmPartB;

  m.setNetCharge(0);
  int noOfElectrons = m.getNumberOfElectrons();
  if(noOfElectrons <= 0 || noOfElectrons % 2 != 0) {
    printf("Error: noOfElectrons = %d. Need even number of electrons.\n", noOfElectrons);
    return -1;
  }
  int noOfOccupiedOrbitals = noOfElectrons / 2;
  printf("noOfElectrons = %5d  ==>  noOfOccupiedOrbitals = %5d\n", noOfElectrons, noOfOccupiedOrbitals);

  int noOfIterationsHF = 25;

  const char* basisSetName = "STO-3G";

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

  ergo_real nuclearRepulsionEnergy = m.getNuclearRepulsionEnergyQuadratic();

  std::vector<Atom> atomList(m.getNoOfAtoms());
  for(int i = 0; i < m.getNoOfAtoms(); i++)
    atomList[i] = m.getAtom(i);

  if(bis.addBasisfuncsForMolecule(m, basisSetName,
				  0, NULL, integralInfo, 0, 0, 0) != 0)
    throw std::runtime_error("bis.addBasisfuncsForMolecule failed.");
  int nBasisFuncs = bis.noOfBasisFuncs;
  printf("nBasisFuncs = %d\n", nBasisFuncs);

  // Get HML overlap matrix
  std::vector<int> permutationHML, inversePermutationHML;
  mat::SizesAndBlocks sizeBlockInfo;
  preparePermutationsHML(bis, sizeBlockInfo, permutationHML, inversePermutationHML, blockSizeHML);
  symmMatrix S_notrunc;
  S_notrunc.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo); 
  if(compute_overlap_matrix_sparse(bis, S_notrunc, 
				   permutationHML) != 0)
    throw std::runtime_error("Error in compute_overlap_matrix_sparse.");

  symmMatrix Ssymm(S_notrunc);
  printf("Ssymm.eucl() = %9.5f\n", (double)Ssymm.eucl(1e-6));

  symmMatrix T_notrunc;
  symmMatrix V_notrunc;

  // Get HML T matrix
  T_notrunc.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo); 
  ergo_real threshold_for_T = 1e-12;
  ergo_real box_size_for_T = 5.8;
  if(compute_T_sparse_linear(bis,
			     integralInfo,
			     threshold_for_T,
			     box_size_for_T,
			     T_notrunc,
			     permutationHML) != 0)
    throw std::runtime_error("Error in compute_T_sparse.");

  // Get HML V matrix
  V_notrunc.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo); 
  ergo_real threshold_for_V = 1e-12;
  ergo_real boxSize_for_V = 8.8;
  ergo_real nuclearRepulsionEnergy_dummy = 0;
  if(compute_V_sparse(bis, 
		      integralInfo,
		      m,
		      threshold_for_V,
		      boxSize_for_V,
		      V_notrunc,
		      permutationHML,
		      nuclearRepulsionEnergy_dummy) != 0)
    throw std::runtime_error("Error in compute_V_sparse.");

  symmMatrix FockMatrix;
  printf("Calling do_HF_HML()\n");
  do_HF_HML(nBasisFuncs,
	    noOfOccupiedOrbitals,
	    sizeBlockInfo,
	    integralInfo,
	    bis,
	    permutationHML,
	    inversePermutationHML,
	    S_notrunc,
	    T_notrunc,
	    V_notrunc,
	    nuclearRepulsionEnergy,
	    noOfIterationsHF,
	    FockMatrix);
  printf("After do_HF_HML()\n");

  int n = nBasisFuncs;
  std::vector<ergo_real> S_full(n*n);
  S_notrunc.fullMatrix(S_full, inversePermutationHML, inversePermutationHML);
  std::vector<ergo_real> F_full(n*n);
  FockMatrix.fullMatrix(F_full, inversePermutationHML, inversePermutationHML);

  // Get eigenvectors and eigenvalues of F (also known as molecular orbitals and orbital energies)
  std::vector<ergo_real> MOs(n*n);
  std::vector<ergo_real> eigv(n);
  get_F_orbs(n, &F_full[0], &S_full[0], &MOs[0], &eigv[0]);

  // Get two-electron integrals
  printf("Computing two-electron integrals...\n");
  four_idx_AO_struct* g_AO = new four_idx_AO_struct;
  for(int p = 0; p < n; p++)
    for(int q = 0; q < n; q++)
      for(int r = 0; r < n; r++)
	for(int s = 0; s < n; s++)
	  g_AO->x[p][q][r][s] = do_2e_integral(p, q, r, s, bis, integralInfo);
  printf("Two-electron integrals done.\n");
  
  // Get two-electron integrals in MO basis
  printf("Computing two-electron integrals in MO basis...\n");
  four_idx_AO_struct* g_MO = new four_idx_AO_struct;
  for(int p = 0; p < n; p++)
    for(int q = 0; q < n; q++)
      for(int r = 0; r < n; r++)
	for(int s = 0; s < n; s++) {
	  ergo_real sum = 0;
	  for(int a = 0; a < n; a++)
	    for(int b = 0; b < n; b++)
	      for(int c = 0; c < n; c++)
		for(int d = 0; d < n; d++)
		  sum += MOs[p*n+a] * MOs[q*n+b] * MOs[r*n+c] * MOs[s*n+d] * g_AO->x[a][b][c][d];
	  g_MO->x[p][q][r][s] = sum;
	}
  printf("Two-electron integrals in MO basis done.\n");

  // Set up matrices A and B needed for response equations
  int n2=n*n;
  std::vector<ergo_real> A(n2*n2);
  std::vector<ergo_real> B(n2*n2);
  get_matrices_A_and_B(n, noOfOccupiedOrbitals, &A[0], &B[0], &eigv[0], g_MO);

  // Set up matrix of double size for linear response eigenvalue equation
  int n2b=2*n2;
  std::vector<ergo_real> M(n2b*n2b);
  // Insert A into upper left corner of M
  for(int p = 0; p < n2; p++)
    for(int q = 0; q < n2; q++)
      M[p*n2b+q] = A[p*n2+q];
  // Insert A into lower right corner of M
  for(int p = 0; p < n2; p++)
    for(int q = 0; q < n2; q++)
      M[(p+n2)*n2b+(q+n2)] = A[p*n2+q];
  // Insert B into upper right corner of M
  for(int p = 0; p < n2; p++)
    for(int q = 0; q < n2; q++)
      M[p*n2b+(q+n2)] = B[p*n2+q];
  // Insert B into lower left corner of M
  for(int p = 0; p < n2; p++)
    for(int q = 0; q < n2; q++)
      M[(p+n2)*n2b+q] = B[p*n2+q];
  // Set up diagonal matrix II with values 1 and -1 on diagonal 
  std::vector<ergo_real> II(n2b*n2b);
  for(int i = 0; i < n2b; i++)
    for(int j = 0; j < n2b; j++) {
      ergo_real value = 0;
      if(i == j) {
	if(i < n2)
	  value = 1;
	else
	  value = -1;
      }
      II[i*n2b+j] = value;
    }

  // Get generalized eigenvalues of M and II
  std::vector<ergo_real> eigvalList(n2b);
  get_all_generalized_eigenvalues(n2b, &M[0], &II[0], &eigvalList[0]);
  std::sort(eigvalList.begin(), eigvalList.begin()+n2b);

  printf("Nonzero eigenvalues of M and II:\n");
  for(int i = 0; i < n2b; i++) {
    if(template_blas_fabs(eigvalList[i]) > 1e-6)
      printf("Eigenvalue %5d: %12.8f\n", i, (double)eigvalList[i]);
  }

  if(moleculeStr == "default") {
    // Verify that computed eigenvalues match reference values that were computed using the standard (less naive) implementation in Ergo.
    ergo_real maxabsdiff = 0;
    ergo_real refValues[] = {0.483554699, 0.556648356, 0.612541560, 0.702754741};
    for(int i = 0; i < 4; i++) {
      ergo_real absdiff = template_blas_fabs(eigvalList[88+i] - refValues[i]);
      if(absdiff > maxabsdiff)
	maxabsdiff = absdiff;
    }
    ergo_real tolerance = 1e-7;
#ifdef PRECISION_SINGLE
    tolerance = 1e-5;
#endif
    printf("Comparison of eigenvalues to ref values, maxabsdiff = %9.5g, tolerance = %9.5g\n", (double)maxabsdiff, (double)tolerance);
    if(maxabsdiff > tolerance)
      throw std::runtime_error("Error: maxabsdiff too large.");
  }

  // FIXME: Solve linear response equations for a particular value of the frequency omega, to compute polarizability values.

  puts("Naive response test succeeded."); 
  report_timing(tmEverything, "Everything");
  return 0;
}

