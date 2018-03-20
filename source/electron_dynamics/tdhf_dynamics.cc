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

/** @file tdhf_dynamics.cc

    @brief Implementation of TDHF electron dynamics, explicitly
    propagating the electron density in time in the precence of a
    time-dependent electric field.

    @author: Elias Rudberg <em>responsible</em>

    This implementation was done to reproduce the results in the article
    by Li et al: Phys. Chem. Chem. Phys., 2005, 7, 233–239.

    DOI link for that article: http://dx.doi.org/10.1039/B415849K
*/

#include <cstring>
#include "tdhf_dynamics.h"
#include "scf_utils.h"
#include "output.h"
#include "pi.h"
#include "integral_matrix_wrappers.h"
#include "integrals_2el_explicit.h"


#if 0
static void print_matrix(int n, const normalMatrix & M, const char* name, std::vector<int> const & inversePermutationHML) {
  std::vector<ergo_real> M_full(n*n);
  M.fullMatrix(M_full, inversePermutationHML, inversePermutationHML);
  printf("matrix '%s':\n", name);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++)
      printf("%13.7f", (double)M_full[i*n+j]);
    printf("\n");
  }
}
#endif


static void get_Fock_matrix_in_ort_basis(const BasisInfoStruct & basisInfo,
					 const IntegralInfo & integralInfo,
					 const JK::ExchWeights & CAM_params,
					 const JK::Params & J_K_params,
					 normalMatrix & F_ort, 
					 const normalMatrix & D_ort, 
					 const triangMatrix & invCholFactor, 
					 const symmMatrix & H_core_Matrix,
					 const mat::SizesAndBlocks & size_block_info,
					 const symmMatrix & refFockMatrix,
					 std::vector<int> const & permutationHML,
					 std::vector<int> const & inversePermutationHML,
					 bool realPart)
{
  Util::TimeMeter timeMeter;

  do_output(LOG_CAT_INFO, LOG_AREA_ED, "get_Fock_matrix_in_ort_basis() starting!");

  normalMatrix Z;
  Z.resetSizesAndBlocks(size_block_info, size_block_info);
  Z = invCholFactor;

  normalMatrix ZT;
  ZT.resetSizesAndBlocks(size_block_info, size_block_info);
  ZT = transpose(Z);

  normalMatrix ZD;
  ZD.resetSizesAndBlocks(size_block_info, size_block_info);
  ZD = Z * D_ort;
  normalMatrix D_nonort;
  D_nonort.resetSizesAndBlocks(size_block_info, size_block_info);
  D_nonort = ZD * ZT;

  // Compute J.
  normalMatrix J;
  J.resetSizesAndBlocks(size_block_info, size_block_info);
  J = 0;
  if(realPart) {
    symmMatrix Jsymm;
    Jsymm.resetSizesAndBlocks(size_block_info, size_block_info);
    symmMatrix D_nonort_symm;
    D_nonort_symm.resetSizesAndBlocks(size_block_info, size_block_info);
    D_nonort_symm = D_nonort;
    D_nonort_symm *= (ergo_real)2;
    if(compute_J_by_boxes_sparse(basisInfo,
				 integralInfo,
				 J_K_params,
				 Jsymm,
				 D_nonort_symm,
				 permutationHML) != 0)
      throw "Error in compute_J_by_boxes_sparse.";  
    J = Jsymm;
  }

  // Compute K.
  normalMatrix K;
  K.resetSizesAndBlocks(size_block_info, size_block_info);
  if(D_nonort.frob() > 1e-13) {
    D_nonort *= (ergo_real)2;
    if(compute_K_by_boxes_sparse_nosymm(basisInfo,
					integralInfo, 
					CAM_params,
					J_K_params,
					K,
					D_nonort,
					permutationHML,
					inversePermutationHML) != 0)
      throw "Error in compute_K_by_boxes_sparse_nonsymm";
  }

  normalMatrix G;
  G.resetSizesAndBlocks(size_block_info, size_block_info);
  G = J + K;

  // OK, now we have G in non-orth basis. We also need H_core.
  normalMatrix F_nonort;
  F_nonort.resetSizesAndBlocks(size_block_info, size_block_info);
  if(realPart)
    F_nonort = H_core_Matrix;
  else
    F_nonort = 0;
  F_nonort += (ergo_real)1 * G;

  // Compare to reference Fock matrix.
  normalMatrix refF;
  refF.resetSizesAndBlocks(size_block_info, size_block_info);
  refF = refFockMatrix;
  ergo_real frobDiff = normalMatrix::frob_diff(F_nonort, refF);
  do_output(LOG_CAT_INFO, LOG_AREA_ED, "get_Fock_matrix_in_ort_basis() frobDiff = %9.4g.", (double)frobDiff);

  // Now transform final result back to orthogonal basis.
  normalMatrix ZTF;
  ZTF.resetSizesAndBlocks(size_block_info, size_block_info);
  ZTF = ZT * F_nonort;
  F_ort = ZTF * Z;

  timeMeter.print(LOG_AREA_ED, "get_Fock_matrix_in_ort_basis");

  do_output(LOG_CAT_INFO, LOG_AREA_ED, "get_Fock_matrix_in_ort_basis() finishing OK!");
}


struct ComplexMatrix {
  normalMatrix re;
  normalMatrix im;
  mat::SizesAndBlocks size_block_info;
  void initialize(mat::SizesAndBlocks szblinfo) {
    size_block_info = szblinfo;
    re.resetSizesAndBlocks(size_block_info, size_block_info);
    im.resetSizesAndBlocks(size_block_info, size_block_info);
  }
  void copy(const ComplexMatrix & X) {
    re = X.re;
    im = X.im;
  }
  void product(const ComplexMatrix & X, const ComplexMatrix & Y) {
    normalMatrix tmp1;
    tmp1.resetSizesAndBlocks(size_block_info, size_block_info);
    normalMatrix tmp2;
    tmp2.resetSizesAndBlocks(size_block_info, size_block_info);
    // First do re part: re = Xre*Yre - Xim*Yim
    tmp1 = X.re * Y.re;
    tmp2 = X.im * Y.im;
    re = tmp1;
    re += (ergo_real)-1 * tmp2;
    // Now do im part: im = Xre*Yim + Xim*Yre
    tmp1 = X.re * Y.im;
    tmp2 = X.im * Y.re;
    im = tmp1 + tmp2;
  }
  void sum(const ComplexMatrix & X, const ComplexMatrix & Y) {
    re = X.re + Y.re;
    im = X.im + Y.im;
  }
  void rescale(ergo_real a) {
    re *= a;
    im *= a;
  }
  void rescale_im(ergo_real a) {
    normalMatrix re_new;
    re_new.resetSizesAndBlocks(size_block_info, size_block_info);
    re_new = im;
    re_new *= -a;
    im = re;
    im *= a;
    re = re_new;
  }
  void add(const ComplexMatrix & X, ergo_real a) {
    re += a * X.re;
    im += a * X.im;
  }
  void do_conjugate() {
    normalMatrix tmp;
    tmp.resetSizesAndBlocks(size_block_info, size_block_info);
    tmp = transpose(re);
    re = tmp;
    tmp = transpose(im);
    im = tmp;
    im *= (ergo_real)-1;
  }
};


#if 0
static void printComplexMatrix(int n, const ComplexMatrix & M, const char* name, 
			       std::vector<int> const & inversePermutationHML) {
  printf("=========== printComplexMatrix ===============\n");
  char s[888];
  sprintf(s, "%s re", name);
  print_matrix(n, M.re, s, inversePermutationHML);
  sprintf(s, "%s im", name);
  print_matrix(n, M.im, s, inversePermutationHML);
}

static void printComplexMatrix2(int n, const ComplexMatrix & M, const char* name, 
				std::vector<int> const & inversePermutationHML) {
  printf("=========== printComplexMatrix2 ===============\n");
  std::vector<ergo_real> M_re_full(n*n);
  M.re.fullMatrix(M_re_full, inversePermutationHML, inversePermutationHML);
  std::vector<ergo_real> M_im_full(n*n);
  M.im.fullMatrix(M_im_full, inversePermutationHML, inversePermutationHML);
  printf("matrix '%s':\n", name);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++)
      printf("%12.8f + %12.8fi            ", (double)M_re_full[i*n+j], (double)M_im_full[i*n+j]);
    printf("\n");
  }
}
#endif


#if 0
/* do_exp_transform_BCH(): compute exp(X) * Y * exp(-X) using BCH
   expansion. */
static void do_exp_transform_BCH(ComplexMatrix & result, 
				 const ComplexMatrix & X, 
				 const ComplexMatrix & Y,
				 mat::SizesAndBlocks size_block_info) {
  result.copy(Y);
  ComplexMatrix prevMainTerm;
  prevMainTerm.initialize(size_block_info);
  prevMainTerm.copy(Y);
  ergo_real preFactor = 1;
  const int n = 10;
  for(int i = 1; i <= n; i++) {
    ComplexMatrix tmp1;
    tmp1.initialize(size_block_info);
    tmp1.product(X, prevMainTerm);
    ComplexMatrix tmp2;
    tmp2.initialize(size_block_info);
    tmp2.product(prevMainTerm, X);
    tmp2.rescale(-1);
    ComplexMatrix commutator;
    commutator.initialize(size_block_info);
    commutator.sum(tmp1, tmp2);
    prevMainTerm.copy(commutator);
    result.add(commutator, preFactor);
    preFactor = preFactor / 2;
  }
}
#endif


static void compute_exp_of_matrix(ComplexMatrix & U, 
				  const ComplexMatrix & X, 
				  mat::SizesAndBlocks size_block_info,
				  int n,
				  std::vector<int> const & inversePermutationHML) {
  Util::TimeMeter timeMeter;
  ComplexMatrix X_to_pow_k;
  X_to_pow_k.initialize(size_block_info);
  X_to_pow_k.re = 1;
  X_to_pow_k.im = 0;
  U.copy(X_to_pow_k);
  ergo_real preFactor = 1;
  const int ORDER_OF_TAYLOR_EXPANSION = 11;
  for(int k = 1; k < ORDER_OF_TAYLOR_EXPANSION; k++) {
    preFactor /= k;
    ComplexMatrix tmp;
    tmp.initialize(size_block_info);
    tmp.product(X_to_pow_k, X);
    X_to_pow_k = tmp;
    U.add(X_to_pow_k, preFactor);
  }
  timeMeter.print(LOG_AREA_ED, "compute_exp_of_matrix");
}


static void compute_U_X_Udagger(ComplexMatrix & result, 
				const ComplexMatrix & U, 
				const ComplexMatrix & X, 
				mat::SizesAndBlocks size_block_info) {
  ComplexMatrix UX;
  UX.initialize(size_block_info);
  UX.product(U, X);
  ComplexMatrix Udagger;
  Udagger.initialize(size_block_info);
  Udagger = U;
  Udagger.do_conjugate();
  result.product(UX, Udagger);  
}


static void get_curr_electric_field(ergo_real* electricField, ergo_real t) {
  const ergo_real omega = 0.10;
  const ergo_real E_max = 0.07;
  const ergo_real phasePhi = 0;
  ergo_real E = 0;
  if(t < 2*pi/omega)
    E = E_max * omega*t/(2*pi);
  else if(t < 4*pi/omega)
    E = E_max;
  else if(t < 6*pi/omega)
    E = E_max * (3 - omega*t/(2*pi));
  else
    E = 0;
  electricField[0] = 0;
  electricField[1] = 0;
  electricField[2] = E * template_blas_sin(omega*t + phasePhi);
}


static ergo_real vectorLength(ergo_real x, ergo_real y, ergo_real z) {
  return template_blas_sqrt(x*x + y*y + z*z);
}


void do_tdhf_dynamics(const BasisInfoStruct & basisInfo,
		      const IntegralInfo & integralInfo,
		      const Molecule & molecule,
		      const Molecule & extraCharges,
		      const SCF::MatOptions& matOpts,
		      const JK::ExchWeights & CAM_params,
		      const JK::Params & J_K_params,
		      const symmMatrix & FockMatrix,
		      const symmMatrix & densityMatrix,
		      const symmMatrix & S_symm,
		      const triangMatrix & invCholFactor,
		      const ED::Params & params)
{
  do_output(LOG_CAT_INFO, LOG_AREA_ED, "do_tdhf_dynamics() starting!");

  // First get density matrix to orthogonal basis.

  normalMatrix Dtest;
  Dtest.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  Dtest = densityMatrix;

  normalMatrix SD;
  SD.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  SD = S_symm * densityMatrix;
  normalMatrix SDS;
  SDS.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  SDS = SD * S_symm;
  symmMatrix D_ort;
  D_ort.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  D_ort = SDS;
  D_ort = transpose(invCholFactor) * D_ort * invCholFactor;
    
  D_ort *= (ergo_real)0.5;
  
  do_output(LOG_CAT_INFO, LOG_AREA_ED, "SCF_restricted::do_electron_dynamics(), D_ort done.");

  // Test D_ort by converting back to non-orthogonal basis and comparing to D.
  symmMatrix D_nonort_test;
  D_nonort_test.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  D_nonort_test = D_ort;
  D_nonort_test = invCholFactor * D_nonort_test * transpose(invCholFactor);
  ergo_real D_frob_diff = symmMatrix::frob_diff(D_nonort_test, densityMatrix);

  do_output(LOG_CAT_INFO, LOG_AREA_ED, 
	    "SCF_restricted::do_electron_dynamics(), D_ort test D_frob_diff = %8.3g.", (double)D_frob_diff);

  // Get matrices needed for computation of instantaneous dipole.
  symmMatrix dipoleMatrixSymm_x, dipoleMatrixSymm_y, dipoleMatrixSymm_z;
  dipoleMatrixSymm_x.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  dipoleMatrixSymm_y.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  dipoleMatrixSymm_z.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  if(compute_operator_matrix_sparse_symm(basisInfo, 
					 1, 0, 0, 
					 dipoleMatrixSymm_x, 
					 matOpts.permutationHML) != 0)
    throw "Error in compute_operator_matrix_sparse_symm.";
  if(compute_operator_matrix_sparse_symm(basisInfo, 
					 0, 1, 0, 
					 dipoleMatrixSymm_y, 
					 matOpts.permutationHML) != 0)
    throw "Error in compute_operator_matrix_sparse_symm.";
  if(compute_operator_matrix_sparse_symm(basisInfo, 
					 0, 0, 1, 
					 dipoleMatrixSymm_z, 
					 matOpts.permutationHML) != 0)
    throw "Error in compute_operator_matrix_sparse_symm.";
  normalMatrix dipoleMatrix_x, dipoleMatrix_y, dipoleMatrix_z;
  dipoleMatrix_x.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  dipoleMatrix_y.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  dipoleMatrix_z.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
  dipoleMatrix_x = dipoleMatrixSymm_x;
  dipoleMatrix_y = dipoleMatrixSymm_y;
  dipoleMatrix_z = dipoleMatrixSymm_z;

  // We use two density matrices for the propagation, curr and
  // prev. Each of them has a real and an imaginary part.
  ComplexMatrix D_ort_curr;
  D_ort_curr.initialize(matOpts.size_block_info);
  ComplexMatrix D_ort_prev;
  D_ort_prev.initialize(matOpts.size_block_info);
  // We start out with zero imaginary part, and the real part for prev
  // is the same as curr.
  D_ort_curr.im.clear();
  D_ort_prev.im.clear();
  D_ort_curr.re = D_ort;
  D_ort_prev.re = D_ort;

  // Now do the dynamics!
  int stepCount = 0;
  ergo_real timestep = params.timestep;
  while(1) {
    Util::TimeMeter timeMeterCurrStep;

    // Before constructing Fock matrix, we need the 1-electron matrix
    // computed for the current electric field, since the field can be
    // time-dependent.
    ergo_real electricField[3];
    ergo_real t = stepCount * timestep;

    if(t > params.max_time)
      break;

    get_curr_electric_field(electricField, t);
    symmMatrix H_core_Matrix_curr;
    H_core_Matrix_curr.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
    ergo_real threshold_integrals_1el = 1e-13;
    int no_of_threads_for_V = 1;
    ergo_real box_size_for_V_and_T = 7.3;
    Util::TimeMeter timeMeterHCore;
    ergo_real nuclearRepulsionEnergy_dummy = 0;
    if(compute_h_core_matrix_sparse(integralInfo, 
				    molecule, 
				    extraCharges,
				    electricField[0],
				    electricField[1],
				    electricField[2],
				    basisInfo, 
				    H_core_Matrix_curr,
				    threshold_integrals_1el,
				    no_of_threads_for_V,
				    box_size_for_V_and_T,
				    nuclearRepulsionEnergy_dummy,
				    matOpts.size_block_info,
				    matOpts.permutationHML) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_ED, "error in compute_h_core_matrix_sparse");
	throw "error in compute_h_core_matrix_sparse";
      }
    timeMeterHCore.print(LOG_AREA_ED, "TDHF dynamics step H_core matrix computation");

    // Get Fock matrix. To do this, first transform density matrix re
    // and im parts to non-orthogonal basis.
    ComplexMatrix F_ort;
    F_ort.initialize(matOpts.size_block_info);

    get_Fock_matrix_in_ort_basis(basisInfo, integralInfo, CAM_params, J_K_params,
				 F_ort.re, D_ort_curr.re, invCholFactor, H_core_Matrix_curr, matOpts.size_block_info, 
				 FockMatrix, matOpts.permutationHML, matOpts.inversePermutationHML, true);
    get_Fock_matrix_in_ort_basis(basisInfo, integralInfo, CAM_params, J_K_params,
				 F_ort.im, D_ort_curr.im, invCholFactor, H_core_Matrix_curr, matOpts.size_block_info, 
				 FockMatrix, matOpts.permutationHML, matOpts.inversePermutationHML, false);
    // Get matrix X = 2*i*timestep*F_ort
    ComplexMatrix X;
    X.initialize(matOpts.size_block_info);
    X.copy(F_ort);
    X.rescale_im(2 * timestep);
#if 0
    // Now compute exp(X)*D_ort*exp(-X)
    ComplexMatrix D_ort_new;
    D_ort_new.initialize(matOpts.size_block_info);
    do_exp_transform_BCH(D_ort_new, X, D_ort_prev, matOpts.size_block_info);
    do_output(LOG_CAT_INFO, LOG_AREA_ED, 
	      "SCF_restricted::do_electron_dynamics(), do_exp_transform_BCH done!");
#else
    const int n = basisInfo.noOfBasisFuncs;
    // Now compute U = exp(X) and U*D_ort*Udagger
    ComplexMatrix U;
    U.initialize(matOpts.size_block_info);
    compute_exp_of_matrix(U, X, matOpts.size_block_info, n, matOpts.inversePermutationHML);
    ComplexMatrix D_ort_new;
    D_ort_new.initialize(matOpts.size_block_info);
    compute_U_X_Udagger(D_ort_new, U, D_ort_prev, matOpts.size_block_info);
#endif

    
    D_ort_prev.copy(D_ort_curr);
    D_ort_curr.copy(D_ort_new);
    
    stepCount++;

    // We want to output some information every "whole" time unit. The
    // step size is typically smaller than 1 a.u. so now we figure out
    // if this timestep is closest to the nearest "whole" time unit.
    ergo_real nearestWholeTime = round((long double)t);
    ergo_real currDiff = template_blas_fabs(t - nearestWholeTime);
    ergo_real prevDiff = template_blas_fabs(t - timestep - nearestWholeTime);
    ergo_real nextDiff = template_blas_fabs(t + timestep - nearestWholeTime);
    if(currDiff <= prevDiff && currDiff <= nextDiff) {
      // Compute instantaneous dipole.
      ergo_real instantaneousDipole_re, instantaneousDipole_im;
      // re
      {
	normalMatrix D_nonort_curr;
	D_nonort_curr.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
	normalMatrix Z;
	Z.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
	Z = invCholFactor;
	normalMatrix ZD;
	ZD.resetSizesAndBlocks(matOpts.size_block_info, matOpts.size_block_info);
	ZD = Z * D_ort_curr.re;
	D_nonort_curr = ZD * transpose(Z);
	ergo_real instantaneousDipole_x = -2 * normalMatrix::trace_ab(D_nonort_curr, dipoleMatrix_x);
	ergo_real instantaneousDipole_y = -2 * normalMatrix::trace_ab(D_nonort_curr, dipoleMatrix_y);
	ergo_real instantaneousDipole_z = -2 * normalMatrix::trace_ab(D_nonort_curr, dipoleMatrix_z);
	// Add contributions for each nucleus in molecule.
	for(int k = 0; k < molecule.getNoOfAtoms(); k++) {
	  const Atom & atom = molecule.getAtom(k);
	  //printf("Atom %d coords: %9.5f %9.5f %9.5f\n", k, atom.coords[0], atom.coords[1], atom.coords[2]);
	  instantaneousDipole_x += atom.charge * atom.coords[0];
	  instantaneousDipole_y += atom.charge * atom.coords[1];
	  instantaneousDipole_z += atom.charge * atom.coords[2];
	}
	instantaneousDipole_re = instantaneousDipole_z;
      }
      // im
      {
	ergo_real instantaneousDipole_x = normalMatrix::trace_ab(D_ort_curr.im, dipoleMatrix_x);
	ergo_real instantaneousDipole_y = normalMatrix::trace_ab(D_ort_curr.im, dipoleMatrix_y);
	ergo_real instantaneousDipole_z = normalMatrix::trace_ab(D_ort_curr.im, dipoleMatrix_z);
	instantaneousDipole_im = 
	  vectorLength(instantaneousDipole_x, instantaneousDipole_y, instantaneousDipole_z);
      }
      do_output(LOG_CAT_INFO, LOG_AREA_ED, "Data for dipole plot: %15.0f %15.10f %15.10f (time was really %15.10f)", 
		(double)t, (double)instantaneousDipole_re, (double)instantaneousDipole_im, (double)t);
    }

    timeMeterCurrStep.print(LOG_AREA_ED, "TDHF dynamics complete step");
  }

  do_output(LOG_CAT_INFO, LOG_AREA_ED, "do_tdhf_dynamics() finished!");
}
