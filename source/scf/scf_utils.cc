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

/** @file scf_utils.cc

    @brief Various utilities used by self-consistent field (SCF)
    code. For example, interface routines converting between HML
    matrix format used in main SCF code and elementwise format used by
    integral code.

    @author: Elias Rudberg <em>responsible</em>
*/

#include <pthread.h>
#include <sstream>

#include "scf_utils.h"
#include "output.h"
#include "integrals_1el.h"
#include "memorymanag.h"
#include "operator_matrix.h"
#include "integrals_1el_kinetic.h"
#include "integrals_1el_potential.h"
#include "utilities.h"
#include "integrals_2el_explicit.h"
#include "basis_func_pair_list.h"
#include "integrals_2el_boxed.h"
#include "integrals_2el_explicit.h"
#include "integrals_2el_layer.h"
#include "density_description_file.h"
#include "mat_acc_extrapolate.h"
#include "basis_func_pair_list_1el.h"
#include "integral_matrix_wrappers.h"
#include "units.h"
#include "densitymanager.h"


/* Flag for skipping DFT stuff, used in QD precision tests. */
//#define SKIP_DFT_FLAG

#ifndef SKIP_DFT_FLAG
#include "dft_common.h"
#include "xc_matrix.h"
#include "xc_matrix_sparse.h"
#endif




template<class Tinvestigator, class Tworker>
static void do_scan_and_report(Tinvestigator investigator, 
			       Tworker worker, 
			       const char* scanName,
			       int nSteps,
			       ergo_real startThresh,
			       ergo_real stepFactor)
{
  investigator.Scan(worker, startThresh, stepFactor, nSteps);
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
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Accuracy scan result for '%s':", scanName);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "      thresh            frob            eucl            maxe            time");
  for(int i = 0; i < nSteps; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "%15.6g %15.6g %15.6g %15.6g %15.6g",
	      (double)threshList[i],
	      (double)errorList_frob[i],
	      (double)errorList_eucl[i],
	      (double)errorList_maxe[i],
	      (double)timeList[i]);
}





class Jworker
{
public:
  Jworker(const symmMatrix & D_,
	  const IntegralInfo & integralInfo_,
	  const BasisInfoStruct & basisInfo_,
	  const triangMatrix & invCholFactor_,
	  bool doInvCholFactorTransformation_,
	  const JK::Params & J_K_params_,
	  std::vector<int> const & permutationHML_) :
    D(D_), 
    integralInfo(integralInfo_), 
    basisInfo(basisInfo_),
    invCholFactor(invCholFactor_),
    doInvCholFactorTransformation(doInvCholFactorTransformation_),
    permutationHML(permutationHML_)
  {
    J_K_params = J_K_params_;
  }
  void ComputeMatrix(ergo_real param,
		     symmMatrix & result) const;
private:
  const symmMatrix & D;
  const IntegralInfo & integralInfo;
  const BasisInfoStruct & basisInfo;
  const triangMatrix & invCholFactor;
  bool doInvCholFactorTransformation;
  JK::Params J_K_params;
  std::vector<int> const & permutationHML;
};

void Jworker::ComputeMatrix(ergo_real param,
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
    throw "Jworker::ComputeMatrix: error in compute_J_by_boxes_sparse";
  if(doInvCholFactorTransformation)
    result = transpose(invCholFactor) * result * invCholFactor;
}

void
do_acc_scan_J(const symmMatrix & D,
	      const IntegralInfo & integralInfo,
	      const BasisInfoStruct & basisInfo,
	      triangMatrix & invCholFactor,
	      bool doInvCholFactorTransformation,
	      const JK::Params & J_K_params,
	      mat::SizesAndBlocks const & matrix_size_block_info,
	      std::vector<int> const & permutationHML,
	      int nSteps,
	      ergo_real startThresh,
	      ergo_real stepFactor)
{
  invCholFactor.readFromFile();
  Jworker worker(D, integralInfo, basisInfo, invCholFactor, doInvCholFactorTransformation, J_K_params, permutationHML);
  MatAccInvestigator<ergo_real, Jworker> investigator(matrix_size_block_info);
  do_scan_and_report(investigator, worker, "J", nSteps, startThresh, stepFactor);
  invCholFactor.writeToFile();
}



class Kworker
{
public:
  Kworker(symmMatrix & D_,
	  const IntegralInfo & integralInfo_,
	  const BasisInfoStruct & basisInfo_,
	  const triangMatrix & invCholFactor_,
	  bool doInvCholFactorTransformation_,
	  const JK::ExchWeights & CAM_params_,
	  const JK::Params & J_K_params_,
	  std::vector<int> const & permutationHML_,
	  std::vector<int> const & inversePermutationHML_) : 
    D(D_), 
    integralInfo(integralInfo_), 
    basisInfo(basisInfo_),
    invCholFactor(invCholFactor_),
    doInvCholFactorTransformation(doInvCholFactorTransformation_),
    CAM_params(CAM_params_),
    permutationHML(permutationHML_),
    inversePermutationHML(inversePermutationHML_)
  {
    J_K_params = J_K_params_;
  }
  void ComputeMatrix(ergo_real param,
		     symmMatrix & result) const;
private:
  symmMatrix & D;
  const IntegralInfo & integralInfo;
  const BasisInfoStruct & basisInfo;
  const triangMatrix & invCholFactor;
  bool doInvCholFactorTransformation;
  const JK::ExchWeights & CAM_params;
  JK::Params J_K_params;
  std::vector<int> const & permutationHML;
  std::vector<int> const & inversePermutationHML;
};

void Kworker::ComputeMatrix(ergo_real param,
			    symmMatrix & result) const
{
  JK::Params J_K_params_tmp = J_K_params;
  J_K_params_tmp.threshold_K = param;
  if(compute_K_by_boxes_sparse(basisInfo, 
			       integralInfo, 
			       CAM_params,
			       J_K_params_tmp, 
			       result, 
			       D,
			       permutationHML,
			       inversePermutationHML) != 0)
    throw "Kworker::ComputeMatrix: error in compute_K_by_boxes_sparse";
  if(doInvCholFactorTransformation)
    result = transpose(invCholFactor) * result * invCholFactor;
}

void
do_acc_scan_K(symmMatrix & D,
	      const IntegralInfo & integralInfo,
	      const BasisInfoStruct & basisInfo,
	      triangMatrix & invCholFactor,
	      bool doInvCholFactorTransformation,
	      const JK::ExchWeights & CAM_params,
	      const JK::Params & J_K_params,
	      mat::SizesAndBlocks const & matrix_size_block_info,
	      std::vector<int> const & permutationHML,
	      std::vector<int> const & inversePermutationHML,
	      int nSteps,
	      ergo_real startThresh,
	      ergo_real stepFactor)
{
  invCholFactor.readFromFile();
  Kworker worker(D, integralInfo, basisInfo, invCholFactor, doInvCholFactorTransformation, CAM_params, J_K_params, permutationHML, inversePermutationHML);
  MatAccInvestigator<ergo_real, Kworker> investigator(matrix_size_block_info);
  do_scan_and_report(investigator, worker, "K", nSteps, startThresh, stepFactor);
  invCholFactor.writeToFile();
}


class Vxc_worker
{
public:
  Vxc_worker(symmMatrix & D_,
	     const IntegralInfo & integralInfo_,
	     const BasisInfoStruct & basisInfo_,
	     const Molecule & molecule_,
	     const Dft::GridParams & gridParams_,
	     int noOfElectrons_,
	     const triangMatrix & invCholFactor_,
	     bool doInvCholFactorTransformation_,
	     mat::SizesAndBlocks const & matrix_size_block_info_,
	     std::vector<int> const & permutationHML_,
	     std::vector<int> const & inversePermutationHML_) : 
    D(D_), 
    integralInfo(integralInfo_), 
    basisInfo(basisInfo_),
    molecule(molecule_),
    gridParams(gridParams_),
    noOfElectrons(noOfElectrons_),
    invCholFactor(invCholFactor_),
    doInvCholFactorTransformation(doInvCholFactorTransformation_),
    matrix_size_block_info(matrix_size_block_info_),
    permutationHML(permutationHML_),
    inversePermutationHML(inversePermutationHML_)
  {
  }
  void ComputeMatrix(ergo_real param,
		     symmMatrix & result) const;
private:
  symmMatrix & D;
  const IntegralInfo & integralInfo;
  const BasisInfoStruct & basisInfo;
  const Molecule & molecule;
  const Dft::GridParams & gridParams;
  int noOfElectrons;
  const triangMatrix & invCholFactor;
  bool doInvCholFactorTransformation;
  mat::SizesAndBlocks const & matrix_size_block_info;
  std::vector<int> const & permutationHML;
  std::vector<int> const & inversePermutationHML;
};

void Vxc_worker::ComputeMatrix(ergo_real param,
			       symmMatrix & result) const
{
  Dft::GridParams gridParams_tmp = gridParams;
  gridParams_tmp.hicuParams.maxError = param;
  result.resetSizesAndBlocks(matrix_size_block_info,
			     matrix_size_block_info);
  /* Remove grid files to make sure new grid is generated. */
  grid_free_files();
  ergo_real dftEnergy = 0;
  Dft::getXC_mt(basisInfo, integralInfo,
		molecule, gridParams_tmp, noOfElectrons,
		D, result, &dftEnergy,
		permutationHML);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Vxc_worker::ComputeMatrix dftEnergy = %25.15f", (double)dftEnergy);
  if(doInvCholFactorTransformation)
    result = transpose(invCholFactor) * result * invCholFactor;
}

void
do_acc_scan_Vxc(symmMatrix & D,
		const IntegralInfo & integralInfo,
		const BasisInfoStruct & basisInfo,
		const Molecule & molecule,
		const Dft::GridParams & gridParams,
		int noOfElectrons,
		triangMatrix & invCholFactor,
		bool doInvCholFactorTransformation,
		mat::SizesAndBlocks const & matrix_size_block_info,
		std::vector<int> const & permutationHML,
		std::vector<int> const & inversePermutationHML,
		int nSteps,
		ergo_real startThresh,
		ergo_real stepFactor)
{
  invCholFactor.readFromFile();
  Vxc_worker worker(D, 
		    integralInfo, 
		    basisInfo, 
		    molecule,
		    gridParams,
		    noOfElectrons,
		    invCholFactor, 
		    doInvCholFactorTransformation, 
		    matrix_size_block_info, 
		    permutationHML, 
		    inversePermutationHML);
  MatAccInvestigator<ergo_real, Vxc_worker> investigator(matrix_size_block_info);
  do_scan_and_report(investigator, worker, "Vxc", nSteps, startThresh, stepFactor);
  invCholFactor.writeToFile();
}





template <typename Tmatrix>
void output_sparsity_template(int n, 
			      const Tmatrix & M, 
			      const char* matrixName, 
			      const char* matrixTypeName) {
  /* Use double to avoid integer overflow. */
  double nnz_as_double = M.nnz();
  ergo_real percent = 100 * nnz_as_double / ((double)n*n);
  ergo_real memUsage = (double)M.memory_usage() / 1000000000;
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "Matrix '%s' (%s): n = %6i, nnz = %12.0f <-> %6.2f %%, mem usage %7.4f G", 
	    matrixName, matrixTypeName, n, nnz_as_double, (double)percent, (double)memUsage);  
}


void output_sparsity(int n, const normalMatrix & M, const char* matrixName)
{
  output_sparsity_template(n, M, matrixName, "normal");
}

void output_sparsity_symm(int n, const symmMatrix & M, const char* matrixName)
{
  output_sparsity_template(n, M, matrixName, " symm ");
}

void output_sparsity_triang(int n, const triangMatrix & M, const char* matrixName)
{
  output_sparsity_template(n, M, matrixName, "triang");
}


static ergo_real
get_trace(int n, const ergo_real* P, const ergo_real* H_core)
{
  ergo_real sum = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      sum += P[i*n+j] * H_core[i*n+j];
  return sum;
}


static int
add_square_matrices(int n, const ergo_real* A, const ergo_real* B, ergo_real* C)
{
  int n2 = n*n;
  for(int i = 0; i < n2; i++)
    C[i] = A[i] + B[i];
  return 0;
}


int
compute_h_core_matrix_simple_dense(const IntegralInfo& integralInfo,
				   const Molecule& molecule,
				   const BasisInfoStruct& basisInfo,
				   symmMatrix & H_core_Matrix_sparse,
				   ergo_real threshold_integrals_1el,
				   int noOfThreadsForV,
				   mat::SizesAndBlocks const & matrix_size_block_info,
				   std::vector<int> const & permutationHML,
				   ergo_real & result_nuclearRepulsionEnergy)
{
  int n = basisInfo.noOfBasisFuncs;
  result_nuclearRepulsionEnergy = molecule.getNuclearRepulsionEnergyQuadratic();
  std::vector<ergo_real> A(n*n);
  if(compute_h_core_matrix_full(integralInfo,
				basisInfo,
				molecule.getNoOfAtoms(),
				molecule.getAtomListPtr(),
				&A[0],
				threshold_integrals_1el) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_h_core_matrix_full");
      return -1;
    }
  int nvalues = n * (n + 1) / 2;
  std::vector<int> rowind(nvalues);
  std::vector<int> colind(nvalues);
  std::vector<ergo_real> values(nvalues);
  // populate vectors
  int count = 0;
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      rowind[count] = i;
      colind[count] = j;
      values[count] = A[i*n+j];
      count++;
    }
  }
  H_core_Matrix_sparse.assign_from_sparse(rowind,
					  colind,
					  values,
					  permutationHML,
					  permutationHML);
  return 0;
}



int
compute_h_core_matrix_sparse(const IntegralInfo& integralInfo,
			     const Molecule& molecule,
			     const Molecule& extraCharges,
			     ergo_real electric_field_x,
			     ergo_real electric_field_y,
			     ergo_real electric_field_z,
			     const BasisInfoStruct& basisInfo,
			     symmMatrix & H_core_Matrix_sparse,
			     ergo_real threshold_integrals_1el,
			     int noOfThreadsForV,
			     ergo_real boxSizeForVT,
			     ergo_real & result_nuclearRepulsionEnergy,
			     mat::SizesAndBlocks const & matrix_size_block_info,
			     std::vector<int> const & permutationHML,
			     int const create_dipole_mtx,
			     std::vector<int> const * const inversePermutationHML,
			     std::string const * const calculation_identifier,
			     std::string const * const method_and_basis_set)
{
  int n = basisInfo.noOfBasisFuncs;
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "entering compute_h_core_matrix_sparse, nAtoms = %i, n = %i, threshold = %g", molecule.getNoOfAtoms(), n, (double)threshold_integrals_1el);

  Util::TimeMeter timeMeterTot;
  
  // Get T

  Util::TimeMeter timeMeterT;

  symmMatrix T(H_core_Matrix_sparse);
  T.clear();

  if(compute_T_sparse_linear(basisInfo,
			     integralInfo,
			     threshold_integrals_1el,
			     boxSizeForVT,
			     T,
			     permutationHML) != 0) {
    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_T_sparse");
    return -1;
  }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "compute_T_sparse_linear returned OK.");
  timeMeterT.print(LOG_AREA_SCF, "compute_T_sparse_linear");

  // Get V
  Util::TimeMeter timeMeterV;
  symmMatrix V;
  V.resetSizesAndBlocks(matrix_size_block_info,
			matrix_size_block_info);
  if(compute_V_sparse(basisInfo,
		      integralInfo, 
		      molecule,
		      threshold_integrals_1el,
		      boxSizeForVT,
		      V,
		      permutationHML,
		      result_nuclearRepulsionEnergy) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_V_sparse");
      return -1;
    }    
  timeMeterV.print(LOG_AREA_SCF, "Computation of V");

  H_core_Matrix_sparse = T + V;
  T.clear();
  V.clear();

  // No truncation of T or V was done so far... and not now!

  // Write to file and read again to reduce memory fragmentation.
  H_core_Matrix_sparse.writeToFile();
  H_core_Matrix_sparse.readFromFile();

  if(create_dipole_mtx == 1) {
    // Output dipole matrices in mtx format
    if ( calculation_identifier == 0 )
      throw "calculation_identifier == 0 when create_dipole_mtx == 1";
    if ( method_and_basis_set == 0 )
      throw "method_and_basis_set == 0 when create_dipole_mtx == 1";
    if ( inversePermutationHML  == 0 )
      throw "inversePermutationHML == 0 when create_dipole_mtx == 1";	  
    {
      // dipole matrix X
      std::stringstream id_ss;
      id_ss << *calculation_identifier << " - dipole matrix X";
      symmMatrix fieldMatrix(H_core_Matrix_sparse);
      fieldMatrix.clear();
      if(compute_operator_matrix_sparse_symm(basisInfo, 
					     1, 0, 0, 
					     fieldMatrix, 
					     permutationHML) != 0)
	return -1;
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Creating mtx file for dipole matrix X");
      write_matrix_in_matrix_market_format( fieldMatrix, 
					    *inversePermutationHML, 
					    "dipole_matrix_x", 
					    id_ss.str(), 
					    *method_and_basis_set );	
    }
    {
      // dipole matrix Y
      std::stringstream id_ss;
      id_ss << *calculation_identifier << " - dipole matrix Y";
      symmMatrix fieldMatrix(H_core_Matrix_sparse);
      fieldMatrix.clear();
      if(compute_operator_matrix_sparse_symm(basisInfo, 
					     0, 1, 0, 
					     fieldMatrix, 
					     permutationHML) != 0)
	return -1;
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Creating mtx file for dipole matrix Y");
      write_matrix_in_matrix_market_format( fieldMatrix, 
					    *inversePermutationHML, 
					    "dipole_matrix_y", 
					    id_ss.str(), 
					    *method_and_basis_set );	
    }
    {
      // dipole matrix Z
      std::stringstream id_ss;
      id_ss << *calculation_identifier << " - dipole matrix Z";
      symmMatrix fieldMatrix(H_core_Matrix_sparse);
      fieldMatrix.clear();
      if(compute_operator_matrix_sparse_symm(basisInfo, 
					     0, 0, 1, 
					     fieldMatrix, 
					     permutationHML) != 0)
	return -1;
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Creating mtx file for dipole matrix Z");
      write_matrix_in_matrix_market_format( fieldMatrix, 
					    *inversePermutationHML, 
					    "dipole_matrix_z", 
					    id_ss.str(), 
					    *method_and_basis_set );	
    }
  } // end if(create_dipole_mtx == 1)

  
  // Add contributions from electric field, if any
  if(electric_field_x != 0)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		"electric_field_x = %22.11f", (double)electric_field_x);
      symmMatrix fieldMatrix(H_core_Matrix_sparse);
      fieldMatrix.clear();
      if(compute_operator_matrix_sparse_symm(basisInfo, 
					     1, 0, 0, 
					     fieldMatrix, 
					     permutationHML) != 0)
	return -1;
      H_core_Matrix_sparse += 
	(ergo_real)(-1.0) * electric_field_x * fieldMatrix;
    }
  if(electric_field_y != 0)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		"electric_field_y = %22.11f", (double)electric_field_y);
      symmMatrix fieldMatrix(H_core_Matrix_sparse);
      fieldMatrix.clear();
      if(compute_operator_matrix_sparse_symm(basisInfo, 
					     0, 1, 0, 
					     fieldMatrix, 
					     permutationHML) != 0)
	return -1;
      H_core_Matrix_sparse += 
	(ergo_real)(-1.0) * electric_field_y * fieldMatrix;
    }
  if(electric_field_z != 0)
    {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
		"electric_field_z = %22.11f", (double)electric_field_z);
      symmMatrix fieldMatrix(H_core_Matrix_sparse);
      fieldMatrix.clear();
      if(compute_operator_matrix_sparse_symm(basisInfo, 
					     0, 0, 1, 
					     fieldMatrix, 
					     permutationHML) != 0)
	return -1;
      H_core_Matrix_sparse += 
	(ergo_real)(-1.0) * electric_field_z * fieldMatrix;
    }

  // Add contribution from extra charges, if any.
  if(extraCharges.getNoOfAtoms() > 0) {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Using 'extra charges', extraCharges.noOfAtoms = %d", extraCharges.getNoOfAtoms());
    ergo_real sumCharge = 0;
    ergo_real maxCharge = extraCharges.getAtom(0).charge;
    ergo_real minCharge = extraCharges.getAtom(0).charge;
    for(int i = 0; i < extraCharges.getNoOfAtoms(); i++) {
      ergo_real currCharge = extraCharges.getAtom(i).charge;
      sumCharge += currCharge;
      if(currCharge > maxCharge)
	maxCharge = currCharge;
      if(currCharge < minCharge)
	minCharge = currCharge;
    }
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Using 'extra charges', sum = %8.4f, min = %8.4f, max = %8.4f", 
	      (double)sumCharge, (double)minCharge, (double)maxCharge);
    ergo_real minDist = template_blas_get_num_limit_max<ergo_real>();
    ergo_real minDistCharge1 = 0;
    ergo_real minDistCharge2 = 0;
    for(int i = 0; i < extraCharges.getNoOfAtoms(); i++)
      for(int j = 0; j < molecule.getNoOfAtoms(); j++) {
	ergo_real dx = extraCharges.getAtom(i).coords[0] - molecule.getAtom(j).coords[0];
	ergo_real dy = extraCharges.getAtom(i).coords[1] - molecule.getAtom(j).coords[1];
	ergo_real dz = extraCharges.getAtom(i).coords[2] - molecule.getAtom(j).coords[2];
	ergo_real dist = template_blas_sqrt(dx*dx+dy*dy+dz*dz);
	if(dist < minDist) {
	  minDist = dist;
	  minDistCharge1 = molecule.getAtom(j).charge;
	  minDistCharge2 = extraCharges.getAtom(i).charge;
	}
      }
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Min dist between any 'extra charge' atom and any atom in main molecule: %9.5f a.u. = %9.5f Angstrom", 
	      (ergo_real)minDist, (ergo_real)(minDist/UNIT_one_Angstrom));
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Charges for min dist: %9.5f (main) and %9.5f (extra)", (double)minDistCharge1, (double)minDistCharge2);
    Util::TimeMeter timeMeter_V_extra;
    symmMatrix V_extra;
    V_extra.resetSizesAndBlocks(matrix_size_block_info,
				matrix_size_block_info);
    ergo_real nuclearRepulsionEnergy_dummy = 0;
    if(compute_V_sparse(basisInfo,
			integralInfo, 
			extraCharges,
			threshold_integrals_1el,
			boxSizeForVT,
			V_extra,
			permutationHML,
			nuclearRepulsionEnergy_dummy) != 0) {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_V_sparse");
      return -1;
    }    
    timeMeter_V_extra.print(LOG_AREA_SCF, "Computation of V_extra");
    H_core_Matrix_sparse += (ergo_real)(1.0) * V_extra;
  }

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "compute_h_core_matrix_sparse ending OK.");
  timeMeterTot.print(LOG_AREA_SCF, "compute_h_core_matrix_sparse");
  
  return 0;
}


int
get_gradient_for_given_mol_and_dens(const IntegralInfo& integralInfo,
				    const Molecule& molecule,
				    const BasisInfoStruct& basisInfo,
				    const symmMatrix & D,
				    ergo_real threshold_integrals_1el,
				    mat::SizesAndBlocks const & matrix_size_block_info,
				    std::vector<int> const & permutationHML,
				    ergo_real* result_gradient_list)
{
  ergo_real boxSize = 1.0;
  if(compute_gradient_of_nucl_and_trDV(basisInfo,
				       integralInfo, 
				       molecule,
				       threshold_integrals_1el,
				       boxSize,
				       D,
				       permutationHML,
				       result_gradient_list) != 0) {
    throw "Error in compute_gradient_of_tr_D_V";
  }
  return 0;
}



/** Saves specified symmetic matrix to a file of specified name.
    @param A the matrix to save. The matrix must be saved to a backing
    store already.
    @param basisInfo the basis set description.
    @param fileName The file that will contain the matrix. It will be
    overwritten without further questions.
    @param inversePermutationHML permutation information needed when using hierarchic matrices.
*/
int
save_symmetric_matrix(symmMatrix& A, const BasisInfoStruct & basisInfo,
                      const char *fileName,
		      std::vector<int> const & inversePermutationHML)
{
  int n = basisInfo.noOfBasisFuncs;
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "save_symmetric_matrix allocating full n*n matrix.");
  ergo_real* potentialMatrixFull = new ergo_real[n*n];
  normalMatrix* tmpMat;
  A.readFromFile();
  tmpMat = new normalMatrix(A);
  A.writeToFile();
  std::vector<ergo_real> fullTmp(n*n);
  tmpMat->fullMatrix(fullTmp,
		     inversePermutationHML,
		     inversePermutationHML);
  std::copy (fullTmp.begin(), fullTmp.end(), potentialMatrixFull);
  delete tmpMat;
  int ret = ddf_writeShellListAndDensityMatricesToFile(basisInfo,
                                                       1,
                                                       &potentialMatrixFull,
                                                       fileName);
  delete [] potentialMatrixFull;
  return ret;
}

int 
add_disturbance_to_matrix(int n, 
			  symmMatrix & A,
			  ergo_real disturbance,
			  int specificElementCount,
			  const int* elementIndexVector,
			  std::vector<int> const & permutationHML)
{
  // Create diagonal matrix with random numbers.
  std::vector<ergo_real> diagonalValues(n);
  std::vector<int> rowind(n);
  std::vector<int> colind(n);
  for(int i = 0; i < n; i++)
    diagonalValues[i] = 0;
  if(specificElementCount == 0)
    {
      // use random numbers
      for(int i = 0; i < n; i++)
	{
	  double randomNumber = (double)rand() / RAND_MAX;
	  // Now randomNumber is between 0 and 1
	  randomNumber *= 2;
	  // Now randomNumber is between 0 and 2
	  randomNumber -= 1;
	  // Now randomNumber is between -1 and 1
	  diagonalValues[i] = randomNumber * disturbance;
	}      
    }
  else
    {
      // use list of specific elements
      for(int i = 0; i < specificElementCount; i++)
	{
	  int idx = elementIndexVector[i];
	  if(idx >= 0 && idx < n)
	    diagonalValues[idx] = disturbance;
	}
    }
  for(int i = 0; i < n; i++)
    {
      rowind[i] = i;
      colind[i] = i;
    }
  symmMatrix B(A);
  B.clear();
  B.assign_from_sparse(rowind,
		       colind,
		       diagonalValues,
		       permutationHML,
		       permutationHML);
  A += (ergo_real)1.0 * B;
  
  return 0;
}



int
get_simple_starting_guess_sparse(int n, 
				 int noOfElectrons, 
				 symmMatrix & densityMatrix)
{
  std::vector<ergo_real> diagonalValues(n);
  std::vector<int> rowind(n);
  std::vector<int> colind(n);
  for(int i = 0; i < n; i++)
    diagonalValues[i] = (ergo_real)noOfElectrons / n;
  for(int i = 0; i < n; i++)
    {
      rowind[i] = i;
      colind[i] = i;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_simple_starting_guess_sparse: calling densityMatrix.assign_from_sparse.");
  densityMatrix.assign_from_sparse(rowind, colind, diagonalValues);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_simple_starting_guess_sparse: densityMatrix.assign_from_sparse returned.");
  return 0;
}


int
get_diag_matrix_from_file(int n, 
			  symmMatrix & M, 
			  const char* fileName,
			  std::vector<int> const & permutationHML)
{
  std::vector<ergo_real> diagonalValues(n);
  std::vector<int> rowind(n);
  std::vector<int> colind(n);
  // First read values from file.
  FILE* f = fopen(fileName, "rt");
  if(f == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error opening file '%s'", fileName);
      return -1;
    }

  const int BUFSIZE = 888888;
  std::vector<char> buf(BUFSIZE);
  memset(&buf[0], 0, BUFSIZE * sizeof(char));

  int nBytes = (int)fread(&buf[0], sizeof(char), BUFSIZE, f);
  fclose(f);
  if(nBytes <= 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error reading file '%s'", fileName);
      return -1;
    }
  if(nBytes >= (BUFSIZE-1000))
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in get_diag_matrix_from_file: input file too large");
      return -1;
    }

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "get_diag_matrix_from_file: File '%s' read OK, nBytes = %i", fileName, nBytes);

  int i = 0;
  char* p = &buf[0];
  while(*p != '\0')
    {
      // Now we expect p to point to the beginning of a new line
      if(i >= n)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error in get_diag_matrix_from_file: too many lines in file.");
	  return -1;
	}

      // Skip blanks
      while(*p == ' ' || *p == '\t')
	p++;

      diagonalValues[i] = atof(p);
      i++;

      // Skip rest of line
      while(*p != '\n' && *p != '\0')
	p++;
      p++;
    } // END WHILE

  for(i = 0; i < n; i++)
    {
      rowind[i] = i;
      colind[i] = i;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_diag_matrix_from_file: calling densityMatrix.assign_from_sparse.");
  M.assign_from_sparse(rowind,
		       colind,
		       diagonalValues,
		       permutationHML,
		       permutationHML);
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_diag_matrix_from_file: densityMatrix.assign_from_sparse returned.");
  
  return 0;
}


int
write_diag_elements_to_file(int n, 
			    const symmMatrix & M, 
			    const char* fileName,
			    std::vector<int> const & permutationHML)
{
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "write_diag_elements_to_file, n = %i", n);

  FILE* f = fopen(fileName, "wt");
  if(f == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_MAIN, "error opening file '%s'.", fileName);
      return -1;
    }
  
  std::vector<int> rowind(n);
  std::vector<int> colind(n);
  std::vector<ergo_real> values(n);

  for(int i = 0; i < n; i++)
    {
      rowind[i] = i;
      colind[i] = i;
    }

  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "write_diag_elements_to_file, calling M.get_values.");
  M.get_values(rowind,
	       colind,
	       values,
	       permutationHML,
	       permutationHML);
  do_output(LOG_CAT_INFO, LOG_AREA_MAIN, "write_diag_elements_to_file, after M.get_values.");

  for(int i = 0; i < n; i++)
    fprintf(f, "%12.6f\n", (double)values[i]);
  fclose(f);
  
  return 0;
}



int 
write_full_matrix(int n, const symmMatrix & M, const char* fileName,
		  std::vector<int> const & inversePermutationHML)
{
  std::vector<ergo_real> buf(n*n);
  M.fullMatrix(buf, inversePermutationHML, inversePermutationHML);
  FILE* f = fopen(fileName, "wb");
  if(f == NULL)
    throw "error in write_full_matrix, in fopen.";
  size_t noOfItemsWritten = fwrite(&buf[0], sizeof(ergo_real), n*n, f);
  fclose(f);  
  if(noOfItemsWritten != (size_t)n*n)
    throw "error in write_full_matrix, in fwrite.";
  return 0;
}


int
write_basis_func_coord_file(const BasisInfoStruct & basisInfo)
{
  FILE* f = fopen("basis_func_coords.m", "wt");
  if(f == NULL)
    throw "write_basis_func_coord_file, in fopen.";
  int n = basisInfo.noOfBasisFuncs;
  fprintf(f, "basis_func_coords = [\n");
  for(int i = 0; i < n; i++)
    fprintf(f, "%15.9f %15.9f %15.9f\n", 
	    (double)basisInfo.basisFuncList[i].centerCoords[0], 
	    (double)basisInfo.basisFuncList[i].centerCoords[1], 
	    (double)basisInfo.basisFuncList[i].centerCoords[2]);
  fprintf(f, "];\n");
  fclose(f);
  return 0;
}


int
write_2el_integral_m_file(const BasisInfoStruct & basisInfo, const IntegralInfo & integralInfo)
{
  int n = basisInfo.noOfBasisFuncs;
  FILE* f = fopen("setup_twoel_integrals.m", "wt");
  if(f == NULL)
    throw "write_2el_integral_m_file, in fopen.";  
  fprintf(f, "twoel_integrals = zeros(%2d, %2d, %2d, %2d);\n", n, n, n, n);
  for(int a = 0; a < n; a++)
    for(int b = 0; b < n; b++)
      for(int c = 0; c < n; c++)
	for(int d = 0; d < n; d++) {
	  ergo_real integralValue = do_2e_integral(a, b, c, d, basisInfo, integralInfo);
	  fprintf(f, "twoel_integrals(%2d, %2d, %2d, %2d) = %15.11f;\n", a+1, b+1, c+1, d+1, (double)integralValue);
	}
  fclose(f);
  return 0;
}



static int
get_2e_matrix_and_energy_simple_HF_sparse(const BasisInfoStruct& basisInfo,
					  const Molecule& molecule,
					  const IntegralInfo& integralInfo, 
					  const JK::ExchWeights & CAM_params,
					  symmMatrix & twoelMatrix_sparse, 
					  symmMatrix & densityMatrix_sparse, 
					  const JK::Params& J_K_params,
					  mat::SizesAndBlocks const & matrix_size_block_info,
					  std::vector<int> const & permutationHML,
					  std::vector<int> const & inversePermutationHML,
					  ergo_real* energy_2el,
					  int get_J_K_matrices,
					  symmMatrix & J_matrix,
					  symmMatrix & K_matrix)
{

  symmMatrix J;
  J.resetSizesAndBlocks(matrix_size_block_info,
			       matrix_size_block_info);
  symmMatrix K;
  K.resetSizesAndBlocks(matrix_size_block_info,
			       matrix_size_block_info);
  // FMM for J
  if(compute_J_by_boxes_sparse(basisInfo,
			       integralInfo,
			       J_K_params,
			       J,
			       densityMatrix_sparse,
			       permutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_J_by_boxes_sparse");
      return -1;
    }
  
  // exchange matrix K
  if(compute_K_by_boxes_sparse(basisInfo,
			       integralInfo,
			       CAM_params,
			       J_K_params,
			       K,
			       densityMatrix_sparse,
			       permutationHML,
			       inversePermutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes_sparse");
      return -1;
    }
  
  twoelMatrix_sparse = J + K;
  if ( get_J_K_matrices == 1 ) {
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Saving J and K matrices as requested.");
    J_matrix = J;
    K_matrix = K;
  }

  *energy_2el = 0.5 * symmMatrix::trace_ab(densityMatrix_sparse, twoelMatrix_sparse);
  
  return 0;
}



static int
get_2e_matrices_and_energy_simple_HF_sparse_unrestricted(const BasisInfoStruct& basisInfo,
							 const Molecule& molecule,
							 const IntegralInfo& integralInfo, 
							 symmMatrix & twoelMatrix_sparse_alpha, 
							 symmMatrix & twoelMatrix_sparse_beta, 
							 symmMatrix & densityMatrix_sparse_alpha, 
							 symmMatrix & densityMatrix_sparse_beta, 
							 const JK::Params& J_K_params,
							 const JK::ExchWeights & CAM_params,
							 mat::SizesAndBlocks const & matrix_size_block_info,
							 std::vector<int> const & permutationHML,
							 std::vector<int> const & inversePermutationHML,
							 ergo_real* energy_2el)
{
  symmMatrix J;
  J.resetSizesAndBlocks(matrix_size_block_info,
			       matrix_size_block_info);
  symmMatrix K_alpha;
  K_alpha.resetSizesAndBlocks(matrix_size_block_info,
				     matrix_size_block_info);
  symmMatrix K_beta;
  K_beta.resetSizesAndBlocks(matrix_size_block_info,
				    matrix_size_block_info);
  // Compute total density matrix
  symmMatrix D_tot(densityMatrix_sparse_alpha);
  D_tot += densityMatrix_sparse_beta;

  // FMM for J
  if(compute_J_by_boxes_sparse(basisInfo,
			       integralInfo,
			       J_K_params,
			       J,
			       D_tot,
			       permutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_J_by_boxes_sparse");
      return -1;
    }
  D_tot.clear();
  
  // exchange matrices
  if(compute_K_by_boxes_sparse(basisInfo,
			       integralInfo,
			       CAM_params,
			       J_K_params,
			       K_alpha,
			       densityMatrix_sparse_alpha,
			       permutationHML,
			       inversePermutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes_sparse");
      return -1;
    }
  if(compute_K_by_boxes_sparse(basisInfo,
			       integralInfo,
			       CAM_params,
			       J_K_params,
			       K_beta,
			       densityMatrix_sparse_beta,
			       permutationHML,
			       inversePermutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes_sparse");
      return -1;
    }
 
  K_alpha *= 2;
  K_beta  *= 2;
  twoelMatrix_sparse_alpha = J + K_alpha;
  twoelMatrix_sparse_beta  = J + K_beta;

  *energy_2el = 0;
  *energy_2el +=  0.5 * symmMatrix::trace_ab(densityMatrix_sparse_alpha, twoelMatrix_sparse_alpha);
  *energy_2el +=  0.5 * symmMatrix::trace_ab(densityMatrix_sparse_beta , twoelMatrix_sparse_beta );
  
  return 0;
}

/** 
Routine for computing the two-electron part of the Fock/KS matrix, in
sparse form. This routine is only used when FMM is enabled.
*/
static int
get_2e_matrix_and_energy_simple_sparse(const BasisInfoStruct& basisInfo,
				       const Molecule& molecule,
				       const IntegralInfo& integralInfo, 
				       symmMatrix & twoelMatrix_sparse, 
				       symmMatrix & densityMatrix_sparse, // should be in memory when this routine is called
				       const JK::Params& J_K_params,
				       const JK::ExchWeights & CAM_params,
				       const Dft::GridParams& gridParams,
				       mat::SizesAndBlocks const & matrix_size_block_info,
				       std::vector<int> const & permutationHML,
				       std::vector<int> const & inversePermutationHML,
				       ergo_real* energy_2el,
				       int do_xc,
				       int noOfElectrons,
				       int get_J_K_Fxc_matrices,
				       symmMatrix & J_matrix,
				       symmMatrix & K_matrix,
				       symmMatrix & Fxc_matrix,
				       SCF_statistics & stats)
{
  {
    // FMM for J
    symmMatrix J;
    J.resetSizesAndBlocks(matrix_size_block_info,
				 matrix_size_block_info);
    stats.start_timer("compute_J_by_boxes_sparse");
    if(compute_J_by_boxes_sparse(basisInfo,
				 integralInfo,
				 J_K_params,
				 J,
				 densityMatrix_sparse,
				 permutationHML) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_J_by_boxes_sparse");
	return -1;
      }
    stats.stop_timer("compute_J_by_boxes_sparse");
    twoelMatrix_sparse = J;
    if ( get_J_K_Fxc_matrices == 1 ) {
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Saving J matrix as requested.");
      J_matrix = J;
    }
  }

  
  if(CAM_params.alpha != 0.0 || CAM_params.computeRangeSeparatedExchange)
    {
      // exchange matrix K
      symmMatrix K;
      K.resetSizesAndBlocks(matrix_size_block_info,
				   matrix_size_block_info);
      stats.start_timer("compute_K_by_boxes_sparse");
      if(compute_K_by_boxes_sparse(basisInfo,
				   integralInfo,
				   CAM_params,
				   J_K_params,
				   K,
				   densityMatrix_sparse,
				   permutationHML,
				   inversePermutationHML) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes_sparse");
	  return -1;
	}
      stats.stop_timer("compute_K_by_boxes_sparse");
      if(!CAM_params.computeRangeSeparatedExchange) {
        twoelMatrix_sparse += CAM_params.alpha * K;
	if ( get_J_K_Fxc_matrices == 1 ) {
	  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Saving K matrix as requested.");
	  K_matrix = CAM_params.alpha * K;
	}
      }
      else {
	twoelMatrix_sparse += K;
	if ( get_J_K_Fxc_matrices == 1 ) {
	  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Saving K matrix as requested.");
	  K_matrix = K;
	}
      }
    }

  *energy_2el = 0.5 * symmMatrix::trace_ab(densityMatrix_sparse, twoelMatrix_sparse);
  if(do_xc)
    {
      ergo_real dftEnergy;
      symmMatrix F_xc;
      // construct matrix Fxc and xc_energy
      
      output_current_memory_usage(LOG_AREA_SCF, "before dft_get_xc_mt");

      if(do_xc == 1) {
        int n = basisInfo.noOfBasisFuncs;
        do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_2e_matrix_and_energy_simple_sparse (do_xc == 1) allocating full n*n matrix.");
        ergo_real* densityMatrix_full = new ergo_real[n*n];
	{
	  std::vector<ergo_real> fullTmp(n*n);
	  densityMatrix_sparse.fullMatrix(fullTmp,
					  inversePermutationHML,
					  inversePermutationHML);
	  std::copy (fullTmp.begin(), fullTmp.end(), densityMatrix_full);
	}
	std::vector<ergo_real> F_xc_full(n*n, 0);
	
        // dft_get_xc_mt is the call to the basic DFT code
        do_output(LOG_CAT_INFO, LOG_AREA_SCF,
                  "calling dft_get_xc_mt, n = %i, noOfElectrons = %i",
                  n, noOfElectrons);

        dft_get_xc_mt(noOfElectrons, densityMatrix_full,
                      &basisInfo, &molecule, gridParams,
                      &F_xc_full[0], &dftEnergy);
      
        output_current_memory_usage(LOG_AREA_SCF, "after  dft_get_xc_mt");
        delete []densityMatrix_full;
      	F_xc.resetSizesAndBlocks(matrix_size_block_info,
					matrix_size_block_info);
	F_xc.assignFromFull(F_xc_full, permutationHML, permutationHML);
      }
      else {
        do_output(LOG_CAT_INFO, LOG_AREA_SCF,
                  "calling sparse Dft::getXC, noOfElectrons = %i",
                  noOfElectrons);
      	F_xc.resetSizesAndBlocks(matrix_size_block_info,
					matrix_size_block_info);

	/* FIXME: what does the comment below mean??  */
        /* Threaded version gathers data in a thread-unsafe way. */
	
	stats.start_timer("getXC_mt");
        Dft::getXC_mt(basisInfo, integralInfo,
                      molecule, gridParams, noOfElectrons,
                      densityMatrix_sparse, F_xc, &dftEnergy,
		      permutationHML);
	stats.stop_timer("getXC_mt");
      
        output_current_memory_usage(LOG_AREA_SCF, "after  dft_get_xc_mt");
      }      
      twoelMatrix_sparse += F_xc;
      if ( get_J_K_Fxc_matrices == 1 ) {
	do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Saving Fxc matrix as requested.");
	Fxc_matrix = F_xc;      
      }
      *energy_2el += dftEnergy;
    }
  
  return 0;
}






static int
get_2e_matrices_and_energy_simple_sparse_unrestricted(const BasisInfoStruct& basisInfo,
                                       const Molecule& molecule,
                                       const IntegralInfo& integralInfo,
				       const JK::ExchWeights & CAM_params,
                                       symmMatrix & twoelMatrix_sparse_alpha,
                                       symmMatrix & twoelMatrix_sparse_beta,
                                       symmMatrix & densityMatrix_sparse_alpha, 
                                       symmMatrix & densityMatrix_sparse_beta,
                                       const JK::Params& J_K_params,
				       const Dft::GridParams& gridParams,
                                       mat::SizesAndBlocks const & matrix_size_block_info,
				       std::vector<int> const & permutationHML,
				       std::vector<int> const & inversePermutationHML,
                                       ergo_real* energy_2el,
                                       int do_xc,
                                       int noOfElectrons)
{
  symmMatrix J;
  J.resetSizesAndBlocks(matrix_size_block_info,
			matrix_size_block_info);
  // Compute total density matrix
  symmMatrix D_tot(densityMatrix_sparse_alpha);
  D_tot += densityMatrix_sparse_beta;

  // FMM for J
  if(compute_J_by_boxes_sparse(basisInfo,
			       integralInfo,
			       J_K_params,
			       J,
			       D_tot,
			       permutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_J_by_boxes_sparse");
      return -1;
    }
  twoelMatrix_sparse_alpha = J;
  twoelMatrix_sparse_beta = J;
  J.clear();

  if(CAM_params.alpha != 0 || CAM_params.computeRangeSeparatedExchange)
    {
      D_tot *= ergo_real(0.5);
      symmMatrix D_spin(ergo_real(0.5)*densityMatrix_sparse_alpha);
      D_spin -= ergo_real(0.5)*densityMatrix_sparse_beta;

      /* Multiply by alpha only if no range-separated exchange is computed.
	 It CAM params are used the hf weight has already been accounted for. */
      ergo_real real_hf_weight = (CAM_params.computeRangeSeparatedExchange == 0)
	? CAM_params.alpha * 2.0 : 2.0;
      symmMatrix K_tot, K_spin;
      K_tot.resetSizesAndBlocks(matrix_size_block_info,
				matrix_size_block_info);
      K_spin.resetSizesAndBlocks(matrix_size_block_info,
				 matrix_size_block_info);

      // exchange matrices
      if(compute_K_by_boxes_sparse(basisInfo,
				   integralInfo,
				   CAM_params,
				   J_K_params,
				   K_tot,
				   D_tot,
				   permutationHML,
				   inversePermutationHML) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes_sparse");
	  return -1;
	}
      twoelMatrix_sparse_alpha += real_hf_weight * K_tot;
      twoelMatrix_sparse_beta  += real_hf_weight * K_tot;
      K_tot.clear();

      if(compute_K_by_boxes_sparse(basisInfo,
				   integralInfo,
				   CAM_params,
				   J_K_params,
				   K_spin,
				   D_spin,
				   permutationHML,
				   inversePermutationHML) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes_sparse");
	  return -1;
	}

      twoelMatrix_sparse_alpha += real_hf_weight * K_spin;
      twoelMatrix_sparse_beta  -= real_hf_weight * K_spin;
    }
  D_tot.clear();
  
  *energy_2el = 0;
  *energy_2el +=  0.5 * symmMatrix::trace_ab(densityMatrix_sparse_alpha, twoelMatrix_sparse_alpha);
  *energy_2el +=  0.5 * symmMatrix::trace_ab(densityMatrix_sparse_beta , twoelMatrix_sparse_beta );

  if(do_xc == 1)
    {
      ergo_real dftEnergy;
      int n = basisInfo.noOfBasisFuncs;
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_2e_matrices_and_energy_simple_sparse_unrestricted allocating full n*n matrices.");
      std::vector<ergo_real> xcMatrix_full_alpha(n*n,0);
      std::vector<ergo_real> xcMatrix_full_beta (n*n,0);
      std::vector<ergo_real> densityMatrix_full_alpha(n*n);
      std::vector<ergo_real> densityMatrix_full_beta (n*n);
      densityMatrix_sparse_alpha.fullMatrix(densityMatrix_full_alpha,
					    inversePermutationHML,
					    inversePermutationHML);
      densityMatrix_sparse_beta .fullMatrix(densityMatrix_full_beta,
					    inversePermutationHML,
					    inversePermutationHML);
      
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "calling dft_get_uxc_mt, n = %i", n);
      

#ifndef SKIP_DFT_FLAG
      dft_get_uxc_mt(noOfElectrons,
                     &densityMatrix_full_alpha[0], &densityMatrix_full_beta[0],
		     &basisInfo, &molecule, gridParams,
		     &xcMatrix_full_alpha[0], &xcMatrix_full_beta[0], &dftEnergy);
#endif

      densityMatrix_full_alpha.reserve(0);
      densityMatrix_full_beta.reserve(0);

      *energy_2el += dftEnergy;

      output_current_memory_usage(LOG_AREA_SCF, "after dft_get_uxc_mt");

      symmMatrix tmp;
      tmp.resetSizesAndBlocks(matrix_size_block_info,
				     matrix_size_block_info);
      tmp.assignFromFull(xcMatrix_full_alpha,
			 permutationHML,
			 permutationHML);
      twoelMatrix_sparse_alpha += tmp;
      tmp.assignFromFull(xcMatrix_full_beta,
			 permutationHML,
			 permutationHML);
      twoelMatrix_sparse_beta  += tmp;
      
    } 
  else
    {
      ergo_real dftEnergy = 0;
      symmMatrix FA_xc, FB_xc;
      do_output(LOG_CAT_INFO, LOG_AREA_SCF,
                "calling sparse Dft::getUXC, noOfElectrons = %i",
                noOfElectrons);
      FA_xc.resetSizesAndBlocks(matrix_size_block_info,
                                matrix_size_block_info);
      FB_xc.resetSizesAndBlocks(matrix_size_block_info,
                                matrix_size_block_info);

      /* FIXME: what does the comment below mean??  */
      /* Threaded version gathers data in a thread-unsafe way. */

      Dft::getUXC_seq(basisInfo, integralInfo,
                     molecule, gridParams, noOfElectrons,
                     densityMatrix_sparse_alpha, densityMatrix_sparse_beta,
                     FA_xc, FB_xc, &dftEnergy,
                     permutationHML);
      
      output_current_memory_usage(LOG_AREA_SCF, "after  dft_get_xc_mt");
      twoelMatrix_sparse_alpha += FA_xc;
      twoelMatrix_sparse_beta += FB_xc;
      *energy_2el += dftEnergy;
    }

  return 0;
}



/** 
General routine for computing the two-electron part of the Fock/KS
matrix.  Both input and output matrices are in sparse form, but full
matrix form may be used in intermediate steps. If FMM is not used,
full matrix format is applied.

@param basisInfo the used basis set.
@param molecule position of atoms (used for eg. XC grid).
@param integralInfo - the integral evaluation object.
@param twoelMatrix_sparse - the evaluation result.

@param densityMatrix_sparse - the density for which 2e matrix is to be
evaluated.

@param J_K_params the settings of the integral evaluation.

@param do_xc whether xc contribution to 2e matrix and energy are to be
added. 1 means that the traditional full matrix code should be called,
2 means that the sparse variant is to be used.

@param energy_2el 2el energy contribution

@param noOfElectrons number of electrons...

@param CAM_params a structure containing parameters needed when functionals like CAMB3LYP are used.
@param gridParams a structure containing parameters for the grid.
@param matrix_size_block_info block sizes etc for hierarchic matrix library.

@param get_J_K_Fxc_matrices flag saying if matrices should be saved for statistics/testing purposes.
If that feature is active, matrices are saved in parameters @param J_matrix @param K_matrix @param Fxc_matrix .
@param stats a structure holding SCF statistics.

@param permutationHML - the permutation of basis functions, needed
for transformation between the dense and sparse formats.
@param inversePermutationHML - the inverse permutation of basis functions, needed
for transformation between the dense and sparse formats.
*/
int
get_2e_matrix_and_energy_sparse(const BasisInfoStruct & basisInfo,
				const Molecule& molecule,
				const IntegralInfo& integralInfo, 
				symmMatrix & twoelMatrix_sparse, 
				symmMatrix & densityMatrix_sparse,
				const JK::Params& J_K_params,
				const JK::ExchWeights & CAM_params,
				const Dft::GridParams& gridParams,
				int do_xc,
				ergo_real* energy_2el,
				int noOfElectrons,
				mat::SizesAndBlocks const & matrix_size_block_info,
				std::vector<int> const & permutationHML,
				std::vector<int> const & inversePermutationHML,
				int get_J_K_Fxc_matrices,
				symmMatrix & J_matrix,
				symmMatrix & K_matrix,
				symmMatrix & Fxc_matrix,
				SCF_statistics & stats)
{
  if(J_K_params.use_naive_fockmat_constr == 1)
    {
      int n = basisInfo.noOfBasisFuncs;
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_2e_matrix_and_energy_sparse naive impl allocating full n*n matrices.");
      ergo_real* twoelMatrixG      = new ergo_real[n*n];
      ergo_real* densityMatrixFull = new ergo_real[n*n];
      {
	std::vector<ergo_real> fullTmp(n*n);
	densityMatrix_sparse.fullMatrix(fullTmp,
					inversePermutationHML,
					inversePermutationHML);
	std::copy (fullTmp.begin(), fullTmp.end(), densityMatrixFull);
      }
      // Use smallest of threshold_J and threshold_K.
      ergo_real threshold_JK = J_K_params.threshold_J;
      if(J_K_params.threshold_K < threshold_JK)
	threshold_JK = J_K_params.threshold_K;
      if(compute_2e_matrix_list_explicit(basisInfo, integralInfo, &twoelMatrixG, &densityMatrixFull, 1, threshold_JK) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_2e_matrix_list_explicit");
	  return -1;
	}
      *energy_2el = get_trace(n, densityMatrixFull, twoelMatrixG)*0.5;
      {
	std::vector<ergo_real> fullTmp(twoelMatrixG, twoelMatrixG + n*n);
	twoelMatrix_sparse.assignFromFull(fullTmp,
					  permutationHML,
					  permutationHML);
      }
      delete []twoelMatrixG;
      delete []densityMatrixFull;
      return 0;
    }
  
  if(do_xc == 0 && J_K_params.use_fmm == 1)
    {
      // Simple Hartree-Fock case, with FMM. 
      // In this case we use a special function to get away with less memory usage.
      if(get_2e_matrix_and_energy_simple_HF_sparse(basisInfo, 
						   molecule,
						   integralInfo, 
						   CAM_params,
						   twoelMatrix_sparse, 
						   densityMatrix_sparse, 
						   J_K_params,
						   matrix_size_block_info,
						   permutationHML,
						   inversePermutationHML,
						   energy_2el,
						   get_J_K_Fxc_matrices,
						   J_matrix,
						   K_matrix) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_2e_matrix_and_energy_simple_HF_sparse");
	  return -1;
	}
      return 0;
    }

  if(J_K_params.use_fmm == 1 || CAM_params.computeRangeSeparatedExchange)
    {
      // Simple case, with FMM. 
      if(get_2e_matrix_and_energy_simple_sparse(basisInfo,
						molecule,
						integralInfo, 
						twoelMatrix_sparse, 
						densityMatrix_sparse, 
						J_K_params,
						CAM_params,
						gridParams,
						matrix_size_block_info,
						permutationHML,
						inversePermutationHML,
						energy_2el,
						do_xc,
						noOfElectrons,
						get_J_K_Fxc_matrices,
						J_matrix,
						K_matrix,
						Fxc_matrix,
						stats) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_2e_matrix_and_energy_simple_sparse");
	  return -1;
	}
      return 0;
    }

  int n = basisInfo.noOfBasisFuncs;
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_2e_matrix_and_energy_sparse allocating full n*n matrices.");
  ergo_real* twoelMatrix_full   = new ergo_real[n*n];
  ergo_real* densityMatrix_full = new ergo_real[n*n];
  {
    std::vector<ergo_real> fullTmp(n*n);
    densityMatrix_sparse.fullMatrix(fullTmp,
				    inversePermutationHML,
				    inversePermutationHML);
    std::copy (fullTmp.begin(), fullTmp.end(), densityMatrix_full);
  }
  {
    // Get J
    if(compute_2e_matrix_coulomb(basisInfo,
				 integralInfo, 
				 twoelMatrix_full, 
				 densityMatrix_full, 
				 J_K_params) != 0)
      {
	do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_2e_matrix_coulomb");
	return -1;
      }
    if(CAM_params.alpha != ergo_real(0.0) || CAM_params.computeRangeSeparatedExchange)
      {
	ergo_real* K   = new ergo_real[n*n];
	// Get K
	if(compute_2e_matrix_exchange(basisInfo, 
				      integralInfo, 
				      CAM_params,
				      K, 
				      densityMatrix_full,
				      J_K_params.threshold_K) != 0)
	  {
	    do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_2e_matrix_exchange");
	    return -1;
	  }
	add_square_matrices(n, twoelMatrix_full, K, twoelMatrix_full);
	delete []K;
      }
  }
  
  *energy_2el = get_trace(n, densityMatrix_full, twoelMatrix_full)*0.5;
      
  if(do_xc)
    {
      ergo_real dftEnergy;
      // construct matrix Fxc and xc_energy
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "calling dft_get_xc_mt, n = %i, noOfElectrons = %i",
		n, noOfElectrons);

      output_current_memory_usage(LOG_AREA_SCF, "before dft_get_xc_mt");

      // call dft_get_xc_mt. The resulting F_xc matrix is added to the existing twoelMatrix.


#ifndef SKIP_DFT_FLAG
      dft_get_xc_mt(noOfElectrons, densityMatrix_full,
                    &basisInfo, &molecule, gridParams,
		    twoelMatrix_full, &dftEnergy);
#endif


      output_current_memory_usage(LOG_AREA_SCF, "after  dft_get_xc_mt");
      
      *energy_2el += dftEnergy;
    }

  {
    std::vector<ergo_real> fullTmp(twoelMatrix_full, twoelMatrix_full + n*n);
    twoelMatrix_sparse.assignFromFull(fullTmp,
				      permutationHML,
				      permutationHML);
  }
  delete [] twoelMatrix_full;
  delete [] densityMatrix_full;

  return 0;
}



int
get_2e_matrices_and_energy_sparse_unrestricted(const BasisInfoStruct & basisInfo, 
					       const Molecule& molecule,
					       const IntegralInfo& integralInfo, 
					       const JK::ExchWeights & CAM_params,
					       symmMatrix & twoelMatrix_sparse_alpha, 
					       symmMatrix & twoelMatrix_sparse_beta, 
					       symmMatrix & densityMatrix_sparse_alpha,
					       symmMatrix & densityMatrix_sparse_beta,
					       const JK::Params& J_K_params,
					       const Dft::GridParams& gridParams,
					       int do_xc,
					       ergo_real* energy_2el,
					       int noOfElectrons,
					       mat::SizesAndBlocks const & matrix_size_block_info,
					       std::vector<int> const & permutationHML,
					       std::vector<int> const & inversePermutationHML)
{
  if(J_K_params.use_naive_fockmat_constr == 1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error: use_naive_fockmat_constr not implemented for unrestricted case.");
      return -1;
    }

  if(do_xc == 0 && J_K_params.use_fmm == 1)
    {
      // Simple Hartree-Fock case, with FMM. 
      // In this case we use a special function to get away with less memory usage.
      if(get_2e_matrices_and_energy_simple_HF_sparse_unrestricted(basisInfo,
								  molecule,
								  integralInfo, 
								  twoelMatrix_sparse_alpha, 
								  twoelMatrix_sparse_beta, 
								  densityMatrix_sparse_alpha, 
								  densityMatrix_sparse_beta, 
								  J_K_params,
								  CAM_params,
								  matrix_size_block_info,
								  permutationHML,
								  inversePermutationHML,
								  energy_2el) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_2e_matrices_and_energy_simple_HF_sparse_unrestricted");
	  return -1;
	}
      return 0;
    }
  if(J_K_params.use_fmm == 1 || CAM_params.computeRangeSeparatedExchange)
    {
      // Simple case, with FMM. 
      if(get_2e_matrices_and_energy_simple_sparse_unrestricted(basisInfo,
						               molecule,
						               integralInfo, 
							       CAM_params,
						twoelMatrix_sparse_alpha,
                                                twoelMatrix_sparse_beta, 
						densityMatrix_sparse_alpha,
                                                densityMatrix_sparse_beta, 
						J_K_params,
					        gridParams,
						matrix_size_block_info,
						permutationHML,
						inversePermutationHML,
						energy_2el,
						do_xc,
						noOfElectrons) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in get_2e_matrices_and_energy_simple_sparse_unrestricted");
	  return -1;
	}
      return 0;
    }

  int n = basisInfo.noOfBasisFuncs;
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "get_2e_matrices_and_energy_sparse_unrestricted allocating full n*n matrices.");
  std::vector<ergo_real> twoelMatrix_full_alpha(n*n);
  std::vector<ergo_real> twoelMatrix_full_beta(n*n);
  std::vector<ergo_real> densityMatrix_full_alpha(n*n);
  std::vector<ergo_real> densityMatrix_full_beta(n*n);
  {
    std::vector<ergo_real> fullTmp(n*n);
    densityMatrix_sparse_alpha.fullMatrix(fullTmp,
					  inversePermutationHML,
					  inversePermutationHML);
    std::copy (fullTmp.begin(), fullTmp.end(), densityMatrix_full_alpha.begin());
    densityMatrix_sparse_beta.fullMatrix(fullTmp,
					 inversePermutationHML,
					 inversePermutationHML);
    std::copy (fullTmp.begin(), fullTmp.end(), densityMatrix_full_beta.begin());
  }

  // get both J and K at the same time
  if(J_K_params.use_differential_density == 1)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error: use_differential_density not impl for unrestricted case.");
      return -1;
    }
  else
    {
      if(CAM_params.computeRangeSeparatedExchange)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, 
		    "error in get_2e_matrices_and_energy_simple_sparse_unrestricted: "
		    "cannot use compute_JK_single_box when CAM params are used.");
	  return -1;
	}

      // Get both J and K at once
      ergo_real* J_list[2];
      ergo_real* K_list[2];
      J_list[0] = new ergo_real[n*n];
      J_list[1] = new ergo_real[n*n];
      K_list[0] = new ergo_real[n*n];
      K_list[1] = new ergo_real[n*n];
      ergo_real* D_list[2];
      D_list[0] = &densityMatrix_full_alpha[0];
      D_list[1] = &densityMatrix_full_beta[0];

      // Use smallest of threshold_J and threshold_K.
      ergo_real threshold_JK = J_K_params.threshold_J;
      if(J_K_params.threshold_K < threshold_JK)
	threshold_JK = J_K_params.threshold_K;

      if(compute_JK_single_box(basisInfo,
			       integralInfo,
			       J_list[0],
			       K_list[0],
			       D_list[0],
			       threshold_JK) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_JK_single_box");
	  return -1;
	}
      if(compute_JK_single_box(basisInfo,
			       integralInfo,
			       J_list[1],
			       K_list[1],
			       D_list[1],
			       threshold_JK) != 0)
	{
	  do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_JK_single_box");
	  return -1;
	}

      ergo_real hf_weight = CAM_params.alpha;
      for(int i = 0; i < n*n; i++)
	{
	  twoelMatrix_full_alpha[i] = J_list[0][i] + J_list[1][i] + 2 * K_list[0][i] * hf_weight;
	  twoelMatrix_full_beta [i] = J_list[0][i] + J_list[1][i] + 2 * K_list[1][i] * hf_weight;
	}
      delete [] J_list[0];
      delete [] J_list[1];
      delete [] K_list[0];
      delete [] K_list[1];
    }


  *energy_2el = 0;
  *energy_2el += 0.5 * get_trace(n, &densityMatrix_full_alpha[0],
                                 &twoelMatrix_full_alpha[0]);
  *energy_2el += 0.5 * get_trace(n, &densityMatrix_full_beta[0],
                                 &twoelMatrix_full_beta[0] );
     
 
  if(do_xc)
    {
      ergo_real dftEnergy = 0;

      std::vector<ergo_real> xcMatrix_full_alpha(n*n);
      std::vector<ergo_real> xcMatrix_full_beta(n*n);
      memset(&xcMatrix_full_alpha[0], 0, n*n*sizeof(ergo_real));
      memset(&xcMatrix_full_beta[0] , 0, n*n*sizeof(ergo_real));

      /* construct xc matrices and xc_energy. */
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "calling dft_get_uxc, n = %i", n);

#ifndef SKIP_DFT_FLAG
      dft_get_uxc_mt(noOfElectrons, 
		     &densityMatrix_full_alpha[0], &densityMatrix_full_beta[0],
                     &basisInfo, &molecule, gridParams,
		     &xcMatrix_full_alpha[0], &xcMatrix_full_beta[0],
                     &dftEnergy);
#endif

      *energy_2el += dftEnergy;
      
      output_current_memory_usage(LOG_AREA_SCF, "after dft_get_uxc_mt");

      for(int i = 0; i < n*n; i++)
        {
          twoelMatrix_full_alpha[i] += xcMatrix_full_alpha[i];
          twoelMatrix_full_beta [i] += xcMatrix_full_beta [i];
        }
    }
  twoelMatrix_sparse_alpha.assignFromFull(twoelMatrix_full_alpha,
                                          permutationHML,
                                          permutationHML);
  twoelMatrix_sparse_beta.assignFromFull(twoelMatrix_full_beta,
                                         permutationHML,
                                         permutationHML);

  return 0;
}

/** Computes  G_c and G_o.
    G_c is defined as J_a + J_b + K_a + K_b.
    G_o is defined as J_a + J_b + K_a.
*/
int
get_2e_matrices_and_energy_restricted_open(const BasisInfoStruct & basisInfo, 
					   const Molecule& molecule,
					   const IntegralInfo& integralInfo, 
					   const JK::ExchWeights & CAM_params,
					   symmMatrix & twoelMatrix_Fc,
					   symmMatrix & twoelMatrix_Fo,
					   symmMatrix & densityMatrix_sparse_alpha,
					   symmMatrix & densityMatrix_sparse_beta,
					   const JK::Params& J_K_params,
					   const Dft::GridParams& gridParams,
					   int do_xc,
					   ergo_real* energy_2el,
					   int noOfElectrons,
					   mat::SizesAndBlocks const & matrix_size_block_info,
					   std::vector<int> const & permutationHML,
					   std::vector<int> const & inversePermutationHML)
{
  if(J_K_params.use_naive_fockmat_constr)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF,
		"error: use_naive_fockmat_constr not implemented for unrestricted case.");
      return -1;
    }

  symmMatrix J_tot, K_alpha, K_beta;
  J_tot.resetSizesAndBlocks(matrix_size_block_info,
				   matrix_size_block_info);
  K_alpha.resetSizesAndBlocks(matrix_size_block_info,
				     matrix_size_block_info);
  K_beta.resetSizesAndBlocks(matrix_size_block_info,
				    matrix_size_block_info);

  // Compute total density matrix
  symmMatrix D_tot(densityMatrix_sparse_alpha);
  D_tot += densityMatrix_sparse_beta;

  // FMM for J
  if(compute_J_by_boxes_sparse(basisInfo,
			       integralInfo,
			       J_K_params,
			       J_tot,
			       D_tot,
			       permutationHML) != 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_J_by_boxes_sparse");
      return -1;
    }
  D_tot.clear();

  if(CAM_params.alpha != ergo_real(0) ||
     CAM_params.computeRangeSeparatedExchange ) {
    // exchange matrices
    if(compute_K_by_boxes_sparse(basisInfo,
                                 integralInfo,
                                 CAM_params,
                                 J_K_params,
                                 K_alpha,
                                 densityMatrix_sparse_alpha,
				 permutationHML,
				 inversePermutationHML) != 0)
      {
        do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes_sparse");
        return -1;
      }

    if(compute_K_by_boxes_sparse(basisInfo,
                                 integralInfo,
                                 CAM_params,
                                 J_K_params,
                                 K_beta,
                                 densityMatrix_sparse_beta,
				 permutationHML,
				 inversePermutationHML) != 0)
      {
        do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error in compute_K_by_boxes_sparse");
        return -1;
      }
    ergo_real hf_weight = CAM_params.computeRangeSeparatedExchange
      ? 1.0 : CAM_params.alpha;
    K_alpha *= hf_weight*2.0;
    K_beta  *= hf_weight*2.0;
  }

  *energy_2el = 0;
  *energy_2el +=  0.5 * symmMatrix::trace_ab(densityMatrix_sparse_alpha, J_tot);
  *energy_2el +=  0.5 * symmMatrix::trace_ab(densityMatrix_sparse_alpha, K_alpha);
  *energy_2el +=  0.5 * symmMatrix::trace_ab(densityMatrix_sparse_beta, J_tot);
  *energy_2el +=  0.5 * symmMatrix::trace_ab(densityMatrix_sparse_beta, K_beta);
  twoelMatrix_Fc  = J_tot;
  twoelMatrix_Fc += ergo_real(0.5)*K_alpha;
  twoelMatrix_Fc += ergo_real(0.5)*K_beta;
  twoelMatrix_Fo  = J_tot;
  twoelMatrix_Fo += K_alpha;

  J_tot.clear(); K_alpha.clear(); K_beta.clear();

  if(do_xc)
    {
      int n = basisInfo.noOfBasisFuncs;
      ergo_real dftEnergy = 0;

      std::vector<ergo_real> density_full_alpha(n*n);
      std::vector<ergo_real> density_full_beta(n*n);
      densityMatrix_sparse_alpha.fullMatrix(density_full_alpha,
					    inversePermutationHML,
					    inversePermutationHML);
      densityMatrix_sparse_beta.fullMatrix(density_full_beta,
					    inversePermutationHML,
					    inversePermutationHML);

      std::vector<ergo_real> xcMatrix_full_alpha(n*n,0);
      std::vector<ergo_real> xcMatrix_full_beta (n*n,0);

      /* construct xc matrices and xc_energy. */
      do_output(LOG_CAT_INFO, LOG_AREA_SCF, "calling dft_get_uxc, n = %i", n);

#ifndef SKIP_DFT_FLAG
      dft_get_uxc_mt(noOfElectrons, 
		     &density_full_alpha[0], &density_full_beta[0],
                     &basisInfo, &molecule, gridParams,
		     &xcMatrix_full_alpha[0], &xcMatrix_full_beta[0], &dftEnergy);
#endif

      *energy_2el += dftEnergy;
      
      output_current_memory_usage(LOG_AREA_SCF, "after dft_get_uxc_mt");
      symmMatrix xcSparse_alpha, xcSparse_beta;
      xcSparse_alpha.resetSizesAndBlocks(matrix_size_block_info,
						matrix_size_block_info);
      xcSparse_beta.resetSizesAndBlocks(matrix_size_block_info,
					       matrix_size_block_info);

      xcSparse_alpha.assignFromFull(xcMatrix_full_alpha,
				    permutationHML, permutationHML);
      xcSparse_beta.assignFromFull(xcMatrix_full_beta,
				   permutationHML, permutationHML);

      twoelMatrix_Fc += ergo_real(0.5)*xcSparse_alpha;
      twoelMatrix_Fc += ergo_real(0.5)*xcSparse_beta;

      twoelMatrix_Fo += xcSparse_alpha;
    }



  return 0;
}







int
compute_FDSminusSDF_sparse(int n, 
			   symmMatrix & F_symm, 
			   symmMatrix & D_symm, 
			   symmMatrix & S_symm, 
			   normalMatrix & result, 
			   ergo_real sparse_threshold)
{
  // FIXME: we should be able to do this without so many temporary objects.
  // We should probably use something like FDS-SDF = FDS - (FDS)^T
  // FIXME: Implement eucl_thresh() for general matrices.

  Util::TimeMeter timeMeter;

  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, "compute_FDSminusSDF_sparse() start.");
  Util::TimeMeter timeMeterWriteAndReadAll;
  std::string sizesStr = mat::FileWritable::writeAndReadAll();
  timeMeterWriteAndReadAll.print(LOG_AREA_DENSFROMF, "FileWritable::writeAndReadAll");
  do_output(LOG_CAT_INFO, LOG_AREA_DENSFROMF, ((std::string)"writeAndReadAll: '" + sizesStr).c_str());

  normalMatrix S(S_symm);
  normalMatrix FD;
  normalMatrix DF;

  S_symm.writeToFile();
  S.writeToFile();


  output_current_memory_usage(LOG_AREA_SCF, "compute_FDSminusSDF_sparse before computing FD (symm-symm multiply)");

  Util::TimeMeter timeMeterFD;

  mat::SingletonForTimings::instance().reset();

  FD  = (ergo_real)1.0 * F_symm * D_symm;

#ifdef MAT_USE_ALLOC_TIMER
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "compute_FDSminusSDF_sparse, F*D mult mat alloc time: %15.5f seconds", 
	    mat::SingletonForTimings::instance().getAccumulatedTime());
#endif

  timeMeterFD.print(LOG_AREA_SCF, "compute_FDSminusSDF_sparse, F*D mult");

  output_current_memory_usage(LOG_AREA_SCF, 
			      "compute_FDSminusSDF_sparse after computing FD");

  output_sparsity(n, FD, "FD before eucl truncation");
  FD.eucl_thresh(sparse_threshold);
  output_sparsity(n, FD, "FD after  eucl truncation");

  F_symm.writeToFile();
  D_symm.writeToFile();

  S.readFromFile();

  result = (ergo_real)1.0 * FD * S; // result now contains FDS

  output_current_memory_usage(LOG_AREA_SCF, 
			      "compute_FDSminusSDF_sparse "
			      "after computing FDS");

  output_sparsity(n, result, "FDS before eucl truncation");
  result.eucl_thresh(sparse_threshold);
  output_sparsity(n, result, "FDS after  eucl truncation");

  // We no longer need FD, set it to zero and free memory
  FD.clear();

  S.writeToFile();
  result.writeToFile();

  // Read D before F to reduce problems with memory fragmentation;
  // we assume D is more dense than F.
  D_symm.readFromFile();
  F_symm.readFromFile();

  output_current_memory_usage(LOG_AREA_SCF, 
			      "compute_FDSminusSDF_sparse before "
			      "computing DF (symm-symm multiply)");
  
  DF  = (ergo_real)1.0 * D_symm * F_symm;
  
  output_current_memory_usage(LOG_AREA_SCF, 
			      "compute_FDSminusSDF_sparse after computing DF");
  
  DF.eucl_thresh(sparse_threshold);
  
  F_symm.writeToFile();
  D_symm.writeToFile();
  
  result.readFromFile();
  S.readFromFile();
  
  result = (ergo_real)1.0 * S * DF + ((ergo_real)(-1.0) * result);
  
  output_current_memory_usage(LOG_AREA_SCF, 
			      "compute_FDSminusSDF_sparse after "
			      "computing result");
  
  S.clear();
  FD.clear();
  DF.clear();
  
  result.eucl_thresh(sparse_threshold);

  // read input matrices back from file
  D_symm.readFromFile();
  F_symm.readFromFile();
  S_symm.readFromFile();

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "compute_FDSminusSDF_sparse ending OK.");
  timeMeter.print(LOG_AREA_SCF, "compute_FDSminusSDF_sparse");
  
  return 0;
}


int 
determine_number_of_electrons_unrestricted(int noOfElectrons, 
					   int alpha_beta_diff, 
					   int* noOfElectrons_alpha, 
					   int* noOfElectrons_beta)
{
  *noOfElectrons_alpha = (noOfElectrons + alpha_beta_diff) / 2;
  *noOfElectrons_beta  = noOfElectrons - *noOfElectrons_alpha;
  int sum = *noOfElectrons_alpha + *noOfElectrons_beta;
  if(sum != noOfElectrons)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "error: n_alpha + n_beta != N");
      return -1;
    }
  int diff = *noOfElectrons_alpha - *noOfElectrons_beta;
  if(diff != alpha_beta_diff)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, 
		"error: n_alpha - n_beta != alpha_beta_diff");
      return -1;
    }
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, 
	    "noOfElectrons_alpha  noOfElectrons_beta  =  %i %i", 
	    *noOfElectrons_alpha, *noOfElectrons_beta);
  return 0;
}


static void
get_mulliken_charges(const symmMatrix & densityMatrix,
		     const symmMatrix & S_symm,
		     const BasisInfoStruct & basisInfo,
		     mat::SizesAndBlocks const & matrix_size_block_info,
		     std::vector<int> const & permutationHML,
		     std::vector<int> const & inversePermutationHML,
		     const Molecule& molecule,
		     std::vector<ergo_real> & resultVector)
{
  int nAtoms = molecule.getNoOfAtoms();
  resultVector.resize(nAtoms);
  for(int i = 0; i < nAtoms; i++)
    resultVector[i] = 0;
  // Create mapping between basis functions and atoms.
  int n = basisInfo.noOfBasisFuncs;
  std::vector<int> atomMapping(n);
  for(int i = 0; i < n; i++) {
    // Check if this basis function is near any atom.
    int foundAtomIndex = -1;
    for(int k = 0; k < nAtoms; k++) {
      ergo_real dx = molecule.getAtom(k).coords[0] - basisInfo.basisFuncList[i].centerCoords[0];
      ergo_real dy = molecule.getAtom(k).coords[1] - basisInfo.basisFuncList[i].centerCoords[1];
      ergo_real dz = molecule.getAtom(k).coords[2] - basisInfo.basisFuncList[i].centerCoords[2];
      ergo_real distance = template_blas_sqrt(dx*dx + dy*dy + dz*dz);
      if(distance < 1e-4)
	foundAtomIndex = k;
    }
    if(foundAtomIndex < 0)
      throw "Error in get_mulliken_charges: all basis functions are not centered on atoms.";
    atomMapping[i] = foundAtomIndex;
  }

  normalMatrix D(densityMatrix);
  normalMatrix S(S_symm);
  // Now get all elements of the matrix D.
  std::vector<int> rowind;
  std::vector<int> colind;
  std::vector<ergo_real> values_D;
  D.get_all_values(rowind,
		   colind,
		   values_D,
		   inversePermutationHML,
		   inversePermutationHML);
  int nvalues = values_D.size();

  // Now get corresponding elements of S.
  std::vector<ergo_real> values_S(nvalues);
  for(int i = 0; i < nvalues; i++)
    values_S[i] = 0;
  S.get_values(rowind,
	       colind,
	       values_S,
	       permutationHML, //inversePermutationHML,
	       permutationHML);//inversePermutationHML);

  // Now go through all elements of D and S and make the corresponding contributions to the result vector.
  for(int i = 0; i < nvalues; i++) {
    int row = rowind[i];
    int col = colind[i];
    int atomIndex1 = atomMapping[row];
    int atomIndex2 = atomMapping[col];
    resultVector[atomIndex1] += 0.5 * values_D[i] * values_S[i];
    resultVector[atomIndex2] += 0.5 * values_D[i] * values_S[i];
  }
}

void
do_mulliken_atomic_charges(const symmMatrix & densityMatrix,
			   const symmMatrix & S_symm,
			   const BasisInfoStruct & basisInfo,
			   mat::SizesAndBlocks const & matrix_size_block_info,
			   std::vector<int> const & permutationHML,
			   std::vector<int> const & inversePermutationHML,
			   const Molecule& molecule)
{
  // Compute Mulliken atomic charges.
  // First get vector with summed electron charge for all atoms.
  std::vector<ergo_real> electronChargeSums;
  get_mulliken_charges(densityMatrix, S_symm, basisInfo, matrix_size_block_info, permutationHML, inversePermutationHML, molecule, electronChargeSums);
  int nAtoms = molecule.getNoOfAtoms();
  ergo_real chargeSum = 0;
  for(int i = 0; i < nAtoms; i++)
    chargeSum += electronChargeSums[i];
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Sum of Mulliken charges (electrons): %15.7f", (double)chargeSum);
  std::vector<ergo_real> atomicCharges(nAtoms);
  for(int i = 0; i < nAtoms; i++)
    atomicCharges[i] = molecule.getAtom(i).charge - electronChargeSums[i];
  for(int i = 0; i < nAtoms; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Mulliken charge of atom %d = %15.7f", i, (double)atomicCharges[i]);    
  chargeSum = 0;
  for(int i = 0; i < nAtoms; i++)
    chargeSum += atomicCharges[i];
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Sum of Mulliken atomic charges: %15.7f", (double)chargeSum);
}

void
do_mulliken_spin_densities(const symmMatrix & spinDensityMatrix,
			   const symmMatrix & S_symm,
			   const BasisInfoStruct & basisInfo,
			   mat::SizesAndBlocks const & matrix_size_block_info,
			   std::vector<int> const & permutationHML,
			   std::vector<int> const & inversePermutationHML,
			   const Molecule& molecule)
{
  std::vector<ergo_real> mullikenCharges;
  get_mulliken_charges(spinDensityMatrix, S_symm, basisInfo, matrix_size_block_info, permutationHML, inversePermutationHML, molecule, mullikenCharges);
  int nAtoms = molecule.getNoOfAtoms();
  for(int i = 0; i < nAtoms; i++)
    do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Mulliken spin density of atom %d = %15.7f", i, (double)mullikenCharges[i]);
  ergo_real chargeSum = 0;
  for(int i = 0; i < nAtoms; i++)
    chargeSum += mullikenCharges[i];
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "Sum of Mulliken atomic spin densities: %15.7f", (double)chargeSum);  
}

static int write_gcube_file_header(FILE* f, 
				   const char* firstLine,
				   const Molecule& m, 
				   double originX, double originY, double originZ, 
				   int NX, double incrX,
				   int NY, double incrY,
				   int NZ, double incrZ) {
  fprintf(f, "%s\n\n", firstLine);
  fprintf(f, "%d  %9.5f  %9.5f  %9.5f\n", m.getNoOfAtoms(), originX, originY, originZ);
  fprintf(f, "%d  %9.5f    0.0    0.0\n", NX, incrX);
  fprintf(f, "%d    0.0  %9.5f    0.0\n", NY, incrY);
  fprintf(f, "%d    0.0    0.0  %9.5f\n", NZ, incrZ);
  for(int i = 0; i < m.getNoOfAtoms(); i++)
    fprintf(f, "%d 0.0 %9.5f %9.5f %9.5f\n", (int)m.getAtom(i).charge, (double)m.getAtom(i).coords[0], (double)m.getAtom(i).coords[1], (double)m.getAtom(i).coords[2]);
  return 0;
}




void get_exp_value_pos_operator(const BasisInfoStruct & basisInfo,
				const Molecule& molecule,
				const symmMatrix & densityMatrix,
				mat::SizesAndBlocks const & matrix_size_block_info,
				std::vector<int> const & permutationHML,
				std::vector<ergo_real> &mean,
				std::vector<ergo_real> &std)

{
  ergo_real expval; 
  int xi = 0, yi = 0, zi = 0;
  mean.clear(); mean.resize(3);
  std.clear(); std.resize(3);

  symmMatrix opMatrix;
  opMatrix.resetSizesAndBlocks(matrix_size_block_info, matrix_size_block_info);

  for(int coord = 0; coord < 3; ++coord)
    {
      switch(coord) 
	{
	case 0: xi = 1; break;
	case 1: yi = 1; break;
	case 2: zi = 1; break;
	default: throw "Error in get_exp_value_pos_operator: wrong value of coord";
	}
      // compute expectated value of position operator
      if(compute_operator_matrix_sparse_symm(basisInfo, xi, yi, zi, opMatrix, permutationHML) != 0)
	throw "Error in compute_operator_matrix_sparse_symm";
      expval = symmMatrix::trace_ab(densityMatrix, opMatrix); // E[X]
      mean[coord] = expval;

      // compute standard deviation
      if(compute_operator_matrix_sparse_symm(basisInfo, 2*xi, 2*yi, 2*zi, opMatrix, permutationHML) != 0)
	throw "Error in compute_operator_matrix_sparse_symm";
      expval = symmMatrix::trace_ab(densityMatrix, opMatrix); // E[X^2]
      std[coord] = template_blas_sqrt(expval - mean[coord]*mean[coord]); // sqrt[ E[X^2] - (E[X])^2 ]

      xi = 0; yi = 0; zi = 0;
    }
}




void 
do_density_images(const BasisInfoStruct & basisInfo,
		  const Molecule& molecule,
		  const ergo_real* densityMatrixFull_tot, 
		  const ergo_real* densityMatrixFull_spin,
		  double output_density_images_boxwidth,
		  const std::string &filename_id)
{

  int nPrims1 = 0;
  int nPrims2 = 0;
  
  DistributionSpecStruct* rho1 = NULL;
  DistributionSpecStruct* rho2 = NULL;
  ergo_real cutoff = 1e-7;//1e-7;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_density_images, calling get_density() with cutoff = %8.4g.", cutoff);

  nPrims1 = get_no_of_primitives_for_density(cutoff, densityMatrixFull_tot , basisInfo);
  nPrims2 = get_no_of_primitives_for_density(cutoff, densityMatrixFull_spin, basisInfo);
  rho1 = new DistributionSpecStruct[nPrims1];
  rho2 = new DistributionSpecStruct[nPrims2];

  Util::TimeMeter timeMeterGetDens;
  int n1 = get_density(basisInfo,
		       densityMatrixFull_tot,
		       cutoff,
		       nPrims1,
		       rho1);      
  int n2 = get_density(basisInfo,
		       densityMatrixFull_spin,
		       cutoff,
		       nPrims2,
		       rho2);
  if(n1 < 0 || n2 < 0)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "aborting do_density_images(): (n1 < 0 || n2 < 0)");
      throw "aborting do_density_images(): (n1 < 0 || n2 < 0)";
    }
  nPrims1 = n1;
  nPrims2 = n2;
  timeMeterGetDens.print(LOG_AREA_SCF, "get_density() calls");

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_density_images, after get_density() calls: nPrims1 = %d, nPrims2 = %d.", nPrims1, nPrims2);

  // get min max coords for molecule
  ergo_real minlist[3];
  ergo_real maxlist[3];
  int coordno;
  for(coordno = 0; coordno < 3; coordno++)
    {
      minlist[coordno] = template_blas_get_num_limit_max<ergo_real>();
      maxlist[coordno] = template_blas_get_num_limit_min<ergo_real>();
    }
  for(int i = 0; i < molecule.getNoOfAtoms(); i++)
    {
      for(coordno = 0; coordno < 3; coordno++)
	{
	  ergo_real curr = molecule.getAtom(i).coords[coordno];
	  if(curr < minlist[coordno])
	    minlist[coordno] = curr;
	  if(curr > maxlist[coordno])
	    maxlist[coordno] = curr;
	}
    }
  // OK, now we have min max coords for molecule.
  // Add some margin.
  ergo_real margin = 4.4;
  for(coordno = 0; coordno < 3; coordno++)
    {
      minlist[coordno] -= margin;
      maxlist[coordno] += margin;
    }

  ergo_real boxWidthApprox = output_density_images_boxwidth; //0.3;
  int Nx = (int)((maxlist[0] - minlist[0]) / boxWidthApprox);
  int Ny = (int)((maxlist[1] - minlist[1]) / boxWidthApprox);
  int Nz = (int)((maxlist[2] - minlist[2]) / boxWidthApprox);

  // Make sure Nx Ny Nz are odd numbers (this seems to be required for the electrostatic potential computation feature in Gabedit).
  if(Nx % 2 == 0) Nx++;
  if(Ny % 2 == 0) Ny++;
  if(Nz % 2 == 0) Nz++;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_density_images: Nx Ny Nz = %5i %5i %5i", Nx, Ny, Nz);
  
  ergo_real boxWidth_x = (maxlist[0] - minlist[0]) / Nx;
  ergo_real boxWidth_y = (maxlist[1] - minlist[1]) / Ny;
  ergo_real boxWidth_z = (maxlist[2] - minlist[2]) / Nz;

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_density_images, preparing to write Gabedit gcube files.");

  // Create Gabedit "gcube" density files.
  std::ostringstream name_f1, name_f2;
  name_f1 << "density" << filename_id << ".gcube";
  name_f2 << "spindensity" << filename_id << ".gcube";

  FILE* f1 = fopen(name_f1.str().c_str(), "wt");
  FILE* f2 = fopen(name_f2.str().c_str(), "wt");
  if(f1 == NULL || f2 == NULL)
    {
      do_output(LOG_CAT_ERROR, LOG_AREA_SCF, "aborting do_density_images: (f1 == NULL || f2 == NULL)");
      throw "aborting do_density_images: (f1 == NULL || f2 == NULL)";
    }
  double originX = minlist[0];
  double originY = minlist[1];
  double originZ = minlist[2];
  double dx = (maxlist[0] - minlist[0]) / Nx;
  double dy = (maxlist[1] - minlist[1]) / Ny;
  double dz = (maxlist[2] - minlist[2]) / Nz;
  // Shift origin coordinates by half of the box width in each direction (it seems Gabedit is interpreting the values in that way).
  originX += 0.5 * dx;
  originY += 0.5 * dy;
  originZ += 0.5 * dz;
  write_gcube_file_header(f1, 
			  "Density file in gcube format, generated by the Ergo program.",
			  molecule, 
			  originX, originY, originZ, 
			  Nx, dx,
			  Ny, dy,
			  Nz, dz);
  write_gcube_file_header(f2, 
			  "Spin density file in gcube format, generated by the Ergo program.",
			  molecule, 
			  originX, originY, originZ, 
			  Nx, dx,
			  Ny, dy,
			  Nz, dz);
  // Now compute density values in all needed boxes.
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_density_images, computing density values needed for Gabedit gcube files.");
  Util::TimeMeter timeMeterGetDensValues;
  std::vector< std::vector< std::vector<double> > > vector_dens(Nx);
  std::vector< std::vector< std::vector<double> > > vector_spin(Nx);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ix = 0; ix < Nx; ix++) {
    vector_dens[ix].resize(Ny);
    vector_spin[ix].resize(Ny);
    for(int iy = 0; iy < Ny; iy++) {
      vector_dens[ix][iy].resize(Nz);
      vector_spin[ix][iy].resize(Nz);
      for(int iz = 0; iz < Nz; iz++) {
	ergo_real minVect[3];
	ergo_real maxVect[3];
	minVect[0] = minlist[0] + ix * (maxlist[0] - minlist[0]) / Nx;
	maxVect[0] = minVect[0] + boxWidth_x;
	minVect[1] = minlist[1] + iy * (maxlist[1] - minlist[1]) / Ny;
	maxVect[1] = minVect[1] + boxWidth_y;
	minVect[2] = minlist[2] + iz * (maxlist[2] - minlist[2]) / Nz;
	maxVect[2] = minVect[2] + boxWidth_z;
	ergo_real dens = integrate_density_in_box_2(nPrims1,
						    rho1,
						    minVect, 
						    maxVect);
	ergo_real spin = integrate_density_in_box_2(nPrims2,
						    rho2,
						    minVect, 
						    maxVect);
	vector_dens[ix][iy][iz] = dens;
	vector_spin[ix][iy][iz] = spin;
      } // END FOR iz
    } // END FOR iy
  } // END FOR ix
  timeMeterGetDensValues.print(LOG_AREA_SCF, "getting density values");
  // Now write to files.
  Util::TimeMeter timeMeterWriteFiles;
  ergo_real sum_dens = 0;
  ergo_real sum_spin = 0;
  for(int ix = 0; ix < Nx; ix++)
    for(int iy = 0; iy < Ny; iy++) {
      int count = 0;
      for(int iz = 0; iz < Nz; iz++) {
	ergo_real dens = vector_dens[ix][iy][iz];
	ergo_real spin = vector_spin[ix][iy][iz];
	sum_dens += dens;
	sum_spin += spin;
	ergo_real volume = dx*dy*dz;
	fprintf(f1, " %9.5f", (double)(dens / volume));
	fprintf(f2, " %9.5f", (double)(spin / volume));
	count++;
	if(count % 6 == 0) {
	  fprintf(f1, "\n");
	  fprintf(f2, "\n");
	}
      }
      if(count % 6 != 0) {
	fprintf(f1, "\n");
	fprintf(f2, "\n");
      }
    }
  timeMeterWriteFiles.print(LOG_AREA_SCF, "writing to gcube files");
  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "sum_dens = %9.5f, sum_spin = %9.5f", sum_dens, sum_spin);

  fclose(f1);
  fclose(f2);

  do_output(LOG_CAT_INFO, LOG_AREA_SCF, "do_density_images done, Gabedit gcube files written.");
}


void
get_hf_weight_and_cam_params(int use_dft, 
			     ergo_real* exch_param_alpha, 
			     ergo_real* exch_param_beta, 
			     ergo_real* exch_param_mu)
{
  if(use_dft)
    {
#ifndef SKIP_DFT_FLAG
      int use_cam_params = 
        fun_get_cam_param(exch_param_alpha, exch_param_beta, exch_param_mu);
      if(!use_cam_params) {
        *exch_param_alpha = fun_get_hf_weight();
	*exch_param_beta = 0.0;
      }
#else
      throw "ERROR in get_hf_weight_and_cam_params: SKIP_DFT_FLAG given and use_dft set.";
#endif
    }
  else
    {
      *exch_param_alpha = 1.0;
      *exch_param_beta = 0;
      *exch_param_mu = 0;
    }
}


void
create_mtx_files_with_different_orderings(const symmMatrix & A,
					  const std::string & calculation_identifier,
					  const std::string & method_and_basis_set,
					  const std::vector<int> & inversePermutationHML,
					  const BasisInfoStruct & basisInfo) {
  // Original ordering
  std::stringstream ss_id_1;
  ss_id_1 << calculation_identifier << " - matrix (original ordering)";
  write_matrix_in_matrix_market_format( A, inversePermutationHML, "S_matrix_original", ss_id_1.str(), method_and_basis_set );
  // HML ordering
  int n = basisInfo.noOfBasisFuncs;
  std::vector<int> permutationNone(n);
  for(int i = 0; i < n; i++)
    permutationNone[i] = i;
  std::stringstream ss_id_2;
  ss_id_2 << calculation_identifier << " - matrix (HML ordering)";
  write_matrix_in_matrix_market_format( A, permutationNone, "S_matrix_HML", ss_id_2.str(), method_and_basis_set );
  // Factor2 ordering
  std::vector<int> permutationFactor2, inversePermutationFactor2;
  std::vector<ergo_real> xcoords(n);
  std::vector<ergo_real> ycoords(n);
  std::vector<ergo_real> zcoords(n);
  for(int i = 0; i < n; i++) {
    xcoords[i] = basisInfo.basisFuncList[i].centerCoords.v[0];
    ycoords[i] = basisInfo.basisFuncList[i].centerCoords.v[1];
    zcoords[i] = basisInfo.basisFuncList[i].centerCoords.v[2];
  }
  getMatrixPermutationOnlyFactor2(xcoords, ycoords, zcoords, 1, 1, permutationFactor2, inversePermutationFactor2);
  std::vector<int> perm(n);
  for(int i = 0; i < n; i++)
    perm[i] = permutationFactor2[ inversePermutationHML[i] ];
  std::stringstream ss_id_3;
  ss_id_3 << calculation_identifier << " - matrix (Factor2 ordering)";
  write_matrix_in_matrix_market_format( A, perm, "S_matrix_Factor2", ss_id_3.str(), method_and_basis_set );
  // Random ordering
  std::vector<int> permRandom(n);
  std::vector<bool> taken(n);
  for(int i = 0; i < n; i++)
    taken[i] = false;
  for(int i = 0; i < n; i++) {
    int randNumber = rand() % (n-i);
    int idx = randNumber;
    for(int k = 0; k < n; k++) {
      if(taken[k] && k <= idx)
	idx++;
    }
    assert(idx >= 0);
    assert(idx < n);
    assert(taken[idx] == false);
    assert(idx >= 0 && idx < n && taken[idx] == false);
    taken[idx] = true;
    permRandom[i] = idx;
  }
  for(int i = 0; i < n; i++)
    assert(taken[i] == true);
  std::stringstream ss_id_4;
  ss_id_4 << calculation_identifier << " - matrix (random ordering)";
  write_matrix_in_matrix_market_format( A, permRandom, "S_matrix_random", ss_id_4.str(), method_and_basis_set );
}


