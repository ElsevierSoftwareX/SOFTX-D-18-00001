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

/** @file mat_acc_extrapolate.h

    @brief Functionality for performing a scan testing different
    threshold values for some kind of matrix computation. Can be used
    for accuracy scans of e.g. Coulomb of HF exchange matrix
    evaluation.

    @author: Elias Rudberg <em>responsible</em>
*/

#ifndef ERGO_MAT_ACC_EXTRAPOLATE_HEADER
#define ERGO_MAT_ACC_EXTRAPOLATE_HEADER


#include <vector>


#include "matrix_utilities.h"



template<class Treal, class Tworker>
  class MatAccInvestigator
{
 public:
  explicit MatAccInvestigator(mat::SizesAndBlocks const & matrix_size_block_info_);
  void Scan(const Tworker & worker,
	    Treal firstParam, 
	    Treal stepFactor, 
	    int nSteps);
  void GetScanResult(Treal* threshList_,
		     Treal* errorList_frob_,
		     Treal* errorList_eucl_,
		     Treal* errorList_maxe_,
		     Treal* timeList_);
 private:
  mat::SizesAndBlocks matrix_size_block_info;
  int nScanSteps;
  Treal baseThresh;
  std::vector<Treal> threshList;
  std::vector<Treal> errorList_frob; // Frobenius norm
  std::vector<Treal> errorList_eucl; // Euclidean norm
  std::vector<Treal> errorList_maxe; // Max element norm
  std::vector<Treal> timeList;
};


template<class Treal, class Tworker>
  MatAccInvestigator<Treal, Tworker>::MatAccInvestigator(mat::SizesAndBlocks const & matrix_size_block_info_)
  : matrix_size_block_info(matrix_size_block_info_)
{}


template<class Treal, class Tworker>
  void MatAccInvestigator<Treal, Tworker>
  ::Scan(const Tworker & worker, 
	 Treal firstParam, 
	 Treal stepFactor, 
	 int nSteps)
{
  nScanSteps = nSteps;
  baseThresh = firstParam;
  threshList.resize(nSteps);
  errorList_frob.resize(nSteps);
  errorList_eucl.resize(nSteps);
  errorList_maxe.resize(nSteps);
  timeList.resize(nSteps);

  // Prepare matrix objects
  symmMatrix accurateMatrix;
  accurateMatrix.resetSizesAndBlocks(matrix_size_block_info,
					    matrix_size_block_info);
  symmMatrix otherMatrix;
  otherMatrix.resetSizesAndBlocks(matrix_size_block_info,
					 matrix_size_block_info);
  symmMatrix errorMatrix;
  errorMatrix.resetSizesAndBlocks(matrix_size_block_info,
					 matrix_size_block_info);

  // Compute "accurate" matrix
  worker.ComputeMatrix(firstParam, accurateMatrix);
  // Compute other matrices and compare them to "accurate" matrix
  Treal currParam = firstParam;
  for(int i = 0; i < nSteps; i++)
    {
      currParam *= stepFactor;
      time_t startTime, endTime;
      time(&startTime);
      worker.ComputeMatrix(currParam, otherMatrix);
      time(&endTime);
      timeList[i] = endTime - startTime;
      threshList[i] = currParam;
      // Compute error matrix
      errorMatrix = otherMatrix;
      errorMatrix += (ergo_real)(-1) * accurateMatrix;
      // Compute different norms of error matrix
      // Frobenius norm
      errorList_frob[i] = errorMatrix.frob();
      // Euclidean norm
      Treal euclAcc = 1e-11;
      errorList_eucl[i] = errorMatrix.eucl(euclAcc);
      // Max element norm
      errorList_maxe[i] = compute_maxabs_sparse(errorMatrix);
    }
  
}


template<class Treal, class Tworker>
  void MatAccInvestigator<Treal, Tworker>
  ::GetScanResult(Treal* threshList_,
		  Treal* errorList_frob_,
		  Treal* errorList_eucl_,
		  Treal* errorList_maxe_,
		  Treal* timeList_)
{
  for(int i = 0; i < nScanSteps; i++)
    {
      threshList_[i] = threshList[i];
      errorList_frob_[i] = errorList_frob[i];
      errorList_eucl_[i] = errorList_eucl[i];
      errorList_maxe_[i] = errorList_maxe[i];
      timeList_      [i] = timeList      [i];
    }
}




#endif
