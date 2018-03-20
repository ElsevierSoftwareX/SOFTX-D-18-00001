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

#include <stdio.h>
#include <unistd.h>
#include <memory>
#include <limits>
#include "matrix_typedefs.h"
#include "matrix_utilities.h"
#include "utilities.h"
#include "config.h" // Needed to get the PRECISION_SINGLE macro

/** @file mat_performance_test.cc Performs some matrix-matrix
    multiplication operations and outputs timings. The point is to
    show how different block sizes affects the performance of
    matrix-matrix multiplication using the hierarchic matrix
    library. */

static void get_random_matrix_full(int n, std::vector<ergo_real> & fullMat) {
  for (int col = 0; col < n; ++col)
    for (int row = 0; row < n; ++row) {
      ergo_real randomNumber = ((ergo_real)rand() / (ergo_real)RAND_MAX);
      fullMat[row + col * n] = randomNumber;
    }
}

static void get_matrix_from_full(normalMatrix & result,
				 mat::SizesAndBlocks sizeBlockInfo,
				 const std::vector<ergo_real> & fullMat) {
  result.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  result.assignFromFull(fullMat);
}

static void verify_mmul_result(int n, 
			       const std::vector<ergo_real> & fullMat_A, 
			       const std::vector<ergo_real> & fullMat_B, 
			       const normalMatrix & C) {
  // Choose some matrix elements to check.
  const int noOfElements = 2*n;
  std::vector<int> rowIdxList(noOfElements);
  std::vector<int> colIdxList(noOfElements);
  std::vector<ergo_real> valuesList(noOfElements);
  for(int i = 0; i < noOfElements; i++) {
    rowIdxList[i] = rand() % n;
    colIdxList[i] = rand() % n;
  }
  C.get_values(rowIdxList, colIdxList, valuesList);
  for(int i = 0; i < noOfElements; i++) {
    int row = rowIdxList[i];
    int col = colIdxList[i];
    ergo_real value = valuesList[i];
    ergo_real sum = 0;
    for(int k = 0; k < n; k++)
      sum += fullMat_A[k*n+row] * fullMat_B[col*n+k];
    ergo_real abs_error = template_blas_fabs(sum - value);
    ergo_real error_tolerance = 1e-7;
#ifdef PRECISION_SINGLE
    error_tolerance = 1e-5;
#endif
    if(abs_error > error_tolerance)
      throw std::runtime_error("ERROR: matrix-matrix multiplication gave wrong result.");
  }
}

int main(int argc, char *argv[])
{
  int n = 50;
  int blockSize_min = 12;
  int blockSize_max = 32;
  int blockSize_step = 2;

  if(argc > 1 && argc != 5) {
    std::cout << "This program expects the following arguments:" << std::endl;
    std::cout << "./mat_performance_test n blockSize_min blockSize_max blockSize_step" << std::endl;
    return -1;
  }

  if(argc == 5) {
    n = atoi(argv[1]);
    blockSize_min = atoi(argv[2]);
    blockSize_max = atoi(argv[3]);
    blockSize_step = atoi(argv[4]);
  }

  std::cout << "n = " << n << std::endl;
  std::cout << "blockSize_min = " << blockSize_min << std::endl;
  std::cout << "blockSize_max = " << blockSize_max << std::endl;
  std::cout << "blockSize_step = " << blockSize_step << std::endl;

  if(n <= 1)
    throw std::runtime_error("ERROR: (n <= 1)");

#ifdef _OPENMP
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  if ( !(env && (defThreads=atoi(env)) > 0) )
    defThreads = 1;
  mat::Params::setNProcs(defThreads);
  mat::Params::setMatrixParallelLevel(1);
  std::cout<<"OpenMP is used, number of threads set to "
           <<mat::Params::getNProcs()<<". Matrix parallel level: "
           <<mat::Params::getMatrixParallelLevel()<<"."<<std::endl;
#endif


  // Set up full matrices A and B.
  std::vector<ergo_real> fullMat_A(n*n);
  get_random_matrix_full(n, fullMat_A);
  std::vector<ergo_real> fullMat_B(n*n);
  get_random_matrix_full(n, fullMat_B);

  // Set up vector with blocksizes to test.
  std::vector<int> blockSizes(444);
  int count = 0;
  blockSizes.at(count) = blockSize_min;
  while(blockSizes.at(count) < blockSize_max) {
    count++;
    blockSizes.at(count) = blockSizes.at(count-1) + blockSize_step;
  }
  int extraBlockSizes[] = {64, 128, 256, 512, 1024, 2048};
  int noOfExtraBlockSizes = sizeof(extraBlockSizes) / sizeof(int);
  for(int i = 0; i < noOfExtraBlockSizes; i++) {
    int extraBlockSize = extraBlockSizes[i];
    if(blockSizes.at(count) < extraBlockSize && extraBlockSize < n)
      blockSizes.at(++count) = extraBlockSize;
  }
  blockSizes.at(++count) = n;
  int noOfBlockSizes = count;
  
  for(int i = 0; i < noOfBlockSizes; i++) {
    int blockSize = blockSizes[i];
    mat::SizesAndBlocks sizeBlockInfo;
    static const int sparseMatrixBlockFactor = 2;
    sizeBlockInfo =
      prepareMatrixSizesAndBlocks(n,
				  blockSize,
				  sparseMatrixBlockFactor,
				  sparseMatrixBlockFactor,
				  sparseMatrixBlockFactor);
    // Create random matrices A and B using the current block size.
    normalMatrix A, B;

    get_matrix_from_full(A, sizeBlockInfo, fullMat_A);
    get_matrix_from_full(B, sizeBlockInfo, fullMat_B);

    // Compute C = A * B
    normalMatrix C;
    Util::TimeMeter tm;
    C.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
    C = A * B;
    double secondsTaken = tm.get_wall_seconds() - tm.get_start_time_wall_seconds();
    std::cout << "C = A * B operation using blockSize "  
	      << std::setw(7) << blockSize << " took "  
	      << std::setw(12) << std::setprecision(5) << std::fixed << secondsTaken << " wall seconds." << std::endl;

    // Verify result.
    verify_mmul_result(n, fullMat_A, fullMat_B, C);
  }

  puts("mat_performance_test finished OK.");
  return 0;
}
