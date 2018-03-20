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

/** @file mat_thread_test.cc Performs some matrix operations and
    outputs timings. The point is to demonstrate the problem with the
    matrix library becoming extremely slow when OpenMP is used with
    gcc. Running this test with OMP_NUM_THREADS=1 and
    OMP_NUM_THREADS=2 gives dramatically different performance
    (threaded runs are about 700 times slower than serial).  */

static void get_random_positive_definite_matrix(symmMatrix & result,
						mat::SizesAndBlocks sizeBlockInfo,
						int n) {
  std::vector<ergo_real> fullMat(n*n);
  for (int col = 0; col < n; ++col)
    for (int row = 0; row < n; ++row) {
      ergo_real randomNumber = (rand() / (ergo_real)RAND_MAX);
      fullMat[row + col * n] = randomNumber;
      fullMat[col + row * n] = randomNumber;
    }
  symmMatrix A;
  A.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  A.assignFromFull(fullMat);
  symmMatrix A2;
  A2.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  Util::TimeMeter tm;
  A2 = (ergo_real)1.0 * A * A;
  double secondsTaken = tm.get_wall_seconds() - tm.get_start_time_wall_seconds();
  printf("A2 = 1.0 * A * A operation took         %6.3f wall seconds.\n", secondsTaken);
  result = A2;
}

int main(int argc, char *argv[])
{
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

  int n = 50;
  mat::SizesAndBlocks sizeBlockInfo;
  static const int sparseMatrixBlockSize = 1, sparseMatrixBlockFactor = 2;
  sizeBlockInfo =
    prepareMatrixSizesAndBlocks(n,
				sparseMatrixBlockSize,
				sparseMatrixBlockFactor,
				sparseMatrixBlockFactor,
				sparseMatrixBlockFactor);

  symmMatrix S;
  get_random_positive_definite_matrix(S, sizeBlockInfo, n);

  // Get inverse Cholesky factor Z
  triangMatrix Z;
  Z.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  Util::TimeMeter tm_inch;
  Z.inch(S, 0, mat::right);
  double secondsTaken_inch = tm_inch.get_wall_seconds() - tm_inch.get_start_time_wall_seconds();
  printf("inch operation took                     %6.3f wall seconds.\n", secondsTaken_inch);

  symmMatrix X;
  X.resetSizesAndBlocks(sizeBlockInfo, sizeBlockInfo);
  X = 1;

  Util::TimeMeter tm;
  X = transpose(Z) * X * Z;
  double secondsTaken = tm.get_wall_seconds() - tm.get_start_time_wall_seconds();
  printf("X = transpose(Z) * X * Z operation took %6.3f wall seconds.\n", secondsTaken);

  puts("mat_thread_test finished OK.");
  return 0;
}
