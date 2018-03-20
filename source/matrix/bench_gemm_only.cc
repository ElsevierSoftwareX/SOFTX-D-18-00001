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

/** @file bench_gemm_only.cc Benchmark of the matrix library with
 * input parameters specifying block sizes, parallel level etc.
 */

#include <fstream>  /* For ifstream */
#include <iomanip> /* For setprecision in fstream */
#include <iostream>
#include <cmath>
#include <stdio.h> /* For FILE */
#include <sys/time.h>
#include "Matrix.h"
#include "Vector.h"
#include "MatrixSymmetric.h"
#include "MatrixTriangular.h"
#include "MatrixGeneral.h"
#include "VectorGeneral.h"
#include "mat_gblas.h"
#include "Lanczos.h"

static double get_wall_seconds() {
  struct timeval tv;
  if(gettimeofday(&tv, NULL) != 0)
    throw std::runtime_error("Error in get_wall_seconds(), in gettimeofday().");
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

using namespace mat;

template<typename real>
int mainFun(int size, int parallel_level, int block_size, 
	    int block_factor_1, int block_factor_2, int block_factor_3) {  
  typedef Matrix<real, real> Mat_1;
  typedef Matrix<real, Mat_1> Mat_2;
  typedef Matrix<real, Mat_2> Mat_3;
  typedef Matrix<real, Mat_3> Mat_4;
  
  typedef Mat_4 matri;
  typedef MatrixGeneral<real, matri> normalMatrix;
  
  try {

    /********** Initialization of SizesAndBlocks                            */
    int nlevels = 4;
    std::vector<int> blockSizes(nlevels);
    blockSizes[3] = 1;
    blockSizes[2] = blockSizes[3] * block_size;
    blockSizes[1] = blockSizes[2] * block_factor_1;
    blockSizes[0] = blockSizes[1] * block_factor_2;
    for(int i = 0; i < nlevels; i++)
      std::cout << "blockSizes[" << i << "] = " << blockSizes[i] << std::endl;
    SizesAndBlocks rows(blockSizes, size);
    SizesAndBlocks cols(blockSizes, size);
    real alpha = 0.77;
    
    {
      normalMatrix A, B, C;
      A.resetSizesAndBlocks(rows,cols);
      B.resetSizesAndBlocks(rows,cols);
      A.random();
      B.random();

      double startTime = get_wall_seconds();
      C = alpha * A * B;
      double secondsTaken = get_wall_seconds() - startTime;
      std::cout << "Operation C = alpha * A * B took " << secondsTaken << " wall seconds." << std::endl;
    }
  }
  catch (Failure e) {
    std::cout << "Failure caught: "<<e.what() << std::endl;
    std::exit(1);
  }
  catch (std::exception e) {
    std::cout << "Exception caught: "<<e.what() << std::endl;
    std::exit(1);
  }
  return 0;
}


int main(int argc,char* argv[]){
  double startTime;

  int do_single = 1;
  int do_double = 1;
  int do_longdouble = 1;
  int size = 100;
  int parallel_level = 2;
  int block_size = 10;
  int block_factor_1 = 8;
  int block_factor_2 = 8;
  int block_factor_3 = 8;
  if(argc != 10) {
    std::cout << "argc != 10, using default parameters." << std::endl;
  }
  else {
    do_single = atoi(argv[1]);
    do_double = atoi(argv[2]);
    do_longdouble = atoi(argv[3]);
    size = atoi(argv[4]);
    parallel_level = atoi(argv[5]);
    block_size = atoi(argv[6]);
    block_factor_1 = atoi(argv[7]);
    block_factor_2 = atoi(argv[8]);
    block_factor_3 = atoi(argv[9]);    
  }
  std::cout << "Parameter values:" << std::endl;
  std::cout << "do_single = " << do_single << std::endl;
  std::cout << "do_double = " << do_double << std::endl;
  std::cout << "do_longdouble = " << do_longdouble << std::endl;
  std::cout << "size = " << size << std::endl;
  std::cout << "parallel_level = " << parallel_level << std::endl;
  std::cout << "block_size = " << block_size << std::endl;
  std::cout << "block_factor_1 = " << block_factor_1 << std::endl;
  std::cout << "block_factor_2 = " << block_factor_2 << std::endl;
  std::cout << "block_factor_3 = " << block_factor_3 << std::endl;
  std::cout << std::endl;

#ifdef _OPENMP
  int defThreads;
  const char *env = getenv("OMP_NUM_THREADS");
  if ( !(env && (defThreads=atoi(env)) > 0) ) {
    defThreads = 1;
  }
  mat::Params::setNProcs(defThreads);
  mat::Params::setMatrixParallelLevel(parallel_level);
  std::cout<<"OpenMP is used, number of threads set to "
	   <<mat::Params::getNProcs()<<". Matrix parallel level: "
	   <<mat::Params::getMatrixParallelLevel()<<"."<<std::endl;
#endif

  if(do_single == 1) {
    std::cout<<"Benchmark of matrix library with single precision:" <<std::endl;
    startTime = get_wall_seconds();
    if (!mainFun<float>(size, parallel_level, block_size, 
			block_factor_1, block_factor_2, block_factor_3))
      std::cout<<"Matrix library benchmark with single precision completed "
	"successfully.\n" 
	       <<"Wall time: "
	       << get_wall_seconds() - startTime << " seconds.\n\n";
  }

  if(do_double == 1) {
    std::cout<<"Benchmark of matrix library with double precision:" <<std::endl;
    startTime = get_wall_seconds();
    if (!mainFun<double>(size, parallel_level, block_size, 
			 block_factor_1, block_factor_2, block_factor_3))
      std::cout<<"Matrix library benchmark with double precision completed "
	"successfully.\n"
	       <<"Wall time: "
	       << get_wall_seconds() - startTime << " seconds.\n\n";
  }

  if(do_longdouble == 1) {
    std::cout<<"Benchmark of matrix library with long double precision:" <<std::endl;
    startTime = get_wall_seconds();
    if (!mainFun<long double>(size, parallel_level, block_size, 
			      block_factor_1, block_factor_2, block_factor_3))
      std::cout<<"Matrix library benchmark with long double precision completed "
	"successfully.\n" 
	       <<"Wall time: "
	       << get_wall_seconds() - startTime << " seconds.\n\n";
  }

  std::exit(0);
}
