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

/** @file bench.cc Benchmark of the matrix library
 *
 * Copyright(c) Emanuel Rubensson 2007
 *
 * @author Emanuel Rubensson  @a responsible @a author
 * @date September 2007
 *
 */

#include <fstream>  /* For ifstream */
#include <iomanip> /* For setprecision in fstream */
#include <iostream>
#include <cmath>
#include <stdio.h> /* For FILE */
#include "Matrix.h"
#include "Vector.h"
#include "MatrixSymmetric.h"
#include "MatrixTriangular.h"
#include "MatrixGeneral.h"
#include "VectorGeneral.h"
#include "mat_gblas.h"
#include "Lanczos.h"

using namespace mat;

template<typename real>
int mainFun(int argc,char* argv[]) {
  
  typedef Matrix<real, real> Mat_1;
  typedef Matrix<real, Mat_1> Mat_2;
  typedef Matrix<real, Mat_2> Mat_3;
  
  typedef Mat_3 matri;
  typedef MatrixSymmetric<real, matri> symmMatrix;
  typedef MatrixGeneral<real, matri> normalMatrix;
  
  try {
    int size = 100;
    switch (argc)
      {
      case 1:
	std::cout<<"No input, using matrix size "<<size<<"."<<std::endl;
	break;
      case 2:
	size = atoi(argv[1]);	
	break;
      default:
	std::cerr<<"Wrong number of input arguments"<<std::endl;
	std::exit(1);
      }

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

    /********** Initialization of SizesAndBlocks                            */
    int nlevels = 3;
    std::vector<int> blockSizes(nlevels);
    blockSizes[nlevels - 1] = 1;
    for (int ind = nlevels - 2; ind >= 0; ind--)
      blockSizes[ind] = blockSizes[ind + 1] * 10;
    SizesAndBlocks rows(blockSizes, size);
    SizesAndBlocks cols(blockSizes, size);
    real alpha = 0.77;
    
    {
      normalMatrix A, B, C;
      A.resetSizesAndBlocks(rows,cols);
      B.resetSizesAndBlocks(rows,cols);
      A.random();
      B.random();

      C = alpha * A * B;
    }
    {
      symmMatrix syA, syB;
      syA.resetSizesAndBlocks(rows,cols);
      syB.resetSizesAndBlocks(rows,cols);
      syA.random();
      syB = alpha * syA * syA;
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
  time_t tm; 
  
  std::cout<<"Benchmark of matrix library with single precision:" <<std::endl;
  time(&tm);
  if (!mainFun<float>(argc,argv))
    std::cout<<"Matrix library benchmark with single precision completed "
      "successfully.\n" 
	     <<"Wall time: "
	     <<((unsigned long)time(NULL))-tm<<" seconds.\n\n";

  std::cout<<"Benchmark of matrix library with double precision:" <<std::endl;
  time(&tm);
  if (!mainFun<double>(argc,argv)) {
    std::cout<<"Matrix library benchmark with double precision completed "
      "successfully.\n"
	     <<"Wall time: "
	     <<((unsigned long)time(NULL))-tm<<" seconds.\n\n";
  }

  std::cout<<"Benchmark of matrix library with long double precision:" <<std::endl;
  time(&tm);
  if (!mainFun<long double>(argc,argv))
    std::cout<<"Matrix library benchmark with long double precision completed "
      "successfully.\n" 
	     <<"Wall time: "
	     <<((unsigned long)time(NULL))-tm<<" seconds.\n\n";

  std::exit(0);
}
